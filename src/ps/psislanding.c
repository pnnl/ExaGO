#include <private/psimpl.h>

/*
  PSConnCompDestroy - Destroys the connected components struct
*/
PetscErrorCode PSConnCompDestroy(PS ps)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFunctionBegin;
  for(i=0; i < ps->nconncomp; i++) {
    ps->conncomp[i].nv = 0;
    ierr = PetscFree(ps->conncomp[i].v);CHKERRQ(ierr);
  }
  ps->nconncomp = 0;
  PetscFunctionReturn(0);
}

/*
  PSIslandCheckandSetRefBus - Checks for active ref. bus and sets one if not available. This
    is a version of PSCheckandSetRefBus that allows setting reference buses on islands.

  Input Parameters:
+ ps - The PS object
- isnum  - island number

  This routine checks if there is an active ref. bus (declared in the data file and has active generators). If it is
  not defined, then the first PV bus is used.
*/
PetscErrorCode PSIslandCheckandSetRefBus(PS ps,PetscInt isnum)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PSBUS          bus;
  PetscBool      pshaslocalref=PETSC_FALSE, pshasref;
  PetscInt       firstlocpvbus=100000000,firstpvbus;

  PetscFunctionBegin;
  if(ps->conncomp[isnum].blackout) PetscFunctionReturn(0);

  for(i=0; i < ps->conncomp[isnum].nv; i++) {
    bus = &ps->bus[ps->busext2intmap[ps->conncomp[isnum].v[i]]];
    if(bus->ide == REF_BUS && bus->ngenON) {
      /* PS has active reference bus */
      pshaslocalref = PETSC_TRUE;
      break;
    }
  }

  ierr = MPI_Allreduce(&pshaslocalref,&pshasref,1,MPIU_BOOL,MPI_LOR,ps->comm->type);CHKERRQ(ierr);

  if(pshasref) PetscFunctionReturn(0);
#if defined DEBUGPS
  ierr = PetscPrintf(PETSC_COMM_SELF,"No active ref. bus in island %d\n",isnum+1);CHKERRQ(ierr);
#endif
  /* No active ref. bus, set first PV bus to ref bus */
  for(i=0; i < ps->conncomp[isnum].nv; i++) {
    bus = &ps->bus[ps->busext2intmap[ps->conncomp[isnum].v[i]]];
    if(bus->ide == PV_BUS && bus->ngenON && !bus->isghost) {
      firstlocpvbus = bus->bus_i;
      break;
    }
  }

  ierr = MPIU_Allreduce(&firstlocpvbus,&firstpvbus,1,MPIU_INT,MPI_MIN,ps->comm->type);CHKERRQ(ierr);

  if(firstpvbus == 100000000) {
#if defined DEBUGPS
      SETERRQ1(PETSC_COMM_SELF,0," ",firstpvbus);
#endif
  } else {
    if (ps->busext2intmap[firstpvbus] != -1) {
      ps->bus[ps->busext2intmap[firstpvbus]].ide = REF_BUS;
      ps->nref++;
      ierr = PetscPrintf(PETSC_COMM_SELF,"Setting bus %d as the new ref. bus\n ",firstpvbus);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/* 
   PSGetConnectedComponents - Does a depth-first search to find the connected components starting from a given bus

   Input Parameters
.  start_bus - The starting bus/node

   Output Parameters
+  visited_nodes -- an array of visited nodes (if a node i is visited then visited_nodes[i] != 0, else it is 0)
-  nnodes_visited -- number of nodes visited (connected components)
*/
PetscErrorCode PSGetConnectedComponents(PSBUS start_bus, PetscInt *visited_nodes, PetscInt *nnodes_visited,PetscInt nrem)
{
  PetscErrorCode ierr;
  const PSLINE   *connlines;
  PSLINE         line;
  PetscInt       nconnlines;
  PSBUS          busf,bust,next_bus;
  PetscInt       i;
  const PSBUS    *connbuses;

  PetscFunctionBegin;

  if(visited_nodes[start_bus->internal_i] == 0) {
    visited_nodes[start_bus->internal_i] = 1;
    *nnodes_visited += 1;
  }
  if(*nnodes_visited == nrem) PetscFunctionReturn(0);

  ierr = PSBUSGetSupportingLines(start_bus,&nconnlines,&connlines);CHKERRQ(ierr);

  for(i=0; i < nconnlines; i++) {
    line = connlines[i];
    if(!line->status) continue;

    ierr = PSLINEGetConnectedBuses(line,&connbuses);CHKERRQ(ierr);
    busf = connbuses[0];
    bust = connbuses[1];

    if(start_bus == busf) next_bus = bust;
    else next_bus = busf;

    if(visited_nodes[next_bus->internal_i] == 0) {
      ierr =  PSGetConnectedComponents(next_bus, visited_nodes, nnodes_visited,nrem);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}

/*
  PSCheckTopology - Checks the network topology and finds if islands exist

*/
PetscErrorCode PSCheckTopology(PS ps)
{
  PetscErrorCode ierr;
  PetscInt       *visited_nodes,nnodes_visited=0,nrem=ps->nbus,i,j;
  PSBUS          start_bus=0;
  PetscInt       Nconncomp=0;
  PetscInt       *scat_conncomp_sizes;
  PetscInt       *scatv_displs,*scatv_sendcounts;
  PetscInt       *scatv_sendbuf;
  PSConnCompgroup *cgroup;
  PetscInt *ncg_packet_glob,pctr=0,ii;
  PSConngroup grp;

  PetscFunctionBegin;
  //  if(ps->comm->size > 1) SETERRQ(PETSC_COMM_SELF,0,"No parallel support for finding islands yet\n");

  if(ps->nconncomp) {
    ierr = PSConnCompDestroy(ps);CHKERRQ(ierr);
  }

  ierr = PetscCalloc1(ps->nbus,&visited_nodes);CHKERRQ(ierr);

  while(nrem != 0) {
    nnodes_visited = 0;
    for(i=0;i < ps->nbus; i++) {
      if(visited_nodes[i] == 0) {
	start_bus = &ps->bus[i];
	break;
      }
    }
    ierr = PSGetConnectedComponents(start_bus,visited_nodes,&nnodes_visited,nrem);CHKERRQ(ierr);

    nrem -= nnodes_visited;
    ps->conncomp[ps->nconncomp].nv = 0;
    ierr = PetscCalloc1(nnodes_visited,&ps->conncomp[ps->nconncomp].v);CHKERRQ(ierr);
    for(i=0;i < ps->nbus; i++) {
      if(visited_nodes[i] == 1) {
	ps->conncomp[ps->nconncomp].v[ps->conncomp[ps->nconncomp].nv++] = ps->bus[i].bus_i;
	visited_nodes[i]++; /* Advance so that the node does not get included in the next iteration */
      }
    }
    ps->nconncomp++;
  }

  ierr = PetscFree(visited_nodes);CHKERRQ(ierr);

  /* At this stage, all processes know about the connected components for their subnetwork */
  /* Prepare to send size information to processor 0 */
  PetscInt ncg_packet_size=0; /* Size of the packet containing the size information */
  PetscInt *ncg_packet;
  PetscInt *tp;

  /* Each packet has the following layout */
  /* [Rank | # of connected componets | size of each connected component | vertices in each connected component] */

  for(i=0; i < ps->nconncomp; i++) {
    ncg_packet_size += ps->conncomp[i].nv + 1; /* Add 1 for its size info */
  }
  ncg_packet_size += 2; /* Add rank and number of connected componenets */

  ierr = PetscCalloc1(ncg_packet_size,&ncg_packet);CHKERRQ(ierr);
  tp = ncg_packet + 2 + ps->nconncomp;
  ncg_packet[0] = ps->comm->rank;
  ncg_packet[1] = ps->nconncomp;
  for(i=0; i < ps->nconncomp; i++) {
    ncg_packet[2+i] = ps->conncomp[i].nv;
    for(j=0; j < ps->conncomp[i].nv; j++) tp[j] = ps->conncomp[i].v[j];
    tp += ps->conncomp[i].nv;
  }

  /* Uncomment for debugging
  for(i=0; i < ncg_packet_size; i++) {
    if(i == ncg_packet_size - 1) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d   \n",ncg_packet[i]);CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d   ",ncg_packet[i]);CHKERRQ(ierr);
    }
  }
  */
  /* Send packet size and number of connected components on each process to root */
  PetscInt sendbuf[2],*recvbuf;
  PetscInt *displs,*reccnt;
  if(!ps->comm->rank) {
    ierr = PetscCalloc1(2*ps->comm->size,&recvbuf);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->comm->size,&displs);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->comm->size,&reccnt);CHKERRQ(ierr);
  }

  sendbuf[0] = ncg_packet_size; /* Size of the packet sent for reduction */
  sendbuf[1] = ps->nconncomp; /* Number of connected components on each process */
  ierr = MPI_Gather(sendbuf,2,MPIU_INT,recvbuf,2,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);
  
  /* Gather packets from each process onto root */
  PetscInt ncg_packet_size_glob=0,nconncomp_glob=0;

  if(!ps->comm->rank) {
    for(i=0; i < ps->comm->size; i++) {
      displs[i] = ncg_packet_size_glob;
      ncg_packet_size_glob += recvbuf[2*i];
      nconncomp_glob += recvbuf[2*i+1];
      reccnt[i] = recvbuf[2*i];
    }
    /* Uncomment for debugging
    ierr = PetscPrintf(PETSC_COMM_SELF,"global packet size = %d, total number of connected components = %d\n",ncg_packet_size_glob,nconncomp_glob);CHKERRQ(ierr);
    */
    ierr = PetscCalloc1(ncg_packet_size_glob,&ncg_packet_glob);CHKERRQ(ierr);
  }
  
  ierr = MPI_Gatherv(ncg_packet,ncg_packet_size,MPIU_INT,ncg_packet_glob,reccnt,displs,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);

  if(!ps->comm->rank) {
    ierr = PetscFree(displs);CHKERRQ(ierr);
    ierr = PetscFree(reccnt);CHKERRQ(ierr);

    /* Uncomment for debugging
    for(i=0; i < ncg_packet_size_glob; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"%d   ",ncg_packet_glob[i]);CHKERRQ(ierr);
    }
    */


    /* Create conncompgroup structs. The number of elements equals to the
       total number of connected components. Each struct will have the rank,
       size of the connected component, and the connected component info. The
       idea is to think of each element of this connected component group struct
       as a vertex on a graph. A vertex is connected to another vertex if they have
       common nodes in the 'c' field of the struct.
    */

    PetscInt *ptr=ncg_packet_glob;
    PetscInt kk,nconncomp_proc,*data_ptr,rank;
    ierr = PetscCalloc1(nconncomp_glob,&cgroup);CHKERRQ(ierr);
    for(i=0; i < ps->comm->size; i++) {
      rank = *ptr;
      nconncomp_proc = *(ptr+1);
      ptr += 2;
      data_ptr = ptr+nconncomp_proc;
      for(kk=0;kk < nconncomp_proc; kk++) { 
	cgroup[pctr].rank = rank;
	cgroup[pctr].nc = *ptr++;
	ierr = PetscCalloc1(cgroup[pctr].nc,&cgroup[pctr].c);CHKERRQ(ierr);
	ierr = PetscMemcpy(cgroup[pctr].c,data_ptr,cgroup[pctr].nc*sizeof(PetscInt));CHKERRQ(ierr);
	data_ptr += cgroup[pctr].nc;
	pctr++;
      }
      ptr = data_ptr;
    }

    /*
    for(i=0; i < pctr; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"ctr = %d rank = %d\n",i,cgroup[i].rank);CHKERRQ(ierr);
      for(kk=0; kk < cgroup[i].nc; kk++) {
	ierr = PetscPrintf(PETSC_COMM_SELF,"%d  ",cgroup[i].c[kk]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    }
    */

    PSConngroupi *gi;
    PetscInt *visited_nodegroups;
    PetscInt nrem=nconncomp_glob,k,nvis=0;
    PSConnCompgroup *start_group;

    nrem = nconncomp_glob;
    grp.n = 0;
    ierr = PetscCalloc1(nconncomp_glob,&grp.ci);CHKERRQ(ierr);
    ierr = PetscCalloc1(nconncomp_glob,&visited_nodegroups);CHKERRQ(ierr);

    while(nrem) {
      gi = &grp.ci[grp.n];
      gi->nc = 0;
      nvis = 0;
      ierr = PetscCalloc1(nconncomp_glob,&gi->cg);CHKERRQ(ierr);
      for(i=0; i < nconncomp_glob; i++) {
	if(visited_nodegroups[i] == 0) {
	  visited_nodegroups[i] += 1;
	  start_group = &cgroup[i];
	  gi->cg[gi->nc++] = start_group;
	  nvis++;
	  break;
	}
      }
      
      PetscInt new_node=1;
      PetscInt l;
      PSConnCompgroup *cgi;
      
      while(new_node) {
	new_node = 0;
	for(i=0; i < nconncomp_glob; i++) {
	  if(visited_nodegroups[i]) continue;
	  for(j=0; j < cgroup[i].nc; j++) {
	    //	    for(k=0; k < start_group->nc; k++) {
	    for(k=0; k < gi->nc; k++) {
	      cgi = gi->cg[k];
	      for(l=0; l < cgi->nc; l++) {
	      //	    ierr = PetscPrintf(PETSC_COMM_SELF,"start_group->c[%d] = %d  cgroup[%d].c[%d] = %d\n",k,start_group->c[k],i,j,cgroup[i].c[j]);CHKERRQ(ierr);
		if(cgi->c[l] == cgroup[i].c[j]) {
		//	      ierr = PetscPrintf(PETSC_COMM_SELF,"Came here\n");CHKERRQ(ierr);
		  visited_nodegroups[i] += 1;
		  gi->cg[gi->nc++] = &cgroup[i];
		  nvis++;
		  /* Break from the loop */
		  l = cgi->nc;
		  j = cgroup[i].nc;
		  k = gi->nc;
		  new_node = 1;
		}
	      }
	    }
	  }
	}
      }
      nrem -= nvis;
      grp.n++;
    }



    /* Uncomment for debugging
    for(i=0; i < grp.n; i++) {
      gi = &grp.ci[i];
      ierr = PetscPrintf(PETSC_COMM_SELF,"Group number = %d\n",i);CHKERRQ(ierr);
      for(j=0; j < gi->nc; j++) {
	for(k=0; k < gi->cg[j]->nc; k++) {
	  ierr = PetscPrintf(PETSC_COMM_SELF,"  %d  ",gi->cg[j]->c[k]);CHKERRQ(ierr);
	}
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    }
    */

    /* Create buffer for scattering sizes */
    ierr = PetscCalloc1(ps->comm->size*(grp.n+1),&scat_conncomp_sizes);CHKERRQ(ierr);
    for(i=0; i < ps->comm->size; i++) {
      scat_conncomp_sizes[(grp.n+1)*i] = grp.n;
    }
    for(i=0; i < grp.n; i++) {
      gi = &grp.ci[i];
      for(j=0; j < gi->nc; j++) {
       	scat_conncomp_sizes[(gi->cg[j]->rank)*(grp.n+1)+1+i] += gi->cg[j]->nc;
      }
    }

    ierr = PetscCalloc1(ps->comm->size,&scatv_displs);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->comm->size,&scatv_sendcounts);CHKERRQ(ierr);
    scatv_displs[0] = 0;
    PetscInt scatv_sendbuf_size = 0;
    for(i=0; i < ps->comm->size; i++) {
      if(i >= 1) scatv_displs[i] = scatv_displs[i-1] + scatv_sendcounts[i-1];
      for(j=0; j < grp.n+1; j++) {
	if(j >= 1) scatv_sendcounts[i] += scat_conncomp_sizes[(grp.n+1)*i+j];
	//ierr = PetscPrintf(PETSC_COMM_SELF,"%d  ",scat_conncomp_sizes[(grp.n+1)*i+j]);CHKERRQ(ierr);
      }
      //ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
      scatv_sendbuf_size += scatv_sendcounts[i];
    }
    //ierr = PetscPrintf(PETSC_COMM_SELF,"scatv_sendbuf_size = %d\n",scatv_sendbuf_size);CHKERRQ(ierr);
    
    ierr = PetscCalloc1(scatv_sendbuf_size,&scatv_sendbuf);CHKERRQ(ierr);
    PetscInt **displs_ptr;
    ierr = PetscCalloc1(ps->comm->size,&displs_ptr);CHKERRQ(ierr);
    for(i=0; i < ps->comm->size;i++) {
      /* Initialize displacement pointers to mark to the beginning of the array for each process */
      displs_ptr[i] = scatv_sendbuf + scatv_displs[i];
    }
    for(i=0; i < grp.n; i++) {
      gi = &grp.ci[i];
      for(j=0; j < gi->nc; j++) {
	ierr = PetscMemcpy(displs_ptr[gi->cg[j]->rank],gi->cg[j]->c,gi->cg[j]->nc*sizeof(PetscInt));CHKERRQ(ierr);
	displs_ptr[gi->cg[j]->rank] += gi->cg[j]->nc;
      }
    } 
    ierr = PetscFree(displs_ptr);CHKERRQ(ierr);

    /* Uncomment for debugging
    for(i=0; i < scatv_sendbuf_size; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"  %d  ",scatv_sendbuf[i]);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    */
    /*
    for(i=0; i < nconncomp_glob; i++) {
      ierr = PetscPrintf(PETSC_COMM_SELF,"Rank = %d, Ncomp = %d  ",cgroup[i].rank,cgroup[i].nc);CHKERRQ(ierr);
      for(kk=0; kk < cgroup[i].nc;kk++) {
	ierr = PetscPrintf(PETSC_COMM_SELF,"%d  ",cgroup[i].c[kk]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\n");CHKERRQ(ierr);
    }
    */
    /*
    ncg_packet[0] = ps->comm->rank;
    ncg_packet[1] = ps->nconncomp;
    PetscInt *tp=ncg_packet + 2 + ps->nconncomp;
    for(i=0; i < ps->nconncomp; i++) {
      ncg_packet[2+i] = ps->conncomp[i].nv;
      for(j=0; j < ps->conncomp[i].nv; j++) tp[j] = ps->conncomp[i].v[j];
      tp += ps->conncomp[i].nv;
    }
    */

    Nconncomp = grp.n;

    ierr = PetscFree(visited_nodegroups);CHKERRQ(ierr);
  }
  
  MPI_Barrier(ps->comm->type);

  ierr = PSConnCompDestroy(ps);CHKERRQ(ierr);
  ps->nconncomp = Nconncomp;
  ierr = MPI_Bcast(&ps->nconncomp,1,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] nconncomp=%d\n",ps->comm->rank,ps->nconncomp);CHKERRQ(ierr);
  
  PetscInt *conncomp_sizes;
  ierr = PetscCalloc1(ps->nconncomp+1,&conncomp_sizes);CHKERRQ(ierr);
  ierr = MPI_Scatter(scat_conncomp_sizes,ps->nconncomp+1,MPIU_INT,conncomp_sizes,ps->nconncomp+1,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);
  for(i=0; i < ps->nconncomp; i++) {
    ps->conncomp[i].nv = conncomp_sizes[i+1];
    //ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] nv[%d] = %d\n",ps->comm->rank,i,ps->conncomp[i].nv);CHKERRQ(ierr);
    ierr = PetscCalloc1(ps->conncomp[i].nv,&ps->conncomp[i].v);CHKERRQ(ierr);    
  }  
  ierr = PetscFree(conncomp_sizes);CHKERRQ(ierr);

  PetscInt *scatv_recvbuf;
  ierr = PetscCalloc1(ps->nbus,&scatv_recvbuf);CHKERRQ(ierr);

  ierr = MPI_Scatterv(scatv_sendbuf,scatv_sendcounts,scatv_displs,MPIU_INT,scatv_recvbuf,ps->nbus,MPIU_INT,0,ps->comm->type);CHKERRQ(ierr);

  PetscInt *ptr=scatv_recvbuf;
  for(i=0; i < ps->nconncomp; i++) {
    ierr = PetscMemcpy(ps->conncomp[i].v,ptr,ps->conncomp[i].nv*sizeof(PetscInt));CHKERRQ(ierr);
    for(j=0; j < ps->conncomp[i].nv; j++) {
      //ierr = PetscPrintf(PETSC_COMM_SELF,"Rank[%d] Group[%d]  Node = %d\n",ps->comm->rank,i,ps->conncomp[i].v[j]);CHKERRQ(ierr);
    }
    ptr += ps->conncomp[i].nv;
  }

  if(!ps->comm->rank) {
    ierr = PetscFree(scatv_displs);CHKERRQ(ierr);
    ierr = PetscFree(scatv_sendcounts);CHKERRQ(ierr);
    ierr = PetscFree(scatv_sendbuf);CHKERRQ(ierr);
    ierr = PetscFree(scat_conncomp_sizes);CHKERRQ(ierr);
    ierr = PetscFree(ncg_packet_glob);CHKERRQ(ierr);
    ierr = PetscFree(recvbuf);CHKERRQ(ierr);

    for(ii=0; ii < pctr; ii++) {
      ierr = PetscFree(cgroup[ii].c);CHKERRQ(ierr);
    }
    ierr = PetscFree(cgroup);CHKERRQ(ierr);

    for(ii=0; ii < grp.n; ii++) {
      ierr = PetscFree(grp.ci[ii].cg);CHKERRQ(ierr);
    }
    ierr = PetscFree(grp.ci);CHKERRQ(ierr);
  }
  ierr = PetscFree(ncg_packet);CHKERRQ(ierr);
  ierr = PetscFree(scatv_recvbuf);CHKERRQ(ierr);

    
  /* Check if the island has an active generator */
  for(i=0; i < ps->nconncomp; i++) {
    PetscInt islandhasgen = 0,k,islandhasgen_glob;
    PSBUS    bus;
    for(k=0; k < ps->conncomp[i].nv; k++) {
      bus = &ps->bus[ps->busext2intmap[ps->conncomp[i].v[k]]];
      if(bus->ngenON) {
	islandhasgen = 1;
	ps->conncomp[i].blackout = 0;
	break;
      }
    }
    
    ierr = MPI_Allreduce(&islandhasgen,&islandhasgen_glob,1,MPIU_INT,MPI_SUM,ps->comm->type);CHKERRQ(ierr);

    if(!islandhasgen_glob) { /* No generators on this island, set all buses isolated on this island */
      PSLOAD load;
      PSLINE line;
      PetscInt nconnlines;
      const PSLINE *connlines;
      PetscInt n;

      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Island %d black out\n",i+1);CHKERRQ(ierr);
      ps->conncomp[i].blackout = 1;
      for(k=0; k < ps->conncomp[i].nv; k++) {
	bus = &ps->bus[ps->busext2intmap[ps->conncomp[i].v[k]]];
	bus->ide = ISOLATED_BUS;
	/* Turn OFF connected lines */
	ierr = PSBUSGetSupportingLines(bus,&nconnlines,&connlines);CHKERRQ(ierr);
	for(n=0; n < nconnlines; n++) {
	  line = connlines[n];
	  line->status = 0;
	}
	
	/* Turn OFF connected loads */
	for(n=0; n < bus->nload; n++) {
	  ierr = PSBUSGetLoad(bus,n,&load);CHKERRQ(ierr);
	  load->status = 0;
	}
      }
    }
  }

  /*
  if(ps->nconncomp > 1) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Found %d islands in the network\n",ps->nconncomp);CHKERRQ(ierr);
    for(i=0; i < ps->nconncomp; i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Buses in island %d:\n",i+1);CHKERRQ(ierr);
      for(j=0; j < ps->conncomp[i].nv; j++) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Bus %d\n",ps->conncomp[i].v[j]);CHKERRQ(ierr);
      }
    }
  }
  */

  for(i=0; i < ps->nconncomp; i++) {
    ierr = PSIslandCheckandSetRefBus(ps,i);CHKERRQ(ierr);
  }

  ierr = PetscFree(visited_nodes);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

