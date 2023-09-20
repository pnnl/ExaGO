# DM, DMNetwork,PetscSection functions used in OPFLOW

OPFLOW makes heavy use of PETSc's functionality for managing unstructured grids. Primarily, it is built on top of the DMNetwork object in PETSc. Additionally, it uses some other objects (PetscSection and DM) from PETSc. This document lists the different functions used by OPFLOW divided by the pre-solve (a one-time set up) stage and solve (repeatedly called during the iterations) stage.

## Pre-solve stage

### DMNetwork functions

1. Create the DMNetwork object 
```
DMNetworkCreate(mpicomm,DM *networkdm)
```

2. Register components
```
DMNetworkRegisterComponent(DM networkdm,const char compname[],size_t compstructsize,PetscInt *compkey)
```

3. Set the number of vertices/edges
```
DMNetworkSetSizes(DM networkdm,PetscInt nsubnet,PetscInt nV[],PetscInt nE[],NCouplesubnet,PetscInt nCoupleEdges[]);
```

4. Set the network connectivity/topology
```
DMNetworkSetEdgeList(DM networkdm,PetscInt edgelist*,PetscInt *edge)
```

5. Create layout (bare graph)
```
DMNetworkLayoutSetUp(DM networkdm)
```

6. Get the edge range (iterator for edges)
```
DMNetworkGetEdgeRange(DM,PetscInt *eStart, PetscInt *eEnd)
```

7.  Add component at vertex/edge point p
```
DMNetworkAddComponent(DM networkdm,PetscInt point,PetscInt compkey,void* compstruct)
```

8. Get the vertex range (iterator for vertices)
```
DMNetworkGetVertexRange(DM networkdm,PetscInt *vStart,PetscInt *vEnd)
```

9. Add variables to vertex/edge point p
```
DMNetworkAddNumVariables(DM networkdm,PetscInt p,PetscInt nx)
```

10. Distribute network (partition and distribute components)
```
DMNetworkDistribute(DM *dm,PetscInt part)
```

11. Get the number of components at a given vertex/edge point p
```
DMNetworkGetNumComponents(DM networkdm,PetscInt p,PetscInt *numComponents)
```
12. Get the jth component at vertex/edge point p
```
DMNetworkGetComponent(DM networkdm,PetscInt p,PetscInt j,PetscInt *key,void *compstruct)
```

13. Get the vertices covering edge e
```
DMNetworkGetConnectedVertices(DM networkdm,PetscInt e,PetscInt *connnodes)
```

14. Get the offset (starting location) for the variables for vertex/edge point p in the local vector.
```
DMNetworkGetVariableOffset(DM networkdm,PetscInt p,PetscInt *startloc)
```

15. Get the global offset (starting location) for the variables for vertex/edge point p in the global vector.
```
DMNetworkGetGlobalOffset(DM networkdm,PetscInt p, PetscInt *startlocglob)
```

16. Get the edges incident at a vertex p
```
DMNetworkGetSupportingEdges(DM networkdm, PetscInt p, PetscInt *nlines, PetscInt *connlines)
```

17. Determine if a vertex is a ghost vertex
```
DMNetworkIsGhostVertex(DM networkdm, PetscBool *isghost)
```

### PetscSection functions
PetscSection is an object associated with the DM. It stores the starting locations for the variables for each vertex/edge point p. Every DM has two sections - a `localsection` for local indices, and a `global` section for global indices.

1. Create PetscSection
```
PetscSectionCreate(MPI_Comm comm,PetscSection *section)
```
2. Set the range for the section
``` 
PetscSectionSetChart(PetscSection section,PetscInt pstart,PetscInt pend)
```

3. Set the degrees of freedom (dof) at point p
```
PetscSectionSetDof(PetscSection section, PetscInt dof)
```
4. Set up section
```
PetscSectionSetUp(PetscSection section)
```

5. Destroy section
```
PetscSectionDestroy(PetscSection section)
```

### DM functions

1. Get the PetscSection associated with the DM
```
DMGetSection(DM dm, PetscSection *section)
```

2. Get the global section for the DM
```
DMGetGlobalSection(DM dm, PetscSection *globalsection)
```

3. Set a section on the DM
```
DMSetSection(DM dm,PetscSection section)
```

## Solve stage

1. Get/Restore the local vector
```
DMGetLocalVector(DM networkdm,Vec &localVec)
DMRestoreLocalVector(DM networkdm,Vec &localVec)
```

2. Scatter values from global vector to local vector
```
DMGlobalToLocalBegin(DM networkdm,Vec globalVec,INSERT_MODE,Vec localVec)
DMGlobalToLocalEnd(DM networkdm,Vec globalVec,INSERT_MODE,Vec localVec)
```

3. Scatter values from local vector to global vector
```
DMLocalToGlobalBegin(DM networkdm,Vec localVec,INSERT_MODE,Vec globalVec)
DMLocalToGlobalEnd(DM networkdm,Vec localVec,INSERT_MODE,Vec globalVec)
```
