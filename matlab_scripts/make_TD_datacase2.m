%% make_TD_case
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

mpopt = mpoption('pf.nr.max_it',50);
mpopt = mpoption(mpopt,'cpf.stop_at','NOSE','cpf.plot.level',2,'cpf.enforce_v_lims',1)
t_case = 'datafiles/case_ACTIVSg200';
%% Test case 1
%d_case = 'case22';
%% Test case 2
d_case = 'case141';

%% Load t_case
mpct = runpf(t_case);
%% Load d_case
mpcdb = runpf(d_case);

t_bdry_buses = [];
d_inj = mpcdb.gen(1,PG);

load_buses = find(mpct.bus(:,PD)~=0 | mpct.bus(:,QD)~=0);
ndnets = ceil(mpct.bus(load_buses,PD)./d_inj);
ii = find(ndnets == 0);
ndnets(ii) = 1;

% Pick the first nn load buses
nn=20;
%load_buses = load_buses(1:nn);
load_buses = [2;6;62;177]; % load buses for IEEE T&D presentation

ctr = 0;
for i=1:length(load_buses)
    ndnets = ceil(mpct.bus(load_buses(i),PD)/d_inj);
    if(~ndnets)
        ndnets = 1;
    end
 %   if(ndnets == 1)
 %       continue;
 %   end
    for k=1:ndnets
        mpcdbs{ctr+1}.datas{k} = mpcdb;
    end
    mpcdbs{ctr+1}.nbdry = ndnets;
    mpcdbs{ctr+1}.bdry_bus = mpct.bus(load_buses(i),1);
%   t_bdry_buses = [t_bdry_buses;[mpct.bus(load_buses(i),1),ndnets]];

    ctr = ctr + 1;
%     if(ctr >= 2)
%         break;
%     end
end
nbdry_buses = length(mpcdbs);

for i=1:nbdry_buses
    for k = 1:mpcdbs{i}.nbdry
        mpcdb = mpcdbs{i}.datas{k};
        %% Set vmin = 0.95 and vmax = 1.05. These will be used in CPF
        %% when enforcing limits
        mpcdb.bus(:,VMIN) = 0.95;
        mpcdb.bus(:,VMAX) = 1.05;

        %% Calculate mpcd.branch resistance and reactance on the mpct KV and MVAbase.
        ref = find(mpcdb.bus(:,BUS_TYPE) == REF);

        %% Make the baseKV for the boundary buses same
        %mpcdb.bus(ref,BASE_KV) = mpct.bus(t_bdry_bus,BASE_KV);

        %% Calculate base impedance for lines at the dist. base MVA
        Zbaseline_dbase = (mpcdb.bus(mpcdb.branch(:,F_BUS),BASE_KV)*1e3).^2/(mpcdb.baseMVA*1e6);

        %% Calculate actual impedance
        Zmpcd_ohm = (mpcdb.branch(:,BR_R) + sqrt(-1)*mpcdb.branch(:,BR_X)).*Zbaseline_dbase;

        %% Calculate base impedance for lines at the trans. base MVA
        Zbaseline_tbase = (mpcdb.bus(mpcdb.branch(:,F_BUS),BASE_KV)*1e3).^2/(mpct.baseMVA*1e6);

        %% Calculate pu impedance at the trans. base MVA
        Zmpcdt_pu = Zmpcd_ohm./Zbaseline_tbase;

        %Zbased = (mpcdb.bus(ref,BASE_KV)*1e3)^2/(mpcdb.baseMVA*1e6);
        %Zmpcd_ohm = (mpcdb.branch(:,BR_R) + sqrt(-1)*mpcdb.branch(:,BR_X))*Zbased;
        %Zbaset = (mpct.bus(t_bdry_bus,BASE_KV)*1e3)^2/(mpct.baseMVA*1e6);
        %Zmpcdt_pu = Zmpcd_ohm/Zbaset;

        %% Put the converted values back in mpcd;
        mpcdb.branch(:,BR_R) = real(Zmpcdt_pu);
        mpcdb.branch(:,BR_X) = imag(Zmpcdt_pu);

        mpcdb.baseMVA = mpct.baseMVA;

        mpcdbs{i}.datas{k} = mpcdb;
    end
end

converged = 0;
tol = 1e-8;
ctr = 0;
max_it = 20;

err = zeros(2*nbdry_buses,1);
for i = 1:nbdry_buses
    tidx = find(mpct.bus(:,BUS_I) == mpcdbs{i}.bdry_bus);
    dinj = zeros(2,1);
    for k = 1:mpcdbs{i}.nbdry
        mpcdb = mpcdbs{i}.datas{k};
        dinj = dinj + [mpcdb.gen(1,PG);mpcdb.gen(1,QG)];
    end
    err(2*i-1:2*i,1) = [mpct.bus(tidx,PD)-dinj(1);mpct.bus(tidx,QD)-dinj(2)];
end

if(norm(err) < tol)
    converged = 1;
end
   
while(~converged & ctr < max_it)
    ctr = ctr + 1;
    %% Run T power flow
    mpct = runpf(mpct,mpopt);

    for i = 1:nbdry_buses
        tidx = find(mpct.bus(:,BUS_I) == mpcdbs{i}.bdry_bus);
        for k = 1:mpcdbs{i}.nbdry
            mpcdb = mpcdbs{i}.datas{k}
            %% Set reference bus voltage = t_boundary bus voltage
            mpcdb.bus(1,VM) = mpct.bus(tidx,VM);
            mpcdb.bus(1,VA) = mpct.bus(tidx,VA);
            mpcdb.gen(1,VG) = mpct.bus(tidx,VM);

            %% Run D power flow
            mpcdbs{i}.datas{k} = runpf(mpcdb,mpopt);
        end
    end
            
    err = zeros(2*nbdry_buses,1);
    for i = 1:nbdry_buses
        tidx = find(mpct.bus(:,BUS_I) == mpcdbs{i}.bdry_bus);
        dinj = zeros(2,1);
        for k = 1:mpcdbs{i}.nbdry
            mpcdb = mpcdbs{i}.datas{k};
            dinj = dinj + [mpcdb.gen(1,PG);mpcdb.gen(1,QG)];
        end
        err(2*i-1:2*i,1) = [mpct.bus(tidx,PD)-dinj(1);mpct.bus(tidx,QD)-dinj(2)];
    end
    
    if(norm(err) < tol)
        converged = 1;
    else
        for i=1:nbdry_buses
            tidx = find(mpct.bus(:,BUS_I) == mpcdbs{i}.bdry_bus);
            dinj = zeros(2,1);
            for k = 1:mpcdbs{i}.nbdry
                mpcdb = mpcdbs{i}.datas{k};
                dinj = dinj + [mpcdb.gen(1,PG);mpcdb.gen(1,QG)];
            end
            mpct.bus(tidx,PD) = dinj(1);
            mpct.bus(tidx,QD) = dinj(2);
        end
    end
end
if(~converged)
    error('Iterative T-D did not converge\n');
end
    
%mpcdb.bus(:,VMIN) = 0.98;
%mpcdb.bus(:,VMAX) = 1.02;
%% Check which bus voltages are below 0.95
% busi = find(mpcdb.bus(:,VM) < 0.95);
% while ~isempty(busi)
%     mpcdb.bus(busi,BS) = mpcdb.bus(busi,BS)+0.1;
%     mpcdb = runpf(mpcdb,mpopt);
%     %% Check which bus voltages are below 0.95
%     busi = find(mpcdb.bus(:,VM) < 0.98);
% end

% mpcdt = mpcdb;
% mpcdt.bus(:,PD) = mpcdt.bus(:,PD)*2.0;
% mpcdt.bus(:,QD) = mpcdt.bus(:,QD)*2.0;
% % 
% mpcd = runcpf(mpcdb,mpcdt,mpopt);

mpcd = mpcdb;

%% Combine the T&D cases
[nbust,colt] = size(mpct.bus);
nbus = nbust;

mpc = struct();
mpc.baseMVA = mpct.baseMVA;
mpc.version = mpct.version;
mpc.bus = mpct.bus;
mpc.gen = mpct.gen;
mpc.branch = mpct.branch;

for i = 1:nbdry_buses
    t_bdry_bus = find(mpct.bus(:,BUS_I) == mpcdbs{i}.bdry_bus);
    
    for k = 1:mpcdbs{i}.nbdry
        [nbust,colt] = size(mpc.bus);
        [ngent,colgt] = size(mpc.gen);
        [nbrancht,colbt] = size(mpc.branch);
        
        mpcd = mpcdbs{i}.datas{k};
        [nbusd,cold] = size(mpcd.bus);
        nbus = nbust+nbusd;

        maxbus = max(mpc.bus(:,1));
        dbusmap = [(1:nbusd)',(maxbus+1:maxbus+nbusd)',[t_bdry_bus,nbust+1:nbus-1]'];
        
        mpc.bus = [mpc.bus(:,1:min(colt,cold));mpcd.bus(:,1:min(cold,colt))];

        %% Renumber distribution bus
        mpc.bus(nbust+1:nbus,1) = dbusmap(mpcd.bus(:,1),2);
        
        [ngend,cold] = size(mpcd.gen);
        ngen = ngent+ngend;
        
        mpc.gen = [mpc.gen(:,1:min(colgt,cold));mpcd.gen(1:min(cold,colgt))];
        mpc.gen(ngent+1:ngen,1) = dbusmap(mpcd.gen(:,1),2);
        
        %% Find reference generator for the distribution and take it out
        refidx = find(mpc.bus(:,2) == REF);
        for kk=1:length(refidx)
            if mpc.bus(refidx(kk),1) > maxbus
                refd = refidx(kk);
                break;
            end
        end
        
        %% Convert REF Bus to PQ bus
        mpc.bus(refd,2) = PQ;
        %% Find reference generator for D and take it out.
        refgend = find(mpc.gen(:,1) == mpc.bus(refd,1));
        mpc.gen(refgend,:) = [];
        
        %% Branches
        [nbranchd,cold] = size(mpcd.branch);
        nbranch = nbrancht + nbranchd;
        mpc.branch = [mpc.branch(:,1:min(colbt,cold));mpcd.branch(:,1:min(cold,colbt))];

        mpc.branch(nbrancht+1:nbranch,F_BUS) = dbusmap(mpcd.branch(:,F_BUS),2);
        mpc.branch(nbrancht+1:nbranch,T_BUS) = dbusmap(mpcd.branch(:,T_BUS),2);

        %% Replace branches with distribution reference bus to t_bdry_bus
        mpc.branch(find(mpc.branch(:,F_BUS) == mpc.bus(refd,1)),F_BUS) = t_bdry_bus;
        mpc.branch(find(mpc.branch(:,T_BUS) == mpc.bus(refd,1)),T_BUS) = t_bdry_bus;
        
        %% Remove reference bus from the distribution part
        mpc.bus(refd,:) = [];
    end
    %% Remove load from t_bdry_bus
    mpc.bus(t_bdry_bus,PD:QD) = 0;
end
    
mpopt.verbose = 2;
mpcout = runpf(mpc,mpopt);
tcasename = strsplit(t_case,'/');
if(length(tcasename) > 1)
    tcasename = tcasename{2};
else
    tcasename = tcasename{1};
end
dcasename = strsplit(d_case,'/');
if(length(dcasename) > 1)
    dcasename = dcasename{2};
else
    dcasename = dcasename{1};
end

nfeeders = 0;
for i = 1:nbdry_buses
    nfeeders = nfeeders + mpcdbs{i}.nbdry;
end

dirname = ['datafiles/dir_',tcasename,'_',dcasename,'_nbdry_',num2str(nbdry_buses),'_nfeeders_',num2str(nfeeders)];
if exist(dirname)
    rmdir(dirname,'s');
end
mkdir(dirname);

%% This is the file for the integrated T-D case
fname = [dirname,'/','combined','.m'];
savecase(fname,mpcout);

%% dyr file for the integrated T-D case
%% Copy-over the tcase dyr file
combinedfilename_dyr = [dirname,'/','combined','.dyr'];
copyfile('datafiles/caseACTIVSg200t_genonly.dyr',combinedfilename_dyr);

fp = fopen(combinedfilename_dyr,'a');
[nbust,colt] = size(mpct.bus);
loadbusds = find(mpc.bus(:,1) > max(mpct.bus(:,1)) & (mpc.bus(:,PD)~=0 | mpc.bus(:,QD)~=0));

%% ZIP
%loaddyrdata = ['''ZIP'',1 ,0.0,0.0,0.0,0.0,100.0,100.0,0.85,/USRWHT'];
%% Composite (ZIP + IM) %% Note disabled motor relay triggering by using -0.7 as threshold.
loaddyrdata = '''COMPLOAD'',1 ,20.0,20.0,20.0,20.0,30.0,60.0,0,1.3000,14.0000,240.0000,0.9000,12.0000,999,999,0,0,0,0,100.0,0.0,1.5,-0.7,20,1.0,1.0,20,0.0,/USRWHT';

for i = 1:length(loadbusds)
    loaddata = [num2str(mpc.bus(loadbusds(i),1)),',',loaddyrdata];
    fprintf(fp,'%s',loaddata);
    fprintf(fp,'\n');
end
fclose(fp);

nbrancht = length(mpct.branch(:,1));
ctr = 0;
for i = 1:nbdry_buses
    td_branch = find(mpc.branch(:,F_BUS) == mpcdbs{i}.bdry_bus);
    td_branch = td_branch(td_branch > nbrancht);

    Pftd = sum(mpc.branch(td_branch,PF));
    Qftd = sum(mpc.branch(td_branch,QF));

    tidx = find(mpct.bus(:,1) == mpcdbs{i}.bdry_bus);
    mpct.bus(tidx,PD) = Pftd;
    mpct.bus(tidx,QD) = Qftd;
    
    for k = 1:mpcdbs{i}.nbdry
        % Network data
        mpcdbs{i}.dfilename{k} = [dirname,'/dcase_tbdry_',num2str(mpcdbs{i}.bdry_bus),'_feeder_',num2str(k),'.m'];
        % Topic for distribution federate (This is used by the C code for
        % assigning a name for the distribution subscription topic
        mpcdbs{i}.dtopic{k} = ['dcase_tbdry_',num2str(mpcdbs{i}.bdry_bus),'_feeder_',num2str(k)];
    end            
end


%% This is the data set for T only case where the load at the t_bdry_bus equals the load injected into the distribution system.
tfilename = [dirname,'/tcase.m'];
savecase(tfilename,mpct);
%% This is the dyr file for the T case
tfilename_dyr = [dirname,'/tcase.dyr'];
copyfile('datafiles/caseACTIVSg200t_genonly.dyr',tfilename_dyr);

tloadbdyrdata = ['''ZIP'',1 ,0.0,0.0,0.0,0.0,100.0,100.0,0.85,/USRWHT'];
fp = fopen(tfilename_dyr,'a');
for i = 1:nbdry_buses
    loaddata = [num2str(mpcdbs{i}.bdry_bus),',',tloadbdyrdata];
    fprintf(fp,'%s',loaddata);
    fprintf(fp,'\n');
end
fclose(fp);

%% Save all the distribution data files
for i = 1:nbdry_buses
    for k = 1:mpcdbs{i}.nbdry
        savecase(mpcdbs{i}.dfilename{k},mpcdbs{i}.datas{k});
        % dyr data
        dfilename_dyr = [dirname,'/dcase_tbdry_',num2str(mpcdbs{i}.bdry_bus),'_feeder_',num2str(k),'.dyr'];
        
        fp = fopen(dfilename_dyr,'w');
        mpcdb = mpcdbs{i}.datas{k};
        loadbusds = find(mpcdb.bus(:,PD)~=0 | mpcdb.bus(:,QD)~=0);

        for j = 1:length(loadbusds)
            loaddata = [num2str(mpcdb.bus(loadbusds(j),1)),',',loaddyrdata];
            fprintf(fp,'%s',loaddata);
            fprintf(fp,'\n');
        end
        fclose(fp);
    end
end
dfilename_event = [dirname,'/dcase_tbdry_',num2str(mpcdbs{1}.bdry_bus),'_feeder_',num2str(1),'.event'];                   
%% 

faultdata = '''FAULT'',1.0,1.2,0.0,99999';
%% Files for both transmission and distribution side faults are created
%% Transmission fault
%% File for both integrated simulation and co-simulation
fault_bus = 5;
combinedfilename_event = [dirname,'/','tfault.event'];
fp = fopen(combinedfilename_event,'w');
fault=[num2str(fault_bus),',',faultdata];
fprintf(fp,'%s',fault);
fclose(fp);

%% Distribution fault
%% Integrated simulation file
combinedfilename_event = [dirname,'/','combined_dfault.event'];
fp = fopen(combinedfilename_event,'w');
fault_bus = max(mpct.bus(:,1)) + 5; % max(mpct.bus(:,1)) gives the offset, 5 is the distribution bus number
fault=[num2str(fault_bus),',',faultdata];
fprintf(fp,'%s',fault);
fclose(fp);
%% Co-simulation file
cosimfilename_event = [dirname,'/','cosim_dfault.event'];
fp = fopen(cosimfilename_event,'w');
fault_bus = 5; % Fault at bus 5 on the distribution feeder
fault=[num2str(fault_bus),',',faultdata];
fprintf(fp,'%s',fault);
fclose(fp);

%copyfile('datafiles/caseACTIVSg200t_141d.event',combinedfilename_event);



%% Create metadata file for tranmsmission federate.
%% This file describes the number of distribution federates, their topics (same as dist file name)
fp = fopen([dirname,'/','metadata.trans'],'w');


fprintf(fp,'%d,%d\n',nbdry_buses,nfeeders);
for i = 1:nbdry_buses
    fprintf(fp,'%d\n',mpcdbs{i}.bdry_bus);
end
for i = 1:nbdry_buses
    for k = 1:mpcdbs{i}.nbdry
        fprintf(fp,'%d,%s\n',mpcdbs{i}.bdry_bus,mpcdbs{i}.dtopic{k});
    end
end
fclose(fp);

fprintf('dirname = %s\n',dirname);


%% TO DO LATER
%% Add solar
sper = 0.1; %% 10% per solar added
loadd = find(mpc.bus(nbust+1:end,PD));

%% Extend gen