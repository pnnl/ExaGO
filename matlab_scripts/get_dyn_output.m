function [t,freq,Vm] = get_dyn_output(netfile,dyrfile,trajfolder)
%GET_DYN_OUTPUT Gets arrays for time, freq, and voltage from a dynamics
%simulation run using DYNXXX
%   This function takes in the network file and dyr file used with DYN and
%   returns the time vector, frequency vector, and bus voltage magnitude
%   vector
% INPUTS:
%   NETFILE - The network data file used in dynamics simulation
%   DYRFILE - The dynamics data file used in dynamics simulation
%  TRAJFOLDER - The folder in which trajectories are saved (typically
%  SA-data)
% NOTE - When running the dynamics simulation, one has to save the
% trajectories either by using option -dyn_ts_save_trajectory or using
% sensitivity analysis (that automatically saves the trajectories)
% NOTE - NEED TO RUN THIS FUNCTION FROM TOP LEVEL TSOPF DIRECTORY!!
% Example run:
%    [t,freq,Vm] = get_dyn_output('case9mod.m','case9mod.dyr','SA-data');

addpath('datafiles');
addpath('visualization/SampleData/DYNdata');

%%% Get the trajectory information in a struct
[results] = ReadTrajectoryBasic(trajfolder);

sz = results.size; %% number of time-points
t = results.t; %% time
x = results.x; %% solutions at each time-step



%% load network file
mpc = loadcase(netfile);

nbus = size(mpc.bus,1); %% Number of buses
bus_i = mpc.bus(:,1);
ngen = size(mpc.gen,1); %% Number of generators
ngenON = find(mpc.gen(:,8) ~= 0);

nvar = ones(nbus,1)*2; %% Initialize number of variables at each bus

dyn_tags = {'GENROU','IEEET1','EXST1','SEXS','TGOV1','STAB1','COMPLOAD','PVD1','CV','ZIP'}; % Tags for models in the dyr file
add_vars = [8,4,4,2,2,3,3,2,2,0]; %% Additional variables to be added at the node if the model is present
len_tags = length(dyn_tags);

%freq_idx = zeros(length(ngenON),1); %% indices for generator frequencies in the x vector
VD_idx = zeros(nbus,1);
VQ_idx = zeros(nbus,1);

mot_slip_idx = zeros(nbus,1);
genidx_arr = [];
%% load dyr file
fd = fopen(dyrfile);
tline = fgetl(fd);

while ischar(tline)
    while(strcmp(tline,'')) % ignore blank lines
       tline = fgetl(fd);
    end
    while(isempty(strfind(tline,'/')))
        tline = [tline,' ',fgetl(fd)];
    end
   % disp(tline);
    tag_found = 0;
    i = 0;
    while(~tag_found)
        i = i+1;
        if(i <= length(add_vars) && ~isempty(strfind(tline,dyn_tags{i})))
            tag_found = 1;
            if(~isempty(strfind(tline,'COMPLOAD')) || ~isempty(strfind(tline,'ZIP')))
                load_model=1;
            else
                load_model = 0;
            end
                
            if(~isempty(strfind(tline,',')))
                a = strsplit(tline,',');
            else
                a = strsplit(tline);
            end
            kk = 1;
            while(strcmp(a{kk},''))
                kk = kk + 1;
            end
            busnum = str2double(a{kk});
            genidx = find(busnum == mpc.gen(ngenON,1));
            if(~load_model && ~isempty(genidx))
                if(length(genidx) > 1)
                    error('Multiple generators at bus %d\n',busnum);
                end
           %     if(mpc.gen(genidx,8))
                    idx = find(busnum == bus_i);
                    nvar(idx) = nvar(idx) + add_vars(i);
                    rr = [];
                    if(~isempty(genidx_arr))
                        rr = find(genidx == genidx_arr);
                    end
                    if(isempty(rr))
                        genidx_arr = [genidx_arr,genidx];
                    
                        gens.idx(genidx) = genidx;
                        gens.busnum(genidx) = idx;
                        gens.model{genidx} = dyn_tags{i};
                    end
                    break;
        %        end
            else
                % Tag for load
                idx = find(busnum == bus_i);
                if(~isempty(strfind(a{2},'COMPLOAD')))
                    pl = str2double(a{4});
                    ip = str2double(a{6});
                    yp = str2double(a{8});
                    if(abs(pl+ip+yp - 100) > 1e-6)
                        mot_slip_idx(idx) = nvar(idx)+1;
                        nvar(idx) = nvar(idx) + add_vars(i);
                    end
                end
                break;
            end
        end
    end
    tline = fgetl(fd);
end

fclose(fd);

if(length(genidx_arr) ~= length(mpc.gen(ngenON,1)))
    %% Some of the generators are modeled as constant voltage sources
    [cvgen,idx1] = setdiff(mpc.gen(ngenON,1),mpc.gen(genidx_arr,1));
    for i=1:length(cvgen)
        idx = find(mpc.gen(cvgen(i),1) == bus_i);
        nvar(idx) = nvar(idx) + 2;
        gens.idx(idx1) = idx1;
        gens.busnum(idx1) = cvgen(i);
        gens.model{idx1} = 'CV';
    end
end

yy = cumsum(nvar);

idx = find(mot_slip_idx);
if(~isempty(idx))
    mot_slip_var_idx = yy(idx-1)+mot_slip_idx(idx);
end
VD_idx(1) = 1;
VD_idx(2:end) = yy(1:end-1)+1;
VQ_idx = VD_idx + 1;

Vm = sqrt(x(VD_idx,:).^2 + x(VQ_idx,:).^2);

ctr = 1;
freq_idx = [];
%% Assume all GENROUs....THIS IS HACK AND WILL NOT WORK IF OTHER GENERATOR MODELS ARE PRESENT
for i = 1:length(ngenON)
    if(strcmp(gens.model{i},'PVD1') | strcmp(gens.model{i},'CV'))
        continue;
    end
    idx_bus = find(mpc.gen(ngenON(i),1) == bus_i);
    
    if ~isempty(idx_bus)
        freq_idx(ctr) = VQ_idx(idx_bus) + 6;
        ctr = ctr + 1;
    end
end

freq = (1+x(freq_idx,:))*60;

rmpath('datafiles');

function [varargout] = ReadTrajectoryBasic(inarg)
%
%   [varargout] = PetscBinaryReadTrajectory(inarg)
%
%  Reads in the trajectory information saved in PETSc binary files and
%  emits as Matlab struct.
%
%  Examples: A = PetscBinaryReadTrajectory('myfile'); read from file 
%

indices = 'int32';
precision = 'float64';
maxsteps = 10000;

t = zeros(1,maxsteps);

for stepnum=1:maxsteps
  filename = sprintf('SA-%06d.bin',stepnum);
  fullname = fullfile(inarg,filename);
  if exist(fullname,'file') ~= 2
    size = stepnum-1;
    break;
  end
  fd = PetscOpenFile(fullname);
  header = double(read(fd,1,indices));  
  
  if isempty(header)
    size = stepnum-1;
    break;
  end

  if  header == 1211214 % Petsc Vec Object
    %% Read state vector   
    m = double(read(fd,1,indices));
    if (stepnum == 1)
      x = zeros(m,maxsteps);
    end
    v = read(fd,m,precision);
    x(:,stepnum) = v;

    %% Read time
    t(stepnum) = read(fd,1,precision);
  end
  % close the reader if we opened it
  close(fd);
end

result.size = size;
if size > 1
  result.t = t(1:size);
  result.x = x(:,1:size);
end
varargout{1} = result;


end

end

