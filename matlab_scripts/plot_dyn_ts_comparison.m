function [output] = plot_dyn_ts_comparison(netfile,dyrfile,eventfile)

addpath('datafiles');
addpath('datafiles/scalability');

filename = ['./DYN2 -netfile ' netfile ' -dyrfile ' dyrfile ' -eventfile ' eventfile ' -dyn_ts_save_trajectory -ts_trajectory_keep_files -ts_event_tol 1e-5'];
pat = {'\s+'};
ts_opt = ' -dyn_ts_final_time 5.0 -dyn_ts_dt 0.00833333 ';
ts_tol = ' -dyn_ts_atol 1e-2 -dyn_ts_rtol 1e-2 -dyn_ts_max_snes_failures 10 -dyn_ts_equation_type 1000';
dt_max = [' -dyn_ts_adapt_dt_max 0.04 -dyn_ts_adapt_dt_min 0.0083333'];
grep_what = [];%' |grep -i "TSStep"';
clear_SA_data = ' rm -rf SA-data';

arkimex{1} = ' -dyn_ts_type arkimex -dyn_ts_arkimex_type 2e -dyn_ts_adapt_type basic ';
arkimex{2} = ' -dyn_ts_type arkimex -dyn_ts_adapt_type basic ';
arkimex{3} = ' -dyn_ts_type arkimex -dyn_ts_arkimex_type 4 -dyn_ts_adapt_type basic ';

rosw{1} = ' -dyn_ts_type rosw -dyn_ts_rosw_type 2m -dyn_ts_adapt_type basic ';
rosw{2} = ' -dyn_ts_type rosw -dyn_ts_adapt_type basic ';
rosw{3} = ' -dyn_ts_type rosw -dyn_ts_rosw_type 4l -dyn_ts_adapt_type basic ';

theta{1} = ' -dyn_ts_type cn -dyn_ts_adapt_type none';
theta{2} = ' -dyn_ts_type cn -dyn_ts_adapt_type basic -dyn_ts_theta_adapt';

semiexplicit{1} = ' -dyn_ts_type rk -dyn_use_semiexplicit -dyn_ts_rk_type 2a -dyn_ts_adapt_type none ';
semiexplicit{2} = ' -dyn_ts_type rk -dyn_use_semiexplicit -dyn_ts_rk_type 3 -dyn_ts_adapt_type none ';

normtype{1} = ' -dyn_ts_adapt_wnormtype 2';
normtype{2} = ' -dyn_ts_adapt_wnormtype INFINITY';
log_summary = ' -log_summary ';

marker = {'c-','k-','kp-','o-','x-','+-','r*-','rs-','rd-.'};

solveropt = {semiexplicit{1} theta{1} theta{2} rosw{1} rosw{2} rosw{3} arkimex{1} arkimex{2} arkimex{3}};

nsolvers = size(solveropt,2);


legend_tag = {{'Alternating Explicit'},{'Fixed-step Trapezoidal','Variable-step Trapezoidal'},{'2^{nd} order ROSW','3^{rd} order ROSW','4^{th} order ROSW'},{'2^{nd} order ARKIMEX','3^{rd} order ARKIMEX','4^{th} order ARKIMEX'}};

title_tag = {'Alternating explicit' 'Implicit trapezoidal'  'Variable-step Rosenbrock' 'Variable-step Implicit RK'};

[VD_idx,VQ_idx,freq_idx,VR_idx] = get_indices(netfile,dyrfile);

t_2 = cell(nsolvers,1);
nsteps_2 = cell(nsolvers,1);
ctr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TWO-NORM %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Two-norm results\n');
for(j=1:nsolvers)
  [status,results] = system(clear_SA_data);  
  options = [ts_opt ts_tol dt_max normtype{1} solveropt{j} grep_what];
  [status,results] = system([filename log_summary options]);
  
  %% Grep total run time
  rline1 = regexp(results,'Time');
  if ~isempty(rline1)
    rline2 = regexp(results(rline1:end),'\n');
    totline = results(rline1:rline1+rline2(1));
    rline3 = regexp(results,'TSStep');
    rline4 = regexp(results(rline3:end),'\n');
    tsstepline = results(rline3:rline3+rline4(1));
   
    totstat = regexp(totline,pat,'split');
    tsstepstat = regexp(tsstepline,pat,'split');
%  if(length(stat{:}) == 16 || length(stat{:}) == 18 || length(stat{:}) == 20 || length(stat{:}) == 22)    
    ctr = ctr+1;
    % Get the solver time
    t_2{j} = str2num(totstat{1}{3});
    nsteps_2{j} = str2num(tsstepstat{1}{2});
    
    [t_vec,Vm_vec,freq_vec,VR_vec] = get_solution_vec(VD_idx,VQ_idx,freq_idx,VR_idx);
    
    
    if(j == 1)
        %% semiexplict
        p =1;
        
        [fmax,imax] = max(max(freq_vec,[],2));
        [vmin,imin] = min(min(Vm_vec,[],2));
        
        figure(4),plot(t_vec,VR_vec(:,:),'linewidth',4);
        h = gca;
        figure(4),set(h,'FontWeight','Bold','FontSize',16);
        figure(4),xlabel('Time (sec)');
        figure(4),ylabel('V_R');
      %  legend('Exciter 1','Exciter 2','Exciter 3');
        
    elseif(j >= 2 && j < 4)
        p = 2;
    elseif(j >= 4 && j < 7)
        p = 3;
    else
        p = 4;
    end
    
%     figure(4),subplot(2,2,p),plot(t_vec,VR_vec(:,:),'linewidth',4);
%     h = gca;
%     figure(4),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
%     figure(4),subplot(2,2,p),xlabel('Time (sec)');
%     figure(4),subplot(2,2,p),ylabel(['Gen. ' num2str(imax) ' V_R']);
%     figure(4),subplot(2,2,p),hold on;
%     figure(4),subplot(2,2,p),title(title_tag{p});
 %legend('Exciter 1','Exciter 2','Exciter 3');
    
    figure(1),subplot(2,2,p),plot(t_vec,(freq_vec(imax,:)+1)*60,marker{j});
    %figure(1),subplot(2,2,p),plot(t_vec,VR_vec(:,:),marker{j}); 
    h = gca;
    figure(1),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
    figure(1),subplot(2,2,p),xlabel('Time (sec)');
    figure(1),subplot(2,2,p),ylabel(['Gen. ' num2str(imax) ' freq (Hz.)']);
    figure(1),subplot(2,2,p),hold on;
    figure(1),subplot(2,2,p),title(title_tag{p});
    
    figure(2),subplot(2,2,p),plot(t_vec,Vm_vec(imin,:),marker{j});
    h = gca;
    figure(2),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
    figure(2),subplot(2,2,p),xlabel('Time (sec)');
    figure(2),subplot(2,2,p),ylabel(['Bus ' num2str(imin) ' Vm (p.u.)']);
    figure(2),subplot(2,2,p),hold on;
    figure(2),subplot(2,2,p),title(title_tag{p});
    
    %%% Time-steps taken by integration schemes
    dt_diff = diff(t_vec);
    diff_idx = find(dt_diff ~= 0);
    figure(3),subplot(2,2,p),plot(t_vec(diff_idx),dt_diff(diff_idx),marker{j});
    h = gca;
    figure(3),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
    figure(3),subplot(2,2,p),xlabel('Time (sec)');
    figure(3),subplot(2,2,p),ylabel('\Delta{t}');
    figure(3),subplot(2,2,p),hold on;
    figure(3),subplot(2,2,p),title(title_tag{p});
    figure(3),subplot(2,2,p),axis([0 5 0 0.05]);
    
  else
    t_2{j} = 0;
    nsteps_2{j} = 'Not converged';
  end
end
figure(1),subplot(2,2,1),legend(legend_tag{1});
figure(1),subplot(2,2,2),legend(legend_tag{2});
figure(1),subplot(2,2,3),legend(legend_tag{3});
figure(1),subplot(2,2,4),legend(legend_tag{4});
figure(2),subplot(2,2,1),legend(legend_tag{1});
figure(2),subplot(2,2,2),legend(legend_tag{2});
figure(2),subplot(2,2,3),legend(legend_tag{3});
figure(2),subplot(2,2,4),legend(legend_tag{4});
figure(3),subplot(2,2,1),legend(legend_tag{1});
figure(3),subplot(2,2,2),legend(legend_tag{2});
figure(3),subplot(2,2,3),legend(legend_tag{3});
figure(3),subplot(2,2,4),legend(legend_tag{4});


t_inf = cell(nsolvers,1);
nsteps_inf = cell(nsolvers,1);
ctr = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INFINITY-NORM %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Infinity-norm results\n');
for(j=1:nsolvers)
   options = [ts_opt ts_tol dt_max normtype{2} solveropt{j} grep_what];
  [status,results] = system([filename log_summary options]);
  
  stat = regexp(results,pat,'split');
  if(length(stat{:}) == 16 || length(stat{:}) == 18 || length(stat{:}) == 20 || length(stat{:}) == 22)
    ctr = ctr + 1;
    t_inf{j} = str2num(stat{1}{4});
    nsteps_inf{j} = str2num(stat{1}{2});
  
    % Get the output file
    [t_vec,Vm_vec,freq_vec] = get_solution_vec(VD_idx,VQ_idx,freq_idx);

    if(j <= 2)
        %% semiexplict
        p =1;
    elseif(j > 2 && j <= 4)
        p = 2;
    elseif(j > 4 && j <= 7)
        p = 3;
    else
        p = 4;
    end

    [fmax,imax] = min(min(freq_vec,[],2));
    figure(4),subplot(2,2,p),plot(t_vec,(freq_vec(imax,:)+1)*60,marker{j});
    h = gca;
    figure(4),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
    figure(4),subplot(2,2,p),xlabel('Time (sec)');
    figure(4),subplot(2,2,p),ylabel(['Gen. ' num2str(imax) ' freq. (Hz)']);
    figure(4),subplot(2,2,p),hold on;
    figure(4),subplot(2,2,p),title(title_tag{p});
    
    [vmin,imin] = min(min(Vm_vec,[],2));
    figure(5),subplot(2,2,p),plot(t_vec,Vm_vec(imin,:)+1,marker{j});
    h = gca;
    figure(5),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
    figure(5),subplot(2,2,p),xlabel('Time (sec)');
    figure(5),subplot(2,2,p),ylabel(['Bus ' num2str(imin) ' Vm (p.u.)']);
    figure(5),subplot(2,2,p),hold on;
    figure(5),subplot(2,2,p),title(title_tag{p});
    
    %%% Time-steps taken by integration schemes
    dt_diff = diff(t_vec);
    diff_idx = find(dt_diff ~= 0);
    figure(6),subplot(2,2,p),plot(t_vec(diff_idx),dt_diff(diff_idx),marker{j});
    h = gca;
    figure(6),subplot(2,2,p),set(h,'FontWeight','Bold','FontSize',16);
    figure(6),subplot(2,2,p),xlabel('Time (sec)');
    figure(6),subplot(2,2,p),ylabel('\Delta{t}');
    figure(6),subplot(2,2,p),hold on;
    figure(6),subplot(2,2,p),title(title_tag{p});
  else
    t_inf{j} = 0;
    nsteps_inf{j} = 'Not converged';
  end
end
figure(4),subplot(2,2,1),legend(legend_tag{1});
figure(4),subplot(2,2,2),legend(legend_tag{2});
figure(4),subplot(2,2,3),legend(legend_tag{3});
figure(4),subplot(2,2,4),legend(legend_tag{4});
figure(5),subplot(2,2,1),legend(legend_tag{1});
figure(5),subplot(2,2,2),legend(legend_tag{2});
figure(5),subplot(2,2,3),legend(legend_tag{3});
figure(5),subplot(2,2,4),legend(legend_tag{4});
figure(6),subplot(2,2,1),legend(legend_tag{1});
figure(6),subplot(2,2,2),legend(legend_tag{2});
figure(6),subplot(2,2,3),legend(legend_tag{3});
figure(6),subplot(2,2,4),legend(legend_tag{4});

output.solver = '';

output.solver = {'2^nd order Explicit RK','3^rd order Explicit RK','Fixed-step Trapezoidal','Variable-step Trapezoidal','2^{nd} order ROSW','3^{rd} order ROSW','4^{th} order ROSW','2^{nd} order ARKIMEX','3^{rd} order ARKIMEX','4^{th} order ARKIMEX'};
output.t_2 = t_2;
output.nsteps_2 = nsteps_2;
output.t_inf = t_inf;
output.nsteps_inf = nsteps_inf;

end

function [t_vec,Vm_vec,freq_vec,VR_vec] = get_solution_vec(VD_idx,VQ_idx,freq_idx,VR_idx)

%%% Get the trajectory information in a struct
[results] = ReadTrajectoryBasic('SA-data');

sz = results.size; %% number of time-points
t_vec = results.t; %% time
x = results.x; %% solutions at each time-step

Vm_vec = sqrt(x(VD_idx,:).^2 + x(VQ_idx,:).^2);
freq_vec = x(freq_idx,:);
VR_vec = x(VR_idx,:);

%% Plot generator frequency deviation
%figure(1),plot(t,x(freq_idx,:)*60+60);
%axis([0 1 -0.01 0.015])

%figure(2),plot(t,Vm);

end

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

function [VD_idx,VQ_idx,freq_idx,VR_idx] = get_indices(netfile,dyrfile)

    %% load network file
    mpc = loadcase(netfile);

    nbus = size(mpc.bus,1); %% Number of buses
    bus_i = mpc.bus(:,1);
    ngen = size(mpc.gen,1); %% Number of generators
    ngenON = find(mpc.gen(:,8) ~= 0);

    nvar = ones(nbus,1)*2; %% Initialize number of variables at each bus

    dyn_tags = {'GENROU','IEEET1','EXST1','TGOV1','STAB1'}; % Tags for models in the dyr file
    add_vars = [8,4,4,2,3]; %% Additiona variables to be added at the node if the model is present
    len_tags = length(dyn_tags);

    freq_idx = zeros(length(ngenON),1); %% indices for generator frequencies in the x vector
    VR_idx = zeros(length(ngenON),1);
    VD_idx = zeros(nbus,1);
    VQ_idx = zeros(nbus,1);

    %% load dyr file
    fd = fopen(dyrfile);
    tline = fgetl(fd);
    while ischar(tline)
        disp(tline);
        for i = 1:len_tags
            if(~isempty(strfind(tline,dyn_tags{i})))
                a = strsplit(tline,',');
                busnum = str2double(a(1));
                genidx = busnum == mpc.gen(:,1);
                if(mpc.gen(genidx,8))
                    idx = find(busnum == bus_i);
                    nvar(idx) = nvar(idx) + add_vars(i);
                    break;
                end
            end
        end
        tline = fgetl(fd);
    end

    fclose(fd);

    yy = cumsum(nvar);
    VD_idx(1) = 1;
    VD_idx(2:end) = yy(1:end-1)+1;
    VQ_idx = VD_idx + 1;

    %% Assume all GENROUs....THIS IS HACK AND WILL NOT WORK IF OTHER GENERATOR MODELS ARE PRESENT
    for i = 1:length(ngenON)
        idx_bus = find(mpc.gen(ngenON(i),1) == bus_i);
        if ~isempty(idx_bus)
            freq_idx(i) = VQ_idx(idx_bus) + 6;
            VR_idx(i) = freq_idx(i) + 4;
        end
    end

end


