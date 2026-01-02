clear; clc;


[mpc, mpopt] = opf_args('pglib_opf_case240_pserc');

% Set solver options
mpopt = mpoption; 

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% add zero columns to bus, gen, branch for multipliers, etc if needed
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
if size(mpc.bus,2) < MU_VMIN
  mpc.bus = [mpc.bus zeros(nb, MU_VMIN-size(mpc.bus,2)) ];
end
if size(mpc.gen,2) < MU_QMIN
  mpc.gen = [ mpc.gen zeros(ng, MU_QMIN-size(mpc.gen,2)) ];
end
if size(mpc.branch,2) < MU_ANGMAX
  mpc.branch = [ mpc.branch zeros(nl, MU_ANGMAX-size(mpc.branch,2)) ];
end
%%-----  convert to internal numbering, remove out-of-service stuff  -----
mpc = ext2int(mpc, mpopt);

load('case240_obj0Init_feasiblepoint(_fromPGSHEobj0FlatSatrt).mat', 'x0');

total_computation_time_mips = 0;
total_computation_time_IPOPT_mumps = 0;
total_computation_time_IPOPT_ma97 = 0;
total_computation_time_KNITRO_IP = 0;


%% mips:
fprintf('---------------- MIPS for 2383-bus OPF ----------------\n');
mpopt = mpoption('exp.use_legacy_core', 1, ...
    'opf.ac.solver', 'MIPS', ...
    'opf.start', 5, ...%input initial point
    'verbose', 3,...
    'mips.max_it', 3000); %'max_it', 150(default)
om = opf_setup(mpc, mpopt);
tic;
[results, success, raw] = nlpopf_solver_For_arbitrary_initpoint(x0, om, mpopt);
total_computation_time_mips = total_computation_time_mips + toc;
f_con_evaluation = raw.output.iterations + 1;
fprintf('objective and constraint evaluation: %d\n',f_con_evaluation);
fprintf('total computation time of mips: %d\n',total_computation_time_mips);


%% ipopt linear_solver: mumps
fprintf('---------------- IPOPT-mumps for 2383-bus OPF result ----------------\n');
mpopt = mpoption('exp.use_legacy_core', 1, ...
    'opf.ac.solver', 'IPOPT', ...
    'opf.start', 5, ... %input initial point
    'ipopt.opts.print_level', 0, ...
    'verbose', 2, ...
    'ipopt.opts.max_iter',3000, ...
    'ipopt.opts.mu_strategy','adaptive', ...
    'ipopt.opts.linear_solver', 'mumps');%mumps
om = opf_setup(mpc, mpopt);
tic;
[results, success, raw] = nlpopf_solver_For_arbitrary_initpoint(x0, om, mpopt);
total_computation_time_IPOPT_mumps = total_computation_time_IPOPT_mumps + toc;
fprintf('total computation time of IPOPT-mumps: %d\n',total_computation_time_mips);


%% ipopt linear_solver: ma97
fprintf('---------------- IPOPT-ma97 for 2383-bus OPF result ----------------\n');
mpopt = mpoption('exp.use_legacy_core', 1, ...
    'opf.ac.solver', 'IPOPT', ...
    'opf.start', 5, ... %input initial point
    'ipopt.opts.print_level', 0, ...
    'verbose', 2, ...
    'ipopt.opts.max_iter',3000, ...
    'ipopt.opts.mu_strategy','adaptive', ...
    'ipopt.opts.hsllib', 'libhsl.dll', ...
    'ipopt.opts.linear_solver', 'ma97'); %ma97
om = opf_setup(mpc, mpopt);
tic;
[results, success, raw] = nlpopf_solver_For_arbitrary_initpoint(x0, om, mpopt);
total_computation_time_IPOPT_ma97 = total_computation_time_IPOPT_ma97 + toc;
fprintf('total computation time of IPOPT-ma97: %d\n',total_computation_time_IPOPT_ma97);
x_ipopt_LOS = results.x;
f_ipopt_LOS = results.f;
% save('x_ipopt_LOS.mat','x_ipopt_LOS','f_ipopt_LOS');

%% knitro: Interior Point Method
fprintf('---------------- KNITRO- Interior Point Method for 2383-bus OPF result ----------------\n');
mpopt = mpoption('exp.use_legacy_core', 1,...
    'opf.start', 5,...
    'opf.ac.solver', 'KNITRO',...
    'verbose',3);
mpopt = mpoption(mpopt, ...
    'knitro.opts.algorithm',    1, ...  % 1/2 Barrier IPM；3 Active-set；4 SQP
    'knitro.maxit',             3000 ...
     );% ref https://www.artelys.com/app/docs/knitro/3_referenceManual/knitromatlabReference.html#setting-options
om = opf_setup(mpc, mpopt);
tic;
[results, success, raw] = nlpopf_solver_For_arbitrary_initpoint(x0, om, mpopt); 
total_computation_time_KNITRO_IP = total_computation_time_KNITRO_IP + toc;
fprintf('total computation time of KNITRO-IP: %d\n',total_computation_time_KNITRO_IP);

%% knitro: SQP Method
fprintf('---------------- KNITRO- SQP Method for 2383-bus OPF result ----------------\n');
mpopt = mpoption('exp.use_legacy_core', 1,...
    'opf.start', 5,...
    'opf.ac.solver', 'KNITRO',...
    'verbose',3);
mpopt = mpoption(mpopt, ...
    'knitro.opts.algorithm',    4, ...  % 1/2 Barrier IPM；3 Active-set；4 SQP
    'knitro.maxit',             3000 ...
     );% ref https://www.artelys.com/app/docs/knitro/3_referenceManual/knitromatlabReference.html#setting-options
om = opf_setup(mpc, mpopt);
tic;
[results, success, raw] = nlpopf_solver_For_arbitrary_initpoint(x0, om, mpopt); 
total_computation_time_KNITRO_IP = total_computation_time_KNITRO_IP + toc;
fprintf('total computation time of KNITRO-SQP: %d\n',total_computation_time_KNITRO_IP);