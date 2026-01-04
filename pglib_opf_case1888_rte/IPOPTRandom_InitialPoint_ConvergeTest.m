clear; clc;

[mpc, mpopt] = opf_args('pglib_opf_case1888_rte');

% Set solver options
mpopt = mpoption; 

% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

% add zero columns to bus, gen, branch for multipliers, etc if needed
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



%% -----------------------------
% Parameters
%% -----------------------------

N = 100; % Number of random initial points
% Initialize storage
All_x_LOS = [];
All_f_LOS = [];

% Statistics variables
index_success = [];
index_failure = [];
failure_status = [];
success_count = 0;
failure_count = 0;
total_iterations = 0;
total_computation_time = 0;


%% -----------------------------
% Main loop
%% -----------------------------
% dim = 4356;
% All_x0  = rand(dim,N); % Completely random initial points: N = 100, dim = number of variables in "case"
% save("All_x0_CompletelyRandom.mat",'All_x0');
load("All_x0_CompletelyRandom.mat",'All_x0');

for j = 1:size(All_x0,2)
    % Set IPOPT solver options
    mpopt = mpoption('exp.use_legacy_core', 1, ...
        'opf.ac.solver', 'IPOPT', ...
        'opf.start', 5, ... % Input the specified initial value
        'ipopt.opts.print_level', 0, ...
        'verbose', 2, ...
        'ipopt.opts.linear_solver', 'ma97',... % or mumps
        'ipopt.opts.hsllib', 'libhsl.dll',...
        'ipopt.opts.mu_strategy','adaptive',...
        'ipopt.opts.max_iter', 3000); 

    % Run OPF
    tic;
    om = opf_setup(mpc, mpopt);
    x0 = All_x0(:,j);
    [results, success, raw] = nlpopf_solver_For_arbitrary_initpoint(x0, om, mpopt);
    total_computation_time = total_computation_time + toc;

    % Save initial point and optimal solution
    x = results.x;
    f = results.f;

    All_f_LOS = [All_f_LOS;f];
    All_x_LOS = [All_x_LOS,x];

    %% -----------------------------
    % Update statistics
    %% -----------------------------
    % nlps_ipopt.m EXITFLAG : exit flag
    %           1 = converged
    %           0 = failed to converge
    %           OUTPUT : output struct with the following fields:
    %           status - see IPOPT documentation for INFO.status
    % https://coin-or.github.io/Ipopt/IpReturnCodes__inc_8h_source.html
    status = raw.output.status;
    if status == 0
        success_count = success_count + 1;
        index_success = [index_success; j];
    else
        failure_count = failure_count + 1;
        index_failure = [index_failure; j];
        failure_status = [failure_status;status];
    end

    total_iterations = total_iterations + 1;

end
% Calculate average computation time for the current case
average_computation_time = total_computation_time / N;
% Store the results for the current case
convergence_results.success_count = success_count;
convergence_results.failure_count = failure_count;
convergence_results.total_iterations = total_iterations;
convergence_results.average_computation_time = average_computation_time;
%% -----------------------------
% Save results
%% -----------------------------
% Display summary for each case
fprintf('Results for case1888rte \n');
fprintf('Success Count: %d\n', success_count);
fprintf('Failure Count: %d\n', failure_count);
fprintf('Total Iterations: %d\n', total_iterations);
fprintf('Average Computation Time: %d\n\n', average_computation_time);
% save('All_x_LOS_f.mat','All_x_LOS','All_f_LOS');
% save('All_success_failure_index.mat','index_success','index_failure','failure_status');



%% Note: status
% {
%    Solve_Succeeded                    = 0,
%    Solved_To_Acceptable_Level         = 1,
%    Infeasible_Problem_Detected        = 2,
%    Search_Direction_Becomes_Too_Small = 3,
%    Diverging_Iterates                 = 4,
%    User_Requested_Stop                = 5,
%    Feasible_Point_Found               = 6,
% 
%    Maximum_Iterations_Exceeded        = -1,
%    Restoration_Failed                 = -2,
%    Error_In_Step_Computation          = -3,
%    Maximum_CpuTime_Exceeded           = -4,
%    Maximum_WallTime_Exceeded          = -5,   
% 
%    Not_Enough_Degrees_Of_Freedom      = -10,
%    Invalid_Problem_Definition         = -11,
%    Invalid_Option                     = -12,
%    Invalid_Number_Detected            = -13,
% 
%    Unrecoverable_Exception            = -100,
%    NonIpopt_Exception_Thrown          = -101,
%    Insufficient_Memory                = -102,
%    Internal_Error                     = -199
%    }



