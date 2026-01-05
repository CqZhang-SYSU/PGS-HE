clear; clc;

%% -----------------------------
% Parameters
%% -----------------------------
mpc = loadcase('case118mod');
N = 10000; % Number of random initial points
abs_tol = 1e-2;
rel_tol = 1e-2;

% Initialize storage
x_LOSs = [];
All_x_LOS = [];
All_f_LOS = [];
All_f_success = [];
All_rand_x0 = [];
sol_classes = [];   % Class index for each initial point
class_ids = [];
next_class_id = 1;

% Statistics variables
index_success = [];
index_failure = [];
index_exceedMaxIter = [];
success_count = 0;
failure_count = 0;
iteration_limit_count = 0;
total_iterations = 0;
total_computation_time = 0;

%% -----------------------------
% Main loop
%% -----------------------------
for j = 1:N
    % Set IPOPT solver options
    mpopt = mpoption('exp.use_legacy_core', 1, ...
                     'opf.ac.solver', 'IPOPT', ...
                     'opf.start', 4, ...
                     'ipopt.opts.print_level', 0, ...
                     'verbose', 1, ...
                     'ipopt.opts.max_iter',3000, ...
                     'ipopt.opts.mu_strategy','adaptive', ...
                     'ipopt.opts.hsllib', 'libhsl.dll', ...
                     'ipopt.opts.linear_solver', 'ma97'); %ma97 or ma57

    % Run OPF
    results = runopf(mpc, mpopt);

    % Save initial point and optimal solution
    x0 = results.x0;
    x = results.x;
    f = results.f;

    All_rand_x0 = [All_rand_x0, x0]; 
   
    All_x_LOS = [All_x_LOS,x];
    %% -----------------------------
    % Classify converged solution
    %% -----------------------------
    if results.success
        found = false;

        if isempty(x_LOSs)
            % First converged solution
            x_LOSs = x;  
            All_f_LOS = f;
            class_ids = 1;             
            sol_classes = [sol_classes, 1];
            next_class_id = 2;         
        else
            % Compare with existing classes
            for k = 1:size(x_LOSs,2)
                dx = x(237:290,:) - x_LOSs(237:290,k);
                % Absolute or relative tolerance check
                if norm(dx,2) <= abs_tol || norm(dx,2) <= rel_tol*max(1,norm(All_x_LOS(237:290,k),2))
                    % Belongs to an existing class
                    sol_classes = [sol_classes, class_ids(k)];
                    found = true;
                    break;
                end
            end

            if ~found
                % New class
                x_LOSs = [x_LOSs, x];     
                All_f_LOS = [All_f_LOS;f];
                class_ids = [class_ids, next_class_id];  
                sol_classes = [sol_classes, next_class_id];
                next_class_id = next_class_id + 1;       
            end
        end
    else
        % Not converged, mark as -1
        sol_classes = [sol_classes, -1];
    end

    %% -----------------------------
    % Update statistics
    %% -----------------------------
    status = results.raw.output.status;
    if status == 0 || status == 1
        success_count = success_count + 1;
        index_success = [index_success; j];
         All_f_success = [All_f_success;f];
    elseif status == -1
        iteration_limit_count = iteration_limit_count + 1;
        index_exceedMaxIter = [index_exceedMaxIter; j];
    else
        failure_count = failure_count + 1;
        index_failure = [index_failure; j];
    end

    total_iterations = total_iterations + 1;

end

%% -----------------------------
% Save results
%% -----------------------------
% save('All_sol_classes_10000initial.mat','sol_classes','class_ids');
% save('All_x0_LOS_f_10000initial.mat','All_rand_x0','All_x_LOS','All_f_LOS');
% save('x_LOS_10000initial.mat','x_LOSs');

