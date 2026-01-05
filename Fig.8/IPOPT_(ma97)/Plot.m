load('All_sol_classes_10000initial.mat','sol_classes','class_ids');
load('All_x0_LOS_f_10000initial.mat','All_rand_x0','All_x_LOS','All_f_LOS');
load('x_LOS_10000initial.mat','x_LOSs');
%% -----------------------------
% 3D Plot: initial points colored by converged class
%% -----------------------------
figure;
num_classes = max(sol_classes);
colors = lines(num_classes);
hold on;
h_legend = gobjects(num_classes*2,1);
legend_entries = cell(num_classes*2,1);
for k = 1:num_classes
    % ----------------------
    % Plot initial points for class k
    % ----------------------
    idx_init = find(sol_classes == k);
    h_legend(2*k-1) = plot(All_rand_x0(266,idx_init), All_rand_x0(267,idx_init), '.', ... %All_rand_x0(268,idx_init)
        'Color', colors(k,:), 'MarkerSize', 8);
    legend_entries{2*k-1} = ' Initial point';
    % ----------------------
    % Plot converged points for class k
    % ----------------------
    idx_conv = find(class_ids == k);
    h_legend(2*k) = plot(x_LOSs(266,idx_conv), x_LOSs(267,idx_conv),'x', ...
        'MarkerEdgeColor', colors(k,:), 'MarkerSize', 10, 'LineWidth', 2);
    legend_entries{2*k} = 'Converged point';

end
legend(h_legend, legend_entries, 'Location', 'best');
xlabel('PG69'); ylabel('PG70'); zlabel('PG72');
title('Initial points and their converged solution of IPOPT-v3.14 "ma97"');
axis equal;
% grid on;


load('All_sol_classes_10000initial.mat', 'sol_classes')
%12963：
[~, index_of_LOS1] = find(sol_classes ==1);
LOS1_num = size(index_of_LOS1,2);
save('LOS1_idx.mat','index_of_LOS1','LOS1_num');

%17798：
[~, index_of_LOS2] = find(sol_classes ==3);
LOS2_num = size(index_of_LOS2,2);
save('LOS2_idx.mat','index_of_LOS2','LOS2_num');
%19570：
[~, index_of_LOS3] = find(sol_classes ==2);
LOS3_num = size(index_of_LOS3,2);
save('LOS3_idx.mat','index_of_LOS3','LOS3_num');