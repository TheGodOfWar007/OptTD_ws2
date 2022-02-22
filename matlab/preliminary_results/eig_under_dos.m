% This file was created to test the hypothesis that the eigenvalues will
% decrease for any general digraph with positive weights under the
% unfluence of a DoS attack. For now the link faliure function is assumed
% to be direct delta => instant link failure thus neglecting transients. A
% transient analysis will be carried out in the future.
close all;
%% IMPORTANT PARAMETERS
NUM_NODES = 6;
W_MIN = 0; %-ve weights will be zeroed out.
W_MAX = 1;
ENABLE_SYMMETRY = false;
n = NUM_NODES; % Alias
cd = color_dict();
%% Initializing a random graph.
A = W_MIN + (W_MAX - W_MIN).*rand(n);
A(A < 0.5) = 0; 
A = A - diag(diag(A));
if ENABLE_SYMMETRY
    A = (A + A')/2;
end
L = diag(sum(A,2)) - A;
G = digraph(A);
% figure(1);
% plot(G);
%% Epsilon bounded perturbation due to DoS attack. Analysis Type 1 - Discrete Variation.
% TYPE 1: Randomly removing links one-by-one and checking eigenvalue
% variation. Also simulating gersgorin disk radius variation.
close all;
A1 = A;
G1 = G;
nz = find((A > 0));
iterations = numel(nz);
L1_seq = zeros(n,n,iterations+1);
A1_seq = zeros(n,n,iterations+1);
L1_seq(:,:,1) = L;
A1_seq(:,:,1) = A;
L1_norm_seq = zeros(1,iterations+1);
L1_norm_seq(1) = norm(L, 'fro');
L1_rank_seq = zeros(1, iterations);
g_disk_r = zeros(1, iterations);
for i = 1:numel(nz)
    nz = find((A1 > 0));
    nz = nz(randperm(numel(nz)));
    rm_idx = nz(1);
    A1(rm_idx) = 0;
    L1 = graph_laplacian(A1);
    L1_seq(:,:,i) = L1;
    A1_seq(:,:,i) = A1;
    L1_norm_seq(i) = norm(L1, 'fro');
    L1_rank_seq(i) = rank(L1);
    g_disk_r(i) = max(diag(L1));
end

%%% Obtaining the eigenvalue sequences using the eigenshuffle function from
%%% MATLAB FEX.
[V1_seq, D1_seq] = eigenshuffle(L1_seq);
D1_real = real(D1_seq);
D1_imag = imag(D1_seq);
D1_abs = abs(D1_seq);

%%% Plotting the variation sequence.
plot_colors = [cd.simulink_blue; cd.simulink_green; cd.simulink_red; cd.simulink_violet; cd.simulink_brown; cd.simulink_cyan];
for i=1:n
    eig_real = D1_real(i,:);
    eig_imag = D1_imag(i,:);
    eig_abs = D1_abs(i,:);
    d_eig_real = horzcat(diff(eig_real),0);
    d_eig_imag = horzcat(diff(eig_imag),0);
    figure(2);
    subplot(2,2,1)
    hold on
    plot(eig_real, eig_imag, 'LineWidth', 1.5, 'Color', plot_colors(i,:),'Marker','.');
    title('2D Plot of Eigenvalue Variation', 'Interpreter','latex');
    grid on;
    subplot(2,2,2)
    hold on
    plot(eig_abs, 'LineWidth',1.5,'Color',plot_colors(i,:),'Marker','.');
    title('$|\lambda|$ v/s Iteration', 'Interpreter','latex');
    grid on;
    subplot(2,2,3);
    hold on;
    plot(L1_norm_seq, eig_abs, 'LineWidth',1.5, 'Color', plot_colors(i,:),'Marker','.');
    title('$|\lambda|$ v/s $\|L\|_{fr}$','Interpreter','latex');
    grid on;
    figure(4);
    subplot(2,2,1)
    hold on
    plot(eig_real, 'LineWidth',1.5,'Color',plot_colors(i,:),'Marker','.');
    title('$Re(\lambda)$ v/s Iteration', 'Interpreter','latex');
    grid on;
    subplot(2,2,2);
    hold on;
    plot(L1_norm_seq, eig_real, 'LineWidth',1.5, 'Color', plot_colors(i,:),'Marker','.');
    title('$Re(\lambda)$ v/s $\|L\|_{fr}$','Interpreter','latex');
    grid on;
end

figure(2);
subplot(2,2,4)
plot(L1_norm_seq, 'LineWidth',1,'Color','k','Marker','.');
title('$\|L\|_{fr}$ v/s Iteration','Interpreter','latex');
grid on;
hold off;
figure(4)
subplot(2,2,3)
hold on;
it_vec = 1:iterations;
plot(it_vec, L1_rank_seq, 'LineWidth',1.5, 'Color', 'b','Marker','.');
title('$rank(L)$ v/s Iterations','Interpreter','latex');
grid on;
g_disk_r_max = max(g_disk_r);
for i=1:iterations
    eig_real = D1_real(:,i);
    eig_imag = D1_imag(:,i);
    figure(3)
    hold all;
    plt = plot(eig_real, eig_imag, '.', 'MarkerSize',10);
    title('Variation in Largest Gersgorin Disk','Interpreter','latex');
    viscircles([g_disk_r(i),0],g_disk_r(i), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', plt.Color);
    xlim([0 Inf]);
    axis equal;
end
%% Epsilon bounded perturbation due to DoS attack. Analysis Type 1 - Continuous Variation.
% TYPE 1: Randomly removing links one-by-one and checking eigenvalue
% variation continuous by linear link removal with time at a sampling
% frequency fs. Set fs according to the speed of the algorithm.
close all;
fs = 100; % The step size is 1/fs.
ts = 1/fs;
A2 = A;
G2 = G;
nz = find((A > 0));
num_nz = numel(nz);
iterations = floor(sum(A2,"all")*fs);
L2_seq = zeros(n,n,iterations);
A2_seq = zeros(n,n,iterations);
L2_seq(:,:,1) = L;
A2_seq(:,:,1) = A;
L2_norm_seq = zeros(1,iterations);
L2_norm_seq(1) = norm(L, 'fro');
g_disk_r = zeros(1, num_nz);
j = 1;
it_per_link = [1];
for i = 1:num_nz
    nz = find((A2 > 0));
    nz = nz(randperm(numel(nz)));
    rm_idx = nz(1);
    % Continuous removal of link at fs.
    while A2(rm_idx) > 0
        A2(rm_idx) = A2(rm_idx) - ts;
        A2(A2 < 0) = 0;
        j = j + 1;
        if j > iterations
            break;
        end
        L2 = graph_laplacian(A2);
        L2_seq(:,:,j) = L2;
        A2_seq(:,:,j) = A2;
        L2_norm_seq(j) = norm(L2, 'fro');
    end
    it_per_link = [it_per_link, j];
    g_disk_r(i) = max(diag(L2));
end

%%% Obtaining the eigenvalue sequences using the eigenshuffle function from
%%% MATLAB FEX.
[V2_seq, D2_seq] = eigenshuffle(L2_seq);
D2_real = real(D2_seq);
D2_imag = imag(D2_seq);
D2_abs = abs(D2_seq);
feidler_eigenvalue = zeros(1,iterations);
%%% Extracting the feidler eigenvalue sequence at every timestep.
for k=1:length(D2_abs)
    vec = D2_abs(:,k);
    feidler_eigenvalue(k) = min(setdiff(vec, min(vec)));
end
feidler_eigenvalue_interx = [];

%%% Plotting the variation sequence.
plot_colors = [cd.simulink_blue; cd.simulink_green; cd.simulink_red; cd.simulink_violet; cd.simulink_brown; cd.simulink_cyan];
for i=1:n
    eig_real = D2_real(i,:);
    eig_imag = D2_imag(i,:);
    eig_abs = D2_abs(i,:);
    d_eig_real = horzcat(diff(eig_real),0);
    d_eig_imag = horzcat(diff(eig_imag),0);
    % Calculating the change points in the feidler eigenvalue.
    feidler_eigenvalue_interx = horzcat(feidler_eigenvalue_interx, InterX([feidler_eigenvalue;1:iterations], [eig_abs; 1:iterations]));
    figure(4);
    subplot(3,2,1)
    hold on
    plot(eig_real, eig_imag, 'LineWidth', 1, 'Color', plot_colors(i,:));
    title('2D Plot of Eigenvalue Variation', 'Interpreter','latex');
    grid on;
    subplot(3,2,2)
    hold on
    plot(eig_abs, 'LineWidth',1,'Color',plot_colors(i,:));
    title('$|\lambda|$ v/s Iteration', 'Interpreter','latex');
    grid on;
    subplot(3,2,3);
    hold on;
    plot(L2_norm_seq, eig_abs, 'LineWidth',1, 'Color', plot_colors(i,:));
    title('$|\lambda|$ v/s $\|L\|_{fr}$','Interpreter','latex');
    grid on;
    figure(6)
    subplot(1,2,1)
    hold all;
    plot(abs(eig_real), 'LineWidth', 1.5, 'Color', plot_colors(i,:));
    title('$|Re(\lambda)|$ v/s Iteration', 'Interpreter','latex');
    grid on;
    subplot(1,2,2)
    hold all;
    plot(abs(eig_imag), 'LineWidth', 1.5, 'Color', plot_colors(i,:));
    title('$|Im(\lambda)|$ v/s Iteration', 'Interpreter','latex');
    grid on;
end

figure(4);
subplot(3,2,4)
plot(L2_norm_seq, 'LineWidth',1,'Color','k');
title('$\|L\|_{fr}$ v/s Iteration','Interpreter','latex');
grid on;
subplot(3,2,5)
plot(horzcat(diff(L2_norm_seq),0), 'LineWidth',1,'Color','k');
title('$\Delta\|L\|_{fr}$ v/s Iteration','Interpreter','latex');
grid on;
subplot(3,2,6)
% plot(1:1:iterations, feidler_eigenvalue, feidler_eigenvalue_interx(1,:), feidler_eigenvalue_interx(2,:), 'bo', 'LineWidth',1);
plot(feidler_eigenvalue,'b', 'LineWidth',1);
title('$|\lambda_2|$ (Feidler Eigenvalue) v/s Iteration','Interpreter','latex');
grid on;
hold off;
g_disk_r_max = max(g_disk_r);
% Changing the default map for better visibility.
choice = 1; % Using Jet. Uncomment the line below to toggle the menu.
% choice = menu('Which ColorOrder do you want?', 'jet', 'random', 'hsv', 'hot', 'cool', 'spring', 'summer',...
% 	'autumn', 'winter', 'lines', 'gray', 'bone', 'copper', 'pink');
newColorMap = color_order(choice, num_nz);
set(gca, 'ColorOrder', newColorMap, 'NextPlot', 'replacechildren','Colormap',newColorMap);

% Now get the new set of default plot colors.
% Verify it changed by printing out the new default color set to the command window.
% Uncomment the line below to check the map.
% newColorOrder = get(gca,'ColorOrder');
for i=1:num_nz
    eig_real = D2_real(:,it_per_link(i):it_per_link(i+1)-1);
    eig_imag = D2_imag(:,it_per_link(i):it_per_link(i+1)-1);
    figure(5)
    hold all;
    plt = plot(eig_real(:), eig_imag(:), '.', 'MarkerSize',6);
    colororder(newColorMap)
    title('Variation in Largest Gersgorin Disk','Interpreter','latex');
    if i > 1
        if abs(g_disk_r(i) - g_disk_r(i-1)) < 1e-6
            continue
        end
    end
    viscircles([g_disk_r(i),0], g_disk_r(i), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', plt.Color);
    xlim([0 Inf]);
    axis equal;
end
% colororder(newColorMap)