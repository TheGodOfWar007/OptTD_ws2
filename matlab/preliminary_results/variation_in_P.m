close all;
%% IMPORTANT PARAMETERS
NUM_NODES = 7;
TOL = 1e-6;
COND_TOL = 1;
W_MIN = 0; %-ve weights will be zeroed out.
W_MAX = 1;
ENABLE_SYMMETRY = false;
kQ = 2;
ALPHA = 0.2;

% Aliases
n = NUM_NODES;

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
Adim = size(A);
Q = diag(repelem(kQ, Adim(1)));
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
g_disk_r = zeros(1, iterations);
i = 1;
while i <= iterations
    nz = find((A1 > 0));
    nz = nz(randperm(numel(nz)));
    rm_idx = nz(1);
    A1(rm_idx) = 0;
    L1 = graph_laplacian(A1);
    L1_seq(:,:,i) = L1;
    A1_seq(:,:,i) = A1;
    L1_norm_seq(i) = norm(L1, 'fro');
    g_disk_r(i) = max(diag(L1));
    [K, cN, rN] = graphK(L1,TOL);
    if cN < COND_TOL
        continue;
    end
    i = i + 1;
end
%%% Obtaining the eigenvalue sequences using the eigenshuffle function from
%%% MATLAB FEX.
% keyboard;
[K1_seq, cN1, rN1] = graphKseq(L1_seq, TOL);
% keyboard
[P1_seq, R1_seq] = lyapPseq(L1_seq, K1_seq, Q, ALPHA);

[V1_seq, D1p_seq] = eigenshuffle(P1_seq);

[V1r_seq, D1r_seq] = eigenshuffle(R1_seq);
D1r_real = real(D1r_seq);
D1r_imag = imag(D1r_seq);
D1r_abs = abs(D1r_seq);

[V1k_seq, D1k_seq] = eigenshuffle(K1_seq);
D1k_real = real(D1k_seq);
D1k_imag = imag(D1k_seq);
D1k_abs = abs(D1k_seq);

R1_norm_seq = getNormSeq(R1_seq, 'fro');
K1_norm_seq = getNormSeq(K1_seq, 'fro');
P1_norm_seq = getNormSeq(P1_seq, 'fro');

%%% Plotting the variation sequence.
plot_colors = [cd.simulink_blue; cd.simulink_green; cd.simulink_red; cd.simulink_violet; cd.simulink_brown; cd.simulink_cyan; cd.old_default_black];

Pseq_dim = size(P1_seq);
Lp_max_ = zeros([1, Pseq_dim(end)]);
Lp_min_ = zeros([1, Pseq_dim(end)]);

[Vq, Dq] = eigenshuffle(Q);
Lq_max = max(Dq);
Lq_min = min(Dq);

Ubound_seq = zeros([1, Pseq_dim(end)]); % -Lq_min/Lp_max
Lbound_seq = zeros([1, Pseq_dim(end)]); % -Lq_max/Lp_min

for i=1:Pseq_dim(end)
    Di = D1p_seq(:,i);
    Di_max = max(Di);
    Di_min = min(Di);
    Lp_max_(i) = Di_max;
    Lp_min_(i) = Di_min;
    Ubound_seq(i) = -Lq_min/Di_max;
    Lbound_seq(i) = -Lq_max/Di_min;
end

for i=1:n
    eig_p = D1p_seq(i,:);

    figure(2);
    subplot(2,2,1)
    hold on
    plot(eig_p, 'LineWidth',1.5,'Color',plot_colors(i,:),'Marker','.');
    title('$|\lambda_P|$ v/s Iteration', 'Interpreter','latex');
    grid on;
    subplot(2,2,2);
    hold on;
    plot(L1_norm_seq, eig_p, 'LineWidth',1.5, 'Color', plot_colors(i,:),'Marker','.');
    title('$|\lambda_P|$ v/s $\|L\|_{fr}$','Interpreter','latex');
    grid on;
end

figure(2);

subplot(2,2,3)
plot(L1_norm_seq, 'LineWidth',1.5,'Color','k','Marker','.');
title('$\|L\|_{fr}$ v/s Iteration','Interpreter','latex');
grid on;
hold off;

subplot(2,2,4)
plot(Lp_max_, 'LineWidth', 1.5, 'LineStyle','-', 'Color', 'b', 'Marker', '.');
hold on;
plot(Lp_min_, 'LineWidth', 1.5, 'LineStyle','--', 'Color' ,'r', 'Marker', '.');
title('$\lambda_{p, max}, \lambda_{p, min}$ v/s Iteration','Interpreter','latex');
hold off;
grid on;

figure(4);

subplot(2,2,1)
plot(L1_norm_seq, Lp_max_, 'LineWidth', 1.5, 'LineStyle', '-', 'Color','b','Marker','.');
hold on;
plot(L1_norm_seq, Lp_min_, 'LineWidth', 1.5, 'LineStyle', '--', 'Color','r','Marker','.');
title('$\lambda_{p, max}, \lambda_{p, min}$ v/s $\|L_i\|_{fr}$','Interpreter','latex');
hold off;
grid on;

subplot(2,2,2)
plot(cN1, 'LineWidth', 1.5, 'Color','k','Marker','.', 'MarkerSize',8);
title('$cond(V_i)$ v/s Iteration (i)','Interpreter','latex');
grid on;

subplot(2,2,3)
plot(rN1, 'LineWidth', 1.5, 'Color','k','Marker','.', 'MarkerSize',8);
title('$rank(L_i)$ v/s Iteration (i)','Interpreter','latex');
grid on;

subplot(2,2,4)
plot(rN1, P1_norm_seq, 'LineWidth', 1.5, 'Color','k','Marker','.', 'MarkerSize',8);
title('$\|P_i\|_{fr}$ v/s $rank(L_i)$','Interpreter','latex');
grid on;

figure(3)
subplot(3,1,1)
plot(P1_norm_seq, 'LineWidth', 1.5, 'Color', 'b', 'Marker', '.');
title('$\|P_i\|_{fr}$ v/s Iteration (i)','Interpreter','latex');
grid on;

subplot(3,1,2)
plot(Ubound_seq, 'LineWidth', 1.5, 'LineStyle', '-', 'Color','b','Marker','.');
hold on;
plot(Lbound_seq, 'LineWidth', 1.5, 'LineStyle', '--', 'Color','r','Marker','.');
title('$-\frac{\lambda_{q, max}}{\lambda_{p, min}}, -\frac{\lambda_{q, min}}{\lambda_{p, max}}$ v/s Iteration (i)','Interpreter','latex');
hold off;
grid on;

subplot(3,1,3)
plot(L1_norm_seq, Ubound_seq, 'LineWidth', 1.5, 'LineStyle', '-', 'Color','b','Marker','.');
hold on;
plot(L1_norm_seq, Lbound_seq, 'LineWidth', 1.5, 'LineStyle', '--', 'Color','r','Marker','.');
title('$-\frac{\lambda_{q, max}}{\lambda_{p, min}}, -\frac{\lambda_{q, min}}{\lambda_{p, max}}$ v/s $\|L_i\|_{fr}$','Interpreter','latex');
hold off;
grid on;

figure(5)
subplot(2,2,1)
plot(L1_norm_seq, R1_norm_seq, 'LineWidth', 1.5, 'Color', 'g', 'Marker', '.');
title('$\|R_i\|_{fr}$ v/s $\|L_i\|_{fr}$','Interpreter','latex');
grid on;

subplot(2,2,2)
plot(L1_norm_seq, K1_norm_seq, 'LineWidth', 1.5, 'Color', 'b', 'Marker', '.');
title('$\|K_i\|_{fr}$ v/s $\|L_i\|_{fr}$','Interpreter','latex');
grid on;

subplot(2,2,3)
plot(L1_norm_seq, P1_norm_seq, 'LineWidth', 1.5, 'Color', 'r', 'Marker', '.');
title('$\|P_i\|_{fr}$ v/s $\|L_i\|_{fr}$','Interpreter','latex');
grid on;

subplot(2,2,4)
plot(R1_norm_seq, P1_norm_seq, 'LineWidth', 1.5, 'Color', 'k', 'Marker', '.');
title('$\|P_i\|_{fr}$ v/s $\|R_i\|_{fr}$','Interpreter','latex');
grid on;

figure(6)
plot(digraph(A));

% %% Epsilon bounded perturbation due to DoS attack. Analysis Type 1 - Continuous Variation.
% % TYPE 1: Randomly removing links one-by-one and checking eigenvalue
% % variation continuous by linear link removal with time at a sampling
% % frequency fs. Set fs according to the speed of the algorithm.
% fs = 100; % The step size is 1/fs.
% ts = 1/fs;
% A2 = A;
% G2 = G;
% nz = find((A > 0));
% num_nz = numel(nz);
% iterations = floor(sum(A2,"all")*fs);
% L2_seq = zeros(n,n,iterations);
% A2_seq = zeros(n,n,iterations);
% L2_seq(:,:,1) = L;
% A2_seq(:,:,1) = A;
% L2_norm_seq = zeros(1,iterations);
% L2_norm_seq(1) = norm(L, 'fro');
% g_disk_r = zeros(1, num_nz);
% j = 1;
% it_per_link = [1];
% for i = 1:num_nz
%     nz = find((A2 > 0));
%     nz = nz(randperm(numel(nz)));
%     rm_idx = nz(1);
%     % Continuous removal of link at fs.
%     while A2(rm_idx) > 0
%         A2(rm_idx) = A2(rm_idx) - ts;
%         A2(A2 < 0) = 0;
%         j = j + 1;
%         if j > iterations
%             break;
%         end
%         L2 = graph_laplacian(A2);
%         L2_seq(:,:,j) = L2;
%         A2_seq(:,:,j) = A2;
%         L2_norm_seq(j) = norm(L2, 'fro');
%     end
%     it_per_link = [it_per_link, j];
%     g_disk_r(i) = max(diag(L2));
% end
% 
% %%% Obtaining the eigenvalue sequences using the eigenshuffle function from
% %%% MATLAB FEX.
% [V2_seq, D2_seq] = eigenshuffle(L2_seq);
% D2_real = real(D2_seq);
% D2_imag = imag(D2_seq);
% D2_abs = abs(D2_seq);
% feidler_eigenvalue = zeros(1,iterations);
% %%% Extracting the feidler eigenvalue sequence at every timestep.
% for k=1:length(D2_abs)
%     vec = D2_abs(:,k);
%     feidler_eigenvalue(k) = min(setdiff(vec, min(vec)));
% end
% feidler_eigenvalue_interx = [];
% 
% %%% Plotting the variation sequence.
% plot_colors = [cd.simulink_blue; cd.simulink_green; cd.simulink_red; cd.simulink_violet; cd.simulink_brown; cd.simulink_cyan];
% for i=1:n
%     eig_p = D2_real(i,:);
%     eig_imag = D2_imag(i,:);
%     eig_abs = D2_abs(i,:);
%     d_eig_real = horzcat(diff(eig_p),0);
%     d_eig_imag = horzcat(diff(eig_imag),0);
%     % Calculating the change points in the feidler eigenvalue.
%     feidler_eigenvalue_interx = horzcat(feidler_eigenvalue_interx, InterX([feidler_eigenvalue;1:iterations], [eig_abs; 1:iterations]));
%     figure(4);
%     subplot(3,2,1)
%     hold on
%     plot(eig_p, eig_imag, 'LineWidth', 1, 'Color', plot_colors(i,:));
%     title('2D Plot of Eigenvalue Variation', 'Interpreter','latex');
%     grid on;
%     subplot(3,2,2)
%     hold on
%     plot(eig_abs, 'LineWidth',1,'Color',plot_colors(i,:));
%     title('$|\lambda|$ v/s Iteration', 'Interpreter','latex');
%     grid on;
%     subplot(3,2,3);
%     hold on;
%     plot(L2_norm_seq, eig_abs, 'LineWidth',1, 'Color', plot_colors(i,:));
%     title('$|\lambda|$ v/s $\|L\|_{fr}$','Interpreter','latex');
%     grid on;
%     figure(6)
%     subplot(1,2,1)
%     hold all;
%     plot(abs(eig_p), 'LineWidth', 1.5, 'Color', plot_colors(i,:));
%     title('$|Re(\lambda)|$ v/s Iteration', 'Interpreter','latex');
%     grid on;
%     subplot(1,2,2)
%     hold all;
%     plot(abs(eig_imag), 'LineWidth', 1.5, 'Color', plot_colors(i,:));
%     title('$|Im(\lambda)|$ v/s Iteration', 'Interpreter','latex');
%     grid on;
% end
% 
% figure(4);
% subplot(3,2,4)
% plot(L2_norm_seq, 'LineWidth',1,'Color','k');
% title('$\|L\|_{fr}$ v/s Iteration','Interpreter','latex');
% grid on;
% subplot(3,2,5)
% plot(horzcat(diff(L2_norm_seq),0), 'LineWidth',1,'Color','k');
% title('$\Delta\|L\|_{fr}$ v/s Iteration','Interpreter','latex');
% grid on;
% subplot(3,2,6)
% % plot(1:1:iterations, feidler_eigenvalue, feidler_eigenvalue_interx(1,:), feidler_eigenvalue_interx(2,:), 'bo', 'LineWidth',1);
% plot(feidler_eigenvalue,'b', 'LineWidth',1);
% title('$|\lambda_2|$ (Feidler Eigenvalue) v/s Iteration','Interpreter','latex');
% grid on;
% hold off;
% g_disk_r_max = max(g_disk_r);
% % Changing the default map for better visibility.
% choice = 1; % Using Jet. Uncomment the line below to toggle the menu.
% % choice = menu('Which ColorOrder do you want?', 'jet', 'random', 'hsv', 'hot', 'cool', 'spring', 'summer',...
% % 	'autumn', 'winter', 'lines', 'gray', 'bone', 'copper', 'pink');
% newColorMap = color_order(choice, num_nz);
% set(gca, 'ColorOrder', newColorMap, 'NextPlot', 'replacechildren','Colormap',newColorMap);
% 
% % Now get the new set of default plot colors.
% % Verify it changed by printing out the new default color set to the command window.
% % Uncomment the line below to check the map.
% % newColorOrder = get(gca,'ColorOrder');
% for i=1:num_nz
%     eig_p = D2_real(:,it_per_link(i):it_per_link(i+1)-1);
%     eig_imag = D2_imag(:,it_per_link(i):it_per_link(i+1)-1);
%     figure(5)
%     hold all;
%     plt = plot(eig_p(:), eig_imag(:), '.', 'MarkerSize',6);
%     colororder(newColorMap)
%     title('Variation in Largest Gersgorin Disk','Interpreter','latex');
%     if i > 1
%         if abs(g_disk_r(i) - g_disk_r(i-1)) < 1e-6
%             continue
%         end
%     end
%     viscircles([g_disk_r(i),0], g_disk_r(i), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', plt.Color);
%     xlim([0 Inf]);
%     axis equal;
% end
% % colororder(newColorMap)