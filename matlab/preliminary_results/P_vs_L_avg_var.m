% This file was created for the sake of analyzing the effect of the link
% distribution in L on the degradation of the performance, as discussed in
% section 1.5 of the 21/11/2021 update of the PDF. Please refer to the
% final version of the Thesis for the same. Note that the measures
% presented in this file are mostly probabilistic and statistical in nature
% and the deterministic nature of the digraph can't be discussed with the
% increased of generality and lack of research in the field of general
% directed graphs.
close all;
%% IMPORTANT PARAMETERS
NUM_NODES = 7;
TOL = 1e-6; % For float comparisons.
COND_TOL = 1; % Skip if condition number < COND_TOL.
% Number of initial connections
MAX_LINKS = NUM_NODES*(NUM_NODES - 1); % Exclude Self Connections
M_MIN = ceil((0.6)*MAX_LINKS); % Min % Initial Connections.
M_MAX = ceil((1)*MAX_LINKS); % Max % Initial Connections.
M_STEP = 1; % Step size for link number increase.
% Range of uniform distribution to choose provisional weights from.
W_MIN = 0.5;
W_MAX = 1.5;
% Range of the norms of the Adj matrix to analyse.
A_F_MIN = 5;
A_F_MAX = 20;
% Hyperparameters
A_F_STEP = 0.6;
NSP_STEP = 1000;

% Attack Params
A_FRAC = 0.7; % Increase this to increase attack strength.
E_NORM = 5; 
ATTACK_MODE = 0;

% SIMULATION MODE FLAGS
ENABLE_SYMMETRY = true;
LINK_DISABLING = true;

% Lyapunov Equation
kQ = 4;
ALPHA = 4;

% Initializing Preliminary variables.
FR_VEC = A_F_MIN:A_F_STEP:A_F_MAX;
M_VEC = M_MIN:M_STEP:M_MAX;

% Aliases
n = NUM_NODES;
e_norm = E_NORM;

cd = color_dict();
%% P1-REDUCTION IN CONVERGENCE RATE BOUNDS WITH NUMBER OF NON-ZERO ELEMENTS
% CONSTANTS: ||A||_F & ||E||_F
%P1 Vars
close all;
A_NORM = 20;
e_norm = A_FRAC*A_NORM;

Dzeta_u_avg = zeros([1, numel(M_VEC)]);
Dzeta_l_avg = zeros([1, numel(M_VEC)]);
rt_avg = Dzeta_l_avg;
r_avg = Dzeta_l_avg;
zeta_u_i_avg = r_avg;
zeta_l_i_avg = r_avg;
zeta_u_f_avg = r_avg;
zeta_l_f_avg = r_avg;

d_zeta_u = zeros([1, NSP_STEP]);
d_zeta_l = d_zeta_u;
zeta_u_i = d_zeta_u;
zeta_l_i = d_zeta_l;
zeta_u_f = d_zeta_u;
zeta_l_f = d_zeta_l;

p_norm_seq = zeros([numel(M_VEC), NSP_STEP]);
lambda2_seq = p_norm_seq;
kap_seq = zeros([1, numel(M_VEC)]);
eig_dist_seq = kap_seq;
eig_dist_i = d_zeta_u;
eig_dist_seq_mat = p_norm_seq;
ubound_act_diff_seq = kap_seq;
ubound_act_diff_i = d_zeta_u;
a_norm_seq = kap_seq;
a_norm_i = d_zeta_u;
e_norm_seq = kap_seq;
e_norm_i = d_zeta_u;
% Norm sequence of the laplacian.
l_norm_seq = kap_seq;
l_norm_i = d_zeta_u;
% Rank lost after attack.
rank_lost_seq = kap_seq;

eig_L_avg = zeros([NUM_NODES, numel(M_VEC)]);
eig_L = zeros([NUM_NODES, NSP_STEP]);
eig_P_avg = eig_L_avg;
eig_P_i = eig_L;

Q = diag(repelem(kQ, NUM_NODES));
[Vq, Dq] = eigenshuffle(Q);
Lq_max_ = max(Dq);
Lq_min_ = min(Dq);
progressbar('Non-Zero Links','Steps per Iteration');
for m = 1:numel(M_VEC)
    rt_m = 0;
    r_m = 0;
    i = 1;
    kap_max = 0;
    while i < NSP_STEP
        A = nrmUMatrix(NUM_NODES, M_VEC(m), A_NORM,"SELF_CONN",false,"w_min",W_MIN, "w_max",W_MAX, "SYMMETRIC",ENABLE_SYMMETRY);
        % PRE-ATTACK
        L = graph_laplacian(A);
        [V_L, D_L] = eigenshuffle(L);
        eig_L(:,i) = D_L;
        [K, cN, rN] = graphK(L,TOL);
        R = L + ALPHA*K;
        a_norm_i(i) = norm(A, 'fro');
        l_norm_i(i) = norm(L, 'fro');
        if cN < COND_TOL
            continue
        end
        try
            P = lyap(-R, Q);
        catch
            continue
        end
        [Vp, Dp] = eigenshuffle(P);
        eig_P_i(:,i) = Dp;
        Lp_max_i = max(Dp);
        Lp_min_i = min(Dp);
        feidler_eig = min(setdiff(D_L, min(D_L)));
        lambda2_seq(m, i) = feidler_eig;
        p_norm_seq(m, i) = norm(P, 'fro');
        % ATTACK
        if LINK_DISABLING
            E = nrmLDAttack(A, e_norm, "SYMMETRIC",ENABLE_SYMMETRY);
        else
            E = nrmUAttack(A, e_norm, "SYMMETRIC",ENABLE_SYMMETRY);
        end
        At = A - E;
        Lt = graph_laplacian(At);
        % Post Attack
        [Kt, cNt, rNt] = graphK(Lt,TOL);
        if cNt < COND_TOL
            continue
        end
        Rt = Lt + ALPHA*Kt;
        Del_R = Rt - R;
        Del_R_nrm = norm(Del_R, 'fro');
        R_nrm = norm(R, 'fro');
        kappa = Del_R_nrm/R_nrm;
        if kappa > kap_max
            kap_max = kappa;
        end
        try
            Pt = lyap(-Rt, Q);
        catch
            continue
        end
        [Vpt, Dpt] = eigenshuffle(Pt);
        Lp_max_f = max(Dpt);
        Lp_min_f = min(Dpt);
        % Change in Bounds
        z_u_i = -Lq_min_/Lp_max_i;
        z_u_f = -Lq_min_/Lp_max_f;
        z_l_i = -Lq_max_/Lp_min_i;
        z_l_f = -Lq_max_/Lp_min_f;
        zeta_l_i(i) = z_l_i;
        zeta_u_i(i) = z_u_i;
        zeta_l_f(i) = z_l_f;
        zeta_u_f(i) = z_u_f;

        d_zeta_u(i) = z_u_f - z_u_i;
        d_zeta_l(i) = z_l_f - z_l_i;
        r_m = r_m + rN/NSP_STEP;
        rt_m = rt_m + rNt/NSP_STEP;

        pt_p_eig_dist = max(max(pdist2(horzcat(Dpt, zeros([numel(Dpt),1])), horzcat(Dpt, zeros([numel(Dp),1])))));
        eig_dist_i(i) = pt_p_eig_dist;
        sep_term = sqrt(3)*(kappa^2)/feidler_eig;
        % A conservative kappa taken for this one.
        eig_dist_ubound = 4*2^(2 - 1/n)*p_norm_seq(m,i)*((1 + sep_term)^(1 - 1/n))*(sep_term^(1/n));
        ubound_act_diff_i(i) = eig_dist_ubound - pt_p_eig_dist;
        eig_dist_seq_mat(m,i) = pt_p_eig_dist;
        e_norm_i(i) = norm(E, 'fro');

        rank_lost_seq(m) = rank_lost_seq(m) + (rank(L) - rank(Lt));

        progressbar([],i/NSP_STEP);
        i = i + 1;
    end
    Dzeta_l_avg(m) = mean(d_zeta_l);
    Dzeta_u_avg(m) = mean(d_zeta_u);
    r_avg(m) = r_m;
    rt_avg(m) = rt_m;
    zeta_u_i_avg(m) = mean(zeta_u_i);
    zeta_l_i_avg(m) = mean(zeta_l_i);
    zeta_u_f_avg(m) = mean(zeta_u_f);
    zeta_l_f_avg(m) = mean(zeta_l_f);
    eig_L_avg(:,m) = mean(eig_L, 2);
    eig_P_avg(:,m) = mean(eig_P_i, 2);
    kap_seq(m) = kap_max;
    eig_dist_seq(m) = mean(eig_dist_i);
    ubound_act_diff_seq(m) = mean(ubound_act_diff_i);
    a_norm_seq(m) = mean(a_norm_i);
    e_norm_seq(m) = mean(e_norm_i);
    l_norm_seq(m) = mean(l_norm_i);

    rank_lost_seq(m) = rank_lost_seq(m)/NSP_STEP;

    progressbar(m/numel(M_VEC),[]);
    progressbar([],0);
end
progressbar(1);

symm_upb_seq = zeros([1, numel(M_VEC)]);
symm_upb_i = zeros([1, NSP_STEP]);
ubound_act_diff_seq2 = symm_upb_seq;
ubound_act_diff_i2 = symm_upb_i;
progressbar('Outer Filter Loop', 'Steps Per Iteration');
for m=1:numel(M_VEC)
    i = 1;
    while i < NSP_STEP
        sep_term = sqrt(3)*(kap_seq(m)^2)/lambda2_seq(m,i);
        symm_upb_i(i) = 4*2^(2 - 1/n)*p_norm_seq(m,i)*((1 + sep_term)^(1 - 1/n))*(sep_term^(1/n));
        ubound_act_diff_i2(i) = symm_upb_i(i) - eig_dist_seq_mat(m,i);
        progressbar([],i/NSP_STEP);
        i = i + 1;
    end
    symm_upb_seq(m) = mean(symm_upb_i);
    ubound_act_diff_seq2(m) = mean(ubound_act_diff_i2);
    progressbar(m/numel(M_VEC),[]);
    progressbar([],0);
end
progressbar(1);

figure(2)

subplot(2,2,1)
plot(M_VEC, zeta_l_i_avg,"LineStyle","--","Color",cd.black,"LineWidth",1.5);
hold on
plot(M_VEC, zeta_u_i_avg,"LineStyle","-","Color",cd.blue,"LineWidth",1.5);
hold off
xlabel("No. of Non-Zero Links")
ylabel("$\zeta_{l}, \zeta_{u}$", "Interpreter","latex", "FontWeight","bold");
title("Convergence Rate Bounds w/o Attack","Interpreter","latex");
grid off

subplot(2,2,2)
plot(M_VEC, zeta_l_f_avg,"LineStyle","--","Color",cd.black,"LineWidth",1.5);
hold on
plot(M_VEC, zeta_u_f_avg,"LineStyle","-","Color",cd.blue,"LineWidth",1.5);
hold off
xlabel("No. of Non-Zero Links")
ylabel("$\tilde{\zeta}_{l}, \tilde{\zeta}_{u}$", "Interpreter","latex", "FontWeight","bold");
title("Convergence Rate Bounds under Attack","Interpreter","latex");
grid off

subplot(2,2,3)
plot(M_VEC, r_avg,"LineStyle","--","Color",cd.magenta,"LineWidth",1.5);
hold on
plot(M_VEC, rt_avg,"LineStyle","-","Color",cd.green,"LineWidth",1.5);
hold off
title("$rank(L), rank(\tilde{L})$ v/s No. of non-zero links","Interpreter","latex", "FontWeight","bold");
grid on

subplot(2,2,4)
plot(M_VEC, abs(eig_L_avg),"LineWidth",1,"Marker",".");
title("$|eig(L)|$ v/s N(non-zero links)", "Interpreter","latex");
grid on

figure(3)
subplot(2,2,1)
plot(M_VEC, Dzeta_u_avg, "LineStyle","-","Color",cd.simulink_blue,"LineWidth",1.5,"Marker",".");
title("$-\Delta \zeta_u$ v/s No. of Non-zero Links", "Interpreter","latex");
grid on
subplot(2,2,2)
plot(M_VEC, Dzeta_l_avg, "LineStyle","--","Color",cd.simulink_green,"LineWidth",1.5,"Marker",".");
title("$-\Delta \zeta_l$ v/s No. of Non-zero Links", "Interpreter","latex");
grid on
subplot(2,2,3)
plot(M_VEC, mean(lambda2_seq,2)',"Color",cd.old_default_black)
title("$\lambda_2$ v/s Non-Zero Links","Interpreter","latex")
grid on
subplot(2,2,4)
plot(M_VEC, kap_seq, "Color", cd.simulink_violet)
title("$\kappa_{max}$ v/s Non-Zero Links","Interpreter","latex")
grid on

figure(4)
subplot(2,1,1)
plot(M_VEC, symm_upb_seq, "Color", cd.simulink_blue, "LineWidth",1.5)
hold on
plot(M_VEC, eig_dist_seq, "Color", cd.simulink_red, "LineStyle","--", "LineWidth",1.5)
title(" Avg UBound on $|\tilde{\lambda}_P - \lambda_P|$ v/s Non-Zero Links (Symm Case)","Interpreter","latex")
grid on
subplot(2,1,2)
plot(M_VEC, eig_dist_seq, "Color", cd.simulink_red, "LineStyle","-", "LineWidth",1.5)
title(" Avg $|\tilde{\lambda}_P - \lambda_P|$ v/s Non-Zero Links (Symm Case)","Interpreter","latex")
grid on

figure(5)
subplot(2,1,1)
plot(M_VEC, ubound_act_diff_seq2, "Color", cd.blue, "LineWidth",1.5)
title("Difference b/w upper bound and actual")
frac = 1;
% sub_vec_ = ubound_act_diff_seq(ceil(frac*numel(M_VEC)/2):end);
% ylim([min(sub_vec_), max(sub_vec_)])
% xlim([1, NUM_NODES^2])
grid on
subplot(2,1,2)
plot(M_VEC, mean(p_norm_seq, 2)',"Color","k", "LineWidth",1.5)
title("Mean $\|P\|_F$ v/s Non-Zero Links", "Interpreter","latex")
grid on

figure(6)
subplot(2,1,1)
plot(M_VEC, eig_P_avg,"LineWidth",1,"Marker",".");
title("$|eig(P)|$ v/s N(non-zero links)", "Interpreter","latex");
grid on
subplot(2,1,2)
plot(M_VEC, l_norm_seq, "LineWidth",1.5)
title("Mean $\|L\|_F$ v/s Non-Zero Links", "Interpreter","latex")
grid on

figure(7)
plot(M_VEC, rank_lost_seq, "LineWidth", 1.5, "Marker",".");
title("Post-Attack Rank Loss", "Interpreter","latex");
xlabel("No. of Non-zero links");
ylabel("rank($L$) - rank($\tilde{L}$)", "Interpreter","latex");
grid on;

figure(8)
plot(M_VEC, a_norm_seq, "Color", "r", "LineWidth",1.5)
hold on
plot(M_VEC, e_norm_seq, "Color","b", "LineStyle","--", "LineWidth",1.5)
title("Verifying $\|A\|_F$", "Interpreter","latex")
grid on
