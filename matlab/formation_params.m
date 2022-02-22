%% Preliminaries
num_bots = 5;
%% Perturbation in Initial Positions.
p_mu = 1;
p_sigma = 1;
init_position_perturbation = p_sigma*randn(num_bots, 3) + p_mu;

%% Formation Config

init_positions = [4, -1, 0; 
                  6, -1, 0;
                  6, 1, 0;
                  4, 1, 0;
                  5, 0, 0];
              
init_positions = init_positions + init_position_perturbation;
              
init_orientations = [0, 0, 0;
                     0, 0, 0;
                     0, 0, 0;
                     0, 0, 0;
                     0, 0, 0];

%% Formation Graph
%%% An Arbitrary Example
% A = [0, 0.9, 0, 0.9, 1;
%      0, 0, 0.9, 0, 0.8;
%      0.9, 0.9, 0, 0, 0.7;
%      0.6, 0, 0.5, 0, 0.9;
%      0.9, 0.9, 0.9, 0.6, 0];

%%% With node 5 as an authority
% A = [0, 0.9, 0, 0.9, 0;
%      0, 0, 0.9, 0, 0;
%      0.9, 0.9, 0, 0, 0;
%      0.6, 0, 0.5, 0, 0;
%      0.9, 0.9, 0.9, 0.6, 0];

%%% With a weakly connected graph.
% A = [0, 0.9, 0 ,0, 0;
%      0, 0, 0.3, 0, 0;
%      0.5, 0.7, 0, 0, 0;
%      0, 0, 0.3, 0, 0.9;
%      0, 0, 0, 0.6, 0];

%%% Disconnected DiGraph
A = [0, 0.9, 0 ,0, 0;
     0, 0, 0.3, 0, 0;
     0.5, 0.7, 0, 0, 0;
     0, 0, 0, 0, 0.9;
     0, 0, 0, 0.6, 0];
A = A - diag(diag(A)); % Removing self connections.
L = diag(sum(A,2)) - A;
G = digraph(A);
plot(G);
%% Test Motion Planner
static_formation = [0.0, 0.0, 0.0;
                    2.0, 0.0, 0.0;
                    2.0, 2.0, 0.0;
                    0.0, 2.0, 0.0;
                    1.0, 1.0, 0.0];
central_bot_idx = 5;
follower_bot_idx = [1,2,3,4];
nf = size(follower_bot_idx);
% Storing transforms for the trajectory generation.
tf_list = zeros(4,4,nf(2));
for i=1:1:nf(2)
    j = follower_bot_idx(i);
    cbi = central_bot_idx;
    trvec = static_formation(j,:) - static_formation(cbi,:);
    tf_list(:,:,i) = trvec2tform(trvec);
end

% Waypoints for a straight line.
% t_step = 2;
% t_s = 5 + t_step;
% t_e = 40;
% V = 1; % m/s
% wp_straight_line = [];
% tp_straight_line = [];
% v_straight_line = [];
% p_t = [4.7; 1];
t_step = 2;
t_s = 0 + t_step;
t_e = 40;
V = 1; % m/s
wp_straight_line = [];
tp_straight_line = [];
v_straight_line = [];
p_t = [4.0; 0];
theta = 45*pi/180; % rad
for t=t_s:t_step:t_e
    tp_straight_line = horzcat(tp_straight_line, t);
    p_t(1) = p_t(1) + V*t_step*cos(theta);
    p_t(2) = p_t(2) + V*t_step*sin(theta);
    v_straight_line = horzcat(v_straight_line, [V*cos(theta); V*sin(theta); 0]);
    wp_straight_line = horzcat(wp_straight_line, vertcat(p_t, 0));
end

% waypoints = horzcat([4;0;0], [4.35; 0.4; 0], [4.7; 1; 0], wp_straight_line);
% time_points = horzcat([0, 2.5, 5], tp_straight_line);
% velocity_conditions = horzcat([0;0;0],[cos(30*pi/180); sin(30*pi/180); 0],[cos(45*pi/180); sin(45*pi/180); 0],v_straight_line);
waypoints = horzcat([4;0;0], wp_straight_line);
time_points = horzcat(0, tp_straight_line);
velocity_conditions = horzcat([0;0;0], v_straight_line);

