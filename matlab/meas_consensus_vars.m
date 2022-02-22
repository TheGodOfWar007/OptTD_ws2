%% Important Global Parameters
FLOAT_ERROR_TOLERANCE = 1e-3;
%% Fetching logged data from Simulink
% Note that, to run this code you should have run the simulink simulation
% first. Otherwise this code will certainly not work.
logged_signals = out.logsout;
% Time Vector
t = out.tout;
dt = vertcat(t(2:end), t(end)) - t; % time_steps.

%% Measuring the agreement time through Group Disagreement Vector.
dis_ag_obj = logged_signals.getElement('dis_ag');
dis_ag = dis_ag_obj.Values.Data;
msk = dis_ag < FLOAT_ERROR_TOLERANCE;
t_agreed = t(msk);
agreement_time = t_agreed(1);
disp('----------------------------------------------------------------------');
disp_str = sprintf("Consensus Achieved in %.7f secs.", agreement_time);
disp(disp_str);
disp('----------------------------------------------------------------------');