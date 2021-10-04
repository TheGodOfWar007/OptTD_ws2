%% Mode Params
ENABLE_WINDOW_MODE = true;
ENABLE_GIF = true;

%% Unpacking Simulink logged data to prepare for plotting.
% Note that, to run this code you should have run the simulink simulation
% first. Otherwise this code will certainly not work.
logged_signals = out.logsout;
%%% Getting the signals.

% Actual States (Actual Trajectory)
s1_obj = logged_signals.getElement('s1');
s2_obj = logged_signals.getElement('s2');
s3_obj = logged_signals.getElement('s3');
s4_obj = logged_signals.getElement('s4');
s5_obj = logged_signals.getElement('s5');

s1 = s1_obj.Values.Data;
s2 = s2_obj.Values.Data;
s3 = s3_obj.Values.Data;
s4 = s4_obj.Values.Data;
s5 = s5_obj.Values.Data;

% Reference States (Ref Trajectory)
ps1_obj = logged_signals.getElement('ps1');
ps2_obj = logged_signals.getElement('ps2');
ps3_obj = logged_signals.getElement('ps3');
ps4_obj = logged_signals.getElement('ps4');
ps5_obj = logged_signals.getElement('ps5');

ps1 = ps1_obj.Values.Data;
ps2 = ps2_obj.Values.Data;
ps3 = ps3_obj.Values.Data;
ps4 = ps4_obj.Values.Data;
ps5 = ps5_obj.Values.Data;

%% Storing Data in a matrix to allow for easy addition of multiple bots when needed.
ps_dim = size(ps1);
s_dim = size(s1);
ps_mat = zeros(ps_dim(1), ps_dim(2), ps_dim(3), num_bots);
s_mat = zeros(s_dim(1), s_dim(2), s_dim(3), num_bots);

%%% Storing in the matrices.
ps_mat(:,:,:,1) = ps1;
ps_mat(:,:,:,2) = ps2;
ps_mat(:,:,:,3) = ps3;
ps_mat(:,:,:,4) = ps4;
ps_mat(:,:,:,5) = ps5;
s_mat(:,:,:,1) = s1;
s_mat(:,:,:,2) = s2;
s_mat(:,:,:,3) = s3;
s_mat(:,:,:,4) = s4;
s_mat(:,:,:,5) = s5;

%% Plot Color List RGB Triplets.
yellow = [1,1,0];
magenta = [1,0,1];
cyan = [0,1,1];
red = [1,0,0];
green = [0,1,0];
blue = [0,0,1];
white = [1,1,1];
black = [0,0,0];
simulink_blue = [0, 0.4470, 0.7410];
simulink_red = [0.8500, 0.3250, 0.0980];
simulink_yellow = [0.9290, 0.6940, 0.1250];
simulink_violet = [0.4940, 0.1840, 0.5560];
simulink_green = [0.4660, 0.6740, 0.1880];
simulink_cyan = [0.3010, 0.7450, 0.9330];
simulink_brown = [0.6350, 0.0780, 0.1840];
old_default_green = [0, 0.5, 0];
old_default_cyan = [0, 0.75, 0.75];
old_default_pink = [0.75, 0, 0.75];
old_default_olive_green = [0.75, 0.75, 0];
old_default_black = [0.25, 0.25, 0.25];
%% Plot Properties and time vector.
% Placeholders
x_idx = 1;
y_idx = 2;
col_vector = 1; % Just a systematic placeholder. 
% Time Vector
t = out.tout;
dt = vertcat(t(2:end), t(end)) - t; % time_steps.

%%% Gif save properties
filename_prefix = 'BotAnimPlot';

% All colour lists should have num_bots elements.
%%% Reference Trajectory.
ref_traj_line_style = '--';
ref_traj_linewidth = 1.5;
ref_traj_colors = [simulink_blue; simulink_green; simulink_red; old_default_black; old_default_pink];

%%% Actual Trajectory
act_traj_line_style = '-';
act_traj_linewidth = 1;
act_traj_colors = ref_traj_colors;

%%% Bot Marker
bot_marker_type = '^';
bot_marker_size = 6; % Default = 6

bot_marker_face_colors = [blue; green; red; black; magenta];

%%% Bot Marker Text Label
text_label_base_str = '\leftarrow B'; % Bot Serial number will be appended to this.

%%% Flags
anim_plot_enable_text = true;

%% Sliding Window Plot Properties.
enable_window_mode = ENABLE_WINDOW_MODE;

%%% Properties,
center_idx = central_bot_idx;
window_size = [4, 4]; %[x.size y.size]


%% Drawing the animated plot.
fig = figure(1);
figAxes = gca;
hold on

% [Note: Animated Plot Method was deleted because axis update method
% provided better control and memory management.]

% USING PLOT AXIS DATA UPDATE METHOD.
filename = '';
if enable_window_mode
    filename = strcat(filename_prefix, '_windowMode.gif');
else
    filename = strcat(filename_prefix, '.gif');
end

% Using num_bots zero matrices to assign num_bots objects.
ref_traj_plots = plot(zeros(num_bots));
hold on
act_traj_plots = plot(zeros(num_bots));
hold on
bot_markers = plot(zeros(num_bots));
hold on
text_array = [];
% Set Properties
for i=1:num_bots
    % Reference Trajectory Properties
    ref_traj_plots(i).Color = ref_traj_colors(i,:);
    ref_traj_plots(i).LineStyle = ref_traj_line_style;
    ref_traj_plots(i).LineWidth = ref_traj_linewidth;
    ref_traj_plots(i).DisplayName = strcat('BOT', string(j), ' Ref Traj');
    % Actual Trajectory Properties
    act_traj_plots(i).Color = act_traj_colors(i,:);
    act_traj_plots(i).LineStyle = act_traj_line_style;
    act_traj_plots(i).LineWidth = act_traj_linewidth;
    act_traj_plots(i).DisplayName = strcat('BOT', string(j), ' Act Traj');
    
    % Bot Marker Properties and first instance.
    bot_markers(i).Marker = bot_marker_type;
    bot_markers(i).MarkerSize = bot_marker_size;
    bot_markers(i).MarkerFaceColor = bot_marker_face_colors(i,:);
    bot_markers(i).XData = [s_mat(x_idx, 1, 1, i)];
    bot_markers(i).YData = [s_mat(y_idx, 1, 1, i)];
    bot_markers(i).DisplayName = strcat('BOT', string(j));
    hold on;
    
    % Set the display text objects.
    text_array = [text_array text()];
    text_array(i).String = strcat(text_label_base_str, string(i));
    
    % Emptying the data of the trajectory traces so that elements can be
    % added later while looping to save memory and avoid the squeeze
    % operation.
    ref_traj_plots(i).XData = [];
    ref_traj_plots(i).YData = [];
    act_traj_plots(i).XData = [];
    act_traj_plots(i).YData = [];
    hold on
    drawnow;
end
% Now plotting using the axis update method.
disp('####### ---- ANIMATED PLOTTER START ---- #######');
disp('----------------------------------------------------------------------');
disp('WARNING: Do not close the figure until it has finished plotting, the')
disp('gif will end at that point. Also do not resize the window otherwise ')
disp('the gif will have partial or zoomed out frames.');
for i=1:length(t)
    for j=1:num_bots
        % Updating Reference Trajectory Trace
        ref_traj_plots(j).XData = horzcat(ref_traj_plots(j).XData, ps_mat(x_idx, 1, i, j));
        ref_traj_plots(j).YData = horzcat(ref_traj_plots(j).YData, ps_mat(y_idx, 1, i, j));
        hold on
        % Updating Actual Trajectory Trace
        act_traj_plots(j).XData = horzcat(act_traj_plots(j).XData, s_mat(x_idx, 1, i, j));
        act_traj_plots(j).YData = horzcat(act_traj_plots(j).YData, s_mat(y_idx, 1, i, j));
        hold on
        % Updating bot marker and text position.
        bot_markers(j).XData = [s_mat(x_idx, 1, i, j)];
        bot_markers(j).YData = [s_mat(y_idx, 1, i, j)];
        if anim_plot_enable_text
            text_array(j).Position = [s_mat(x_idx, 1, i, j), s_mat(y_idx, 1, i, j)];
        end
    end
    
    if enable_window_mode
       center_x = s_mat(x_idx, 1, i, center_idx);
       center_y = s_mat(y_idx, 1, i, center_idx);
       figAxes.XLim = [center_x-window_size(1)/2, center_x+window_size(1)/2];
       figAxes.YLim = [center_y-window_size(2)/2, center_y+window_size(2)/2];
    end
    
    % Making/Updating the gif.
    if ENABLE_GIF
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im,256);

        if i == 1 
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append'); 
        end 
    end
    
    drawnow;
    pause(dt(i));
end
disp('----------------------------------------------------------------------');
disp('The Animated plot has finished plotting. You can close the figure now.');

hold off
