%% Spider plots
% loads decoding results, e.g. PFC_validationAccuracy_WT is a matrix with
% 100 rows (cross-validation repititions) and 10 columns (decoding
% performance parameters returned by ft_subspacediscriminant)
% Spider plots created with the following code:
% Víctor Martínez-Cagigal (2021). Polygonal Plot (https://www.mathworks.com/matlabcentral/fileexchange/62200-polygonal-plot), MATLAB Central File Exchange. Retrieved April 27, 2021.
%%
clear all
close all
%% DMTP %%
load 'F:\MD\results\DMTP\1sSD_5sdelay\Classifier\DMTP_local'
operant_all(1,1,:) = PFC_validationAccuracy_WT(:,1);
operant_all(2,1,:) = dHC_validationAccuracy_WT(:,1);
operant_all(3,1,:) = vHC_validationAccuracy_WT(:,1);
operant_all(4,1,:) = MD_validationAccuracy_WT(:,1);

operant_all(1,3,:) = PFC_validationAccuracy_perm_WT(:,1);
operant_all(2,3,:) = dHC_validationAccuracy_perm_WT(:,1);
operant_all(3,3,:) = vHC_validationAccuracy_perm_WT(:,1);
operant_all(4,3,:) = MD_validationAccuracy_perm_WT(:,1);

load 'F:\MD\results\DMTP\1sSD_5sdelay\Classifier\DMTP_conn'
operant_all(5,1,:) = PFC_dHC_validationAccuracy_WT(:,1);
operant_all(6,1,:) = PFC_vHC_validationAccuracy_WT(:,1);
operant_all(7,1,:) = PFC_MD_validationAccuracy_WT(:,1);
operant_all(8,1,:) = vHC_dHC_validationAccuracy_WT(:,1);

operant_all(5,3,:) = PFC_dHC_validationAccuracy_perm_WT(:,1);
operant_all(6,3,:) = PFC_vHC_validationAccuracy_perm_WT(:,1);
operant_all(7,3,:) = PFC_MD_validationAccuracy_perm_WT(:,1);
operant_all(8,3,:) = vHC_dHC_validationAccuracy_perm_WT(:,1);

%% DNMTP %%
% clear all
load 'F:\MD\results\DNMTP\DNMTP_local'
operant_all(1,2,:) = PFC_validationAccuracy_WT(:,1);
operant_all(2,2,:) = dHC_validationAccuracy_WT(:,1);
operant_all(3,2,:) = vHC_validationAccuracy_WT(:,1);
operant_all(4,2,:) = MD_validationAccuracy_WT(:,1);

operant_all(1,4,:) = PFC_validationAccuracy_perm_WT(:,1);
operant_all(2,4,:) = dHC_validationAccuracy_perm_WT(:,1);
operant_all(3,4,:) = vHC_validationAccuracy_perm_WT(:,1);
operant_all(4,4,:) = MD_validationAccuracy_perm_WT(:,1);

load 'F:\MD\results\DNMTP\DNMTP_conn'
operant_all(5,2,:) = PFC_dHC_validationAccuracy_WT(:,1);
operant_all(6,2,:) = PFC_vHC_validationAccuracy_WT(:,1);
operant_all(7,2,:) = PFC_MD_validationAccuracy_WT(:,1);
operant_all(8,2,:) = vHC_dHC_validationAccuracy_WT(:,1);

operant_all(5,4,:) = PFC_dHC_validationAccuracy_perm_WT(:,1);
operant_all(6,4,:) = PFC_vHC_validationAccuracy_perm_WT(:,1);
operant_all(7,4,:) = PFC_MD_validationAccuracy_perm_WT(:,1);
operant_all(8,4,:) = vHC_dHC_validationAccuracy_perm_WT(:,1);

opt_area.err        = 'sem';
opt_area.FaceAlpha  = 0.3;
opt_area.Color      = ['r';'g';'r';'g'];
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_lines.Labels    = false;
opt_lines.Legend    = {'DMTP','DNMTP'};
opt_lines.Color = ['r';'g';'r';'g'];
                  
opt_axes.Background = 'w';
minvalue = 0.4;
maxvalue = 1;
opt_axes.Ticks = minvalue:0.1:maxvalue;

opt_axes.Labels     = {'Local PFC','Local dHC','Local vHC','Local MD',...
    'PFC-\newlinedHC','PFC-vHC','PFC-MD','vHC-dHC'};

figure
% subplot(1,3,2)
polygonplot(operant_all,opt_axes,opt_lines,opt_area,minvalue,maxvalue);
% sgtitle('DNMTP')






%% T-Maze
% clear all
load 'F:\MD\results\T_Maze\Classifier\T_Maze_local'
T_Maze_all(1,1,:) = PFC_validationAccuracy_WT(:,1);
T_Maze_all(2,1,:) = dHC_validationAccuracy_WT(:,1);
T_Maze_all(3,1,:) = vHC_validationAccuracy_WT(:,1);
T_Maze_all(4,1,:) = MD_validationAccuracy_WT(:,1);

T_Maze_all(1,3,:) = PFC_validationAccuracy_perm_WT(:,1);
T_Maze_all(2,3,:) = dHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(3,3,:) = vHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(4,3,:) = MD_validationAccuracy_perm_WT(:,1);

load 'F:\MD\results\T_Maze\Classifier\T_Maze_conn'
T_Maze_all(5,1,:) = PFC_dHC_validationAccuracy_WT(:,1);
T_Maze_all(6,1,:) = PFC_vHC_validationAccuracy_WT(:,1);
T_Maze_all(7,1,:) = PFC_MD_validationAccuracy_WT(:,1);
T_Maze_all(8,1,:) = vHC_dHC_validationAccuracy_WT(:,1);

T_Maze_all(5,3,:) = PFC_dHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(6,3,:) = PFC_vHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(7,3,:) = PFC_MD_validationAccuracy_perm_WT(:,1);
T_Maze_all(8,3,:) = vHC_dHC_validationAccuracy_perm_WT(:,1);

load 'F:\MD\results\T_Maze\Classifier\T_Maze_local_30s'
T_Maze_all(1,2,:) = PFC_validationAccuracy_WT(:,1);
T_Maze_all(2,2,:) = dHC_validationAccuracy_WT(:,1);
T_Maze_all(3,2,:) = vHC_validationAccuracy_WT(:,1);
T_Maze_all(4,2,:) = MD_validationAccuracy_WT(:,1);

T_Maze_all(1,4,:) = PFC_validationAccuracy_perm_WT(:,1);
T_Maze_all(2,4,:) = dHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(3,4,:) = vHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(4,4,:) = MD_validationAccuracy_perm_WT(:,1);

load 'F:\MD\results\T_Maze\Classifier\30s_conn'
T_Maze_all(5,2,:) = PFC_dHC_validationAccuracy_WT(:,1);
T_Maze_all(6,2,:) = PFC_vHC_validationAccuracy_WT(:,1);
T_Maze_all(7,2,:) = PFC_MD_validationAccuracy_WT(:,1);
T_Maze_all(8,2,:) = vHC_dHC_validationAccuracy_WT(:,1);

T_Maze_all(5,4,:) = PFC_dHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(6,4,:) = PFC_vHC_validationAccuracy_perm_WT(:,1);
T_Maze_all(7,4,:) = PFC_MD_validationAccuracy_perm_WT(:,1);
T_Maze_all(8,4,:) = vHC_dHC_validationAccuracy_perm_WT(:,1);




opt_area.err        = 'sem';
opt_area.FaceAlpha  = 0.3;
opt_area.Color      = ['k';'b';'k';'b'];
opt_lines.LineWidth = 2;
opt_lines.LineStyle = '-';
opt_lines.Marker    = 'none';
opt_lines.Labels    = false;
opt_lines.Legend    = {'5s','30s'};
opt_lines.Color = ['k';'b';'k';'b'];
                  
opt_axes.Background = 'w';
minvalue = 0.4;
maxvalue = 1;
opt_axes.Ticks = minvalue:0.1:maxvalue;

opt_axes.Labels     = {'Local PFC','Local dHC','Local vHC','Local MD',...
    'PFC-\newlinedHC','PFC-vHC','PFC-MD','vHC-dHC'};

figure
% subplot(1,3,3)
polygonplot(T_Maze_all,opt_axes,opt_lines,opt_area,minvalue,maxvalue);
% title('T-Maze')


