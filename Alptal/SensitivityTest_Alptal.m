clear
close all

load LWsub_Alptal_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_AlbSurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Alptal_Fdiff.mat
LW_in_bc_CLM_Fdiff = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_Fdiff = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Alptal_SWC.mat
LW_in_bc_CLM_SWC = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_SWC = LW_in_bc_SP;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Alptal.mat','LWR_Atmosphere_Alptal_2004',...
    'LWR_Atmosphere_Alptal_2005','LWR_Atmosphere_Alptal_2006',...
    'LWR_Atmosphere_Alptal_2007')


%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
CandidCoral = [1 0.44 0.32];
EcstaticEmerald = [0.17 0.52 0.5];


%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
simperiod_2004 = 1:2346;    % from DataPrep_Alptal.m
start_2004 = 121*24-length(simperiod_2004)+1;    % from DataPrep_Alptal.m

% evaluation period
EP_2004 = start_2004+18:1728;   % 48 days
EP_2005 = 25:1752;              % 72 days
EP_2006 = 25:1872;              % 77 days
EP_2007 = 25:2256;              % 93 days

% gaps
    % 2004
gap_2004_1 = 601:624;           % 25 January 2004
gap_2004_2 = 1009:1032;         % 11 February 2004
gap_2004_3 = 1201:1224;         % 19 February 2004
gap_2004_4 = 1633:1728;         % 8 - 11 March 2004
    % 2005
gap_2005_1 = 553:624;           % 23 - 25 January 2005
gap_2005_2 = 769:792;           % 1 February 2005
gap_2005_3 = 1057:1272;         % 13 - 21 February 2005
gap_2005_4 = 1609:1656;         % 8 - 9 March 2005
    % 2006
gap_2006 = 937:1032;            % 8 - 11 February 2006           
    % 2007
gap_2007_1 = 49:96;             % 2 - 3 January 2007
gap_2007_2 = 1873:2016;         % 19 - 24 March 2007

% validation data
LW_sub_val_eval04 = vertcat(LW_in_bc_1h_all(EP_2004(1):gap_2004_1(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_4(end)+1:EP_2004(end),1));
LW_sub_val_eval05 = vertcat(LW_in_bc_1h_all(EP_2005(1):gap_2005_1(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_1(end)+1:gap_2005_2(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_2(end)+1:gap_2005_3(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_3(end)+1:gap_2005_4(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_4(end)+1:EP_2005(end),2));
LW_sub_val_eval06 = vertcat(LW_in_bc_1h_all(EP_2006(1):gap_2006(1)-1,3),...
    LW_in_bc_1h_all(gap_2006(end)+1:EP_2006(end),3));
LW_sub_val_eval07 = vertcat(LW_in_bc_1h_all(EP_2007(1):gap_2007_1(1)-1,4),...
    LW_in_bc_1h_all(gap_2007_1(end)+1:gap_2007_2(1)-1,4),...
    LW_in_bc_1h_all(gap_2007_2(end)+1:EP_2007(end),4));

% surface albedo
    % CLM
LW_sub_CLM_eval_AlbSurf04 = vertcat(LW_in_bc_CLM_AlbSurf(EP_2004(1):gap_2004_1(1)-1,1,:),...
    LW_in_bc_CLM_AlbSurf(gap_2004_1(end)+1:gap_2004_2(1)-1,1,:),...
    LW_in_bc_CLM_AlbSurf(gap_2004_2(end)+1:gap_2004_3(1)-1,1,:),...
    LW_in_bc_CLM_AlbSurf(gap_2004_3(end)+1:gap_2004_4(1)-1,1,:),...
    LW_in_bc_CLM_AlbSurf(gap_2004_4(end)+1:EP_2004(end),1,:));
LW_sub_CLM_eval_AlbSurf04 = squeeze(LW_sub_CLM_eval_AlbSurf04);
LW_sub_CLM_eval_AlbSurf05 = vertcat(LW_in_bc_CLM_AlbSurf(EP_2005(1):gap_2005_1(1)-1,2,:),...
    LW_in_bc_CLM_AlbSurf(gap_2005_1(end)+1:gap_2005_2(1)-1,2,:),...
    LW_in_bc_CLM_AlbSurf(gap_2005_2(end)+1:gap_2005_3(1)-1,2,:),...
    LW_in_bc_CLM_AlbSurf(gap_2005_3(end)+1:gap_2005_4(1)-1,2,:),...
    LW_in_bc_CLM_AlbSurf(gap_2005_4(end)+1:EP_2005(end),2,:));
LW_sub_CLM_eval_AlbSurf05 = squeeze(LW_sub_CLM_eval_AlbSurf05);
LW_sub_CLM_eval_AlbSurf06 = vertcat(LW_in_bc_CLM_AlbSurf(EP_2006(1):gap_2006(1)-1,3,:),...
    LW_in_bc_CLM_AlbSurf(gap_2006(end)+1:EP_2006(end),3,:));
LW_sub_CLM_eval_AlbSurf06 = squeeze(LW_sub_CLM_eval_AlbSurf06);
LW_sub_CLM_eval_AlbSurf07 = vertcat(LW_in_bc_CLM_AlbSurf(EP_2007(1):gap_2007_1(1)-1,4,:),...
    LW_in_bc_CLM_AlbSurf(gap_2007_1(end)+1:gap_2007_2(1)-1,4,:),...
    LW_in_bc_CLM_AlbSurf(gap_2007_2(end)+1:EP_2007(end),4,:));
LW_sub_CLM_eval_AlbSurf07 = squeeze(LW_sub_CLM_eval_AlbSurf07);
    % SNOWPACK
LW_sub_SP_eval_AlbSurf04 = vertcat(LW_in_bc_SP_AlbSurf(EP_2004(1):gap_2004_1(1)-1,1,:),...
    LW_in_bc_SP_AlbSurf(gap_2004_1(end)+1:gap_2004_2(1)-1,1,:),...
    LW_in_bc_SP_AlbSurf(gap_2004_2(end)+1:gap_2004_3(1)-1,1,:),...
    LW_in_bc_SP_AlbSurf(gap_2004_3(end)+1:gap_2004_4(1)-1,1,:),...
    LW_in_bc_SP_AlbSurf(gap_2004_4(end)+1:EP_2004(end),1,:));
LW_sub_SP_eval_AlbSurf04 = squeeze(LW_sub_SP_eval_AlbSurf04);
LW_sub_SP_eval_AlbSurf05 = vertcat(LW_in_bc_SP_AlbSurf(EP_2005(1):gap_2005_1(1)-1,2,:),...
    LW_in_bc_SP_AlbSurf(gap_2005_1(end)+1:gap_2005_2(1)-1,2,:),...
    LW_in_bc_SP_AlbSurf(gap_2005_2(end)+1:gap_2005_3(1)-1,2,:),...
    LW_in_bc_SP_AlbSurf(gap_2005_3(end)+1:gap_2005_4(1)-1,2,:),...
    LW_in_bc_SP_AlbSurf(gap_2005_4(end)+1:EP_2005(end),2,:));
LW_sub_SP_eval_AlbSurf05 = squeeze(LW_sub_SP_eval_AlbSurf05);
LW_sub_SP_eval_AlbSurf06 = vertcat(LW_in_bc_SP_AlbSurf(EP_2006(1):gap_2006(1)-1,3,:),...
    LW_in_bc_SP_AlbSurf(gap_2006(end)+1:EP_2006(end),3,:));
LW_sub_SP_eval_AlbSurf06 = squeeze(LW_sub_SP_eval_AlbSurf06);
LW_sub_SP_eval_AlbSurf07 = vertcat(LW_in_bc_SP_AlbSurf(EP_2007(1):gap_2007_1(1)-1,4,:),...
    LW_in_bc_SP_AlbSurf(gap_2007_1(end)+1:gap_2007_2(1)-1,4,:),...
    LW_in_bc_SP_AlbSurf(gap_2007_2(end)+1:EP_2007(end),4,:));
LW_sub_SP_eval_AlbSurf07 = squeeze(LW_sub_SP_eval_AlbSurf07);


% diffuse fraction
    % CLM
LW_sub_CLM_eval_Fdiff04 = vertcat(LW_in_bc_CLM_Fdiff(EP_2004(1):gap_2004_1(1)-1,1,:),...
    LW_in_bc_CLM_Fdiff(gap_2004_1(end)+1:gap_2004_2(1)-1,1,:),...
    LW_in_bc_CLM_Fdiff(gap_2004_2(end)+1:gap_2004_3(1)-1,1,:),...
    LW_in_bc_CLM_Fdiff(gap_2004_3(end)+1:gap_2004_4(1)-1,1,:),...
    LW_in_bc_CLM_Fdiff(gap_2004_4(end)+1:EP_2004(end),1,:));
LW_sub_CLM_eval_Fdiff04 = squeeze(LW_sub_CLM_eval_Fdiff04);
LW_sub_CLM_eval_Fdiff05 = vertcat(LW_in_bc_CLM_Fdiff(EP_2005(1):gap_2005_1(1)-1,2,:),...
    LW_in_bc_CLM_Fdiff(gap_2005_1(end)+1:gap_2005_2(1)-1,2,:),...
    LW_in_bc_CLM_Fdiff(gap_2005_2(end)+1:gap_2005_3(1)-1,2,:),...
    LW_in_bc_CLM_Fdiff(gap_2005_3(end)+1:gap_2005_4(1)-1,2,:),...
    LW_in_bc_CLM_Fdiff(gap_2005_4(end)+1:EP_2005(end),2,:));
LW_sub_CLM_eval_Fdiff05 = squeeze(LW_sub_CLM_eval_Fdiff05);
LW_sub_CLM_eval_Fdiff06 = vertcat(LW_in_bc_CLM_Fdiff(EP_2006(1):gap_2006(1)-1,3,:),...
    LW_in_bc_CLM_Fdiff(gap_2006(end)+1:EP_2006(end),3,:));
LW_sub_CLM_eval_Fdiff06 = squeeze(LW_sub_CLM_eval_Fdiff06);
LW_sub_CLM_eval_Fdiff07 = vertcat(LW_in_bc_CLM_Fdiff(EP_2007(1):gap_2007_1(1)-1,4,:),...
    LW_in_bc_CLM_Fdiff(gap_2007_1(end)+1:gap_2007_2(1)-1,4,:),...
    LW_in_bc_CLM_Fdiff(gap_2007_2(end)+1:EP_2007(end),4,:));
LW_sub_CLM_eval_Fdiff07 = squeeze(LW_sub_CLM_eval_Fdiff07);
    % SNOWPACK
LW_sub_SP_eval_Fdiff04 = vertcat(LW_in_bc_SP_Fdiff(EP_2004(1):gap_2004_1(1)-1,1,:),...
    LW_in_bc_SP_Fdiff(gap_2004_1(end)+1:gap_2004_2(1)-1,1,:),...
    LW_in_bc_SP_Fdiff(gap_2004_2(end)+1:gap_2004_3(1)-1,1,:),...
    LW_in_bc_SP_Fdiff(gap_2004_3(end)+1:gap_2004_4(1)-1,1,:),...
    LW_in_bc_SP_Fdiff(gap_2004_4(end)+1:EP_2004(end),1,:));
LW_sub_SP_eval_Fdiff04 = squeeze(LW_sub_SP_eval_Fdiff04);
LW_sub_SP_eval_Fdiff05 = vertcat(LW_in_bc_SP_Fdiff(EP_2005(1):gap_2005_1(1)-1,2,:),...
    LW_in_bc_SP_Fdiff(gap_2005_1(end)+1:gap_2005_2(1)-1,2,:),...
    LW_in_bc_SP_Fdiff(gap_2005_2(end)+1:gap_2005_3(1)-1,2,:),...
    LW_in_bc_SP_Fdiff(gap_2005_3(end)+1:gap_2005_4(1)-1,2,:),...
    LW_in_bc_SP_Fdiff(gap_2005_4(end)+1:EP_2005(end),2,:));
LW_sub_SP_eval_Fdiff05 = squeeze(LW_sub_SP_eval_Fdiff05);
LW_sub_SP_eval_Fdiff06 = vertcat(LW_in_bc_SP_Fdiff(EP_2006(1):gap_2006(1)-1,3,:),...
    LW_in_bc_SP_Fdiff(gap_2006(end)+1:EP_2006(end),3,:));
LW_sub_SP_eval_Fdiff06 = squeeze(LW_sub_SP_eval_Fdiff06);
LW_sub_SP_eval_Fdiff07 = vertcat(LW_in_bc_SP_Fdiff(EP_2007(1):gap_2007_1(1)-1,4,:),...
    LW_in_bc_SP_Fdiff(gap_2007_1(end)+1:gap_2007_2(1)-1,4,:),...
    LW_in_bc_SP_Fdiff(gap_2007_2(end)+1:EP_2007(end),4,:));
LW_sub_SP_eval_Fdiff07 = squeeze(LW_sub_SP_eval_Fdiff07);

% soil water content
    % CLM
LW_sub_CLM_eval_SWC04 = vertcat(LW_in_bc_CLM_SWC(EP_2004(1):gap_2004_1(1)-1,1,:),...
    LW_in_bc_CLM_SWC(gap_2004_1(end)+1:gap_2004_2(1)-1,1,:),...
    LW_in_bc_CLM_SWC(gap_2004_2(end)+1:gap_2004_3(1)-1,1,:),...
    LW_in_bc_CLM_SWC(gap_2004_3(end)+1:gap_2004_4(1)-1,1,:),...
    LW_in_bc_CLM_SWC(gap_2004_4(end)+1:EP_2004(end),1,:));
LW_sub_CLM_eval_SWC04 = squeeze(LW_sub_CLM_eval_SWC04);
LW_sub_CLM_eval_SWC05 = vertcat(LW_in_bc_CLM_SWC(EP_2005(1):gap_2005_1(1)-1,2,:),...
    LW_in_bc_CLM_SWC(gap_2005_1(end)+1:gap_2005_2(1)-1,2,:),...
    LW_in_bc_CLM_SWC(gap_2005_2(end)+1:gap_2005_3(1)-1,2,:),...
    LW_in_bc_CLM_SWC(gap_2005_3(end)+1:gap_2005_4(1)-1,2,:),...
    LW_in_bc_CLM_SWC(gap_2005_4(end)+1:EP_2005(end),2,:));
LW_sub_CLM_eval_SWC05 = squeeze(LW_sub_CLM_eval_SWC05);
LW_sub_CLM_eval_SWC06 = vertcat(LW_in_bc_CLM_SWC(EP_2006(1):gap_2006(1)-1,3,:),...
    LW_in_bc_CLM_SWC(gap_2006(end)+1:EP_2006(end),3,:));
LW_sub_CLM_eval_SWC06 = squeeze(LW_sub_CLM_eval_SWC06);
LW_sub_CLM_eval_SWC07 = vertcat(LW_in_bc_CLM_SWC(EP_2007(1):gap_2007_1(1)-1,4,:),...
    LW_in_bc_CLM_SWC(gap_2007_1(end)+1:gap_2007_2(1)-1,4,:),...
    LW_in_bc_CLM_SWC(gap_2007_2(end)+1:EP_2007(end),4,:));
LW_sub_CLM_eval_SWC07 = squeeze(LW_sub_CLM_eval_SWC07);
    % SNOWPACK
LW_sub_SP_eval_SWC04 = vertcat(LW_in_bc_SP_SWC(EP_2004(1):gap_2004_1(1)-1,1,:),...
    LW_in_bc_SP_SWC(gap_2004_1(end)+1:gap_2004_2(1)-1,1,:),...
    LW_in_bc_SP_SWC(gap_2004_2(end)+1:gap_2004_3(1)-1,1,:),...
    LW_in_bc_SP_SWC(gap_2004_3(end)+1:gap_2004_4(1)-1,1,:),...
    LW_in_bc_SP_SWC(gap_2004_4(end)+1:EP_2004(end),1,:));
LW_sub_SP_eval_SWC04 = squeeze(LW_sub_SP_eval_SWC04);
LW_sub_SP_eval_SWC05 = vertcat(LW_in_bc_SP_SWC(EP_2005(1):gap_2005_1(1)-1,2,:),...
    LW_in_bc_SP_SWC(gap_2005_1(end)+1:gap_2005_2(1)-1,2,:),...
    LW_in_bc_SP_SWC(gap_2005_2(end)+1:gap_2005_3(1)-1,2,:),...
    LW_in_bc_SP_SWC(gap_2005_3(end)+1:gap_2005_4(1)-1,2,:),...
    LW_in_bc_SP_SWC(gap_2005_4(end)+1:EP_2005(end),2,:));
LW_sub_SP_eval_SWC05 = squeeze(LW_sub_SP_eval_SWC05);
LW_sub_SP_eval_SWC06 = vertcat(LW_in_bc_SP_SWC(EP_2006(1):gap_2006(1)-1,3,:),...
    LW_in_bc_SP_SWC(gap_2006(end)+1:EP_2006(end),3,:));
LW_sub_SP_eval_SWC06 = squeeze(LW_sub_SP_eval_SWC06);
LW_sub_SP_eval_SWC07 = vertcat(LW_in_bc_SP_SWC(EP_2007(1):gap_2007_1(1)-1,4,:),...
    LW_in_bc_SP_SWC(gap_2007_1(end)+1:gap_2007_2(1)-1,4,:),...
    LW_in_bc_SP_SWC(gap_2007_2(end)+1:EP_2007(end),4,:));
LW_sub_SP_eval_SWC07 = squeeze(LW_sub_SP_eval_SWC07);

%--------------------------  concatenate years  --------------------------%
LW_sub_val_eval = vertcat(LW_sub_val_eval04,LW_sub_val_eval05,...
    LW_sub_val_eval06,LW_sub_val_eval07);
LW_sub_CLM_eval_AlbSurf = vertcat(LW_sub_CLM_eval_AlbSurf04,LW_sub_CLM_eval_AlbSurf05,...
    LW_sub_CLM_eval_AlbSurf06,LW_sub_CLM_eval_AlbSurf07);
LW_sub_CLM_eval_Fdiff = vertcat(LW_sub_CLM_eval_Fdiff04,LW_sub_CLM_eval_Fdiff05,...
    LW_sub_CLM_eval_Fdiff06,LW_sub_CLM_eval_Fdiff07);
LW_sub_CLM_eval_SWC = vertcat(LW_sub_CLM_eval_SWC04,LW_sub_CLM_eval_SWC05,...
    LW_sub_CLM_eval_SWC06,LW_sub_CLM_eval_SWC07);
LW_sub_SP_eval_AlbSurf = vertcat(LW_sub_SP_eval_AlbSurf04,LW_sub_SP_eval_AlbSurf05,...
    LW_sub_SP_eval_AlbSurf06,LW_sub_SP_eval_AlbSurf07);
LW_sub_SP_eval_Fdiff = vertcat(LW_sub_SP_eval_Fdiff04,LW_sub_SP_eval_Fdiff05,...
    LW_sub_SP_eval_Fdiff06,LW_sub_SP_eval_Fdiff07);
LW_sub_SP_eval_SWC = vertcat(LW_sub_SP_eval_SWC04,LW_sub_SP_eval_SWC05,...
    LW_sub_SP_eval_SWC06,LW_sub_SP_eval_SWC07);


%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%
% CLM
RMSE_CLM = nan(11,3); MBD_CLM = nan(11,3);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_AlbSurf(:,j),LW_sub_val_eval);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_AlbSurf(:,j),LW_sub_val_eval);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_Fdiff(:,j),LW_sub_val_eval);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_Fdiff(:,j),LW_sub_val_eval);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_SWC(:,j),LW_sub_val_eval);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_SWC(:,j),LW_sub_val_eval);
end
RMSE_CLM_change = nan(size(RMSE_CLM));
for i=1:3
    RMSE_CLM_change(:,i) = RMSE_CLM(:,i) - RMSE_CLM(6,i);
end
MBD_CLM_change = nan(size(MBD_CLM));
for i=1:3
    MBD_CLM_change(:,i) = MBD_CLM(:,i) - MBD_CLM(6,i);
end

% SNOWPACK
RMSE_SP = nan(11,3); MBD_SP = nan(11,3);
for j=1:11
    RMSE_SP(j,1) = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval_AlbSurf(:,j),LW_sub_val_eval);
    MBD_SP(j,1) = MBD(length(LW_sub_val_eval),LW_sub_SP_eval_AlbSurf(:,j),LW_sub_val_eval);
    RMSE_SP(j,2) = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval_Fdiff(:,j),LW_sub_val_eval);
    MBD_SP(j,2) = MBD(length(LW_sub_val_eval),LW_sub_SP_eval_Fdiff(:,j),LW_sub_val_eval);
    RMSE_SP(j,3) = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval_SWC(:,j),LW_sub_val_eval);
    MBD_SP(j,3) = MBD(length(LW_sub_val_eval),LW_sub_SP_eval_SWC(:,j),LW_sub_val_eval);
end
RMSE_SP_change = nan(size(RMSE_SP));
for i=1:3
    RMSE_SP_change(:,i) = RMSE_SP(:,i) - RMSE_SP(6,i);
end
MBD_SP_change = nan(size(MBD_SP));
for i=1:3
    MBD_SP_change(:,i) = MBD_SP(:,i) - MBD_SP(6,i);
end

%------------------------  longwave  enhancement  ------------------------%
LW_atm_eval = vertcat(LWR_Atmosphere_Alptal_2004,LWR_Atmosphere_Alptal_2005,...
    LWR_Atmosphere_Alptal_2006,LWR_Atmosphere_Alptal_2007);

LWE_sub_val_eval = LW_sub_val_eval./LW_atm_eval;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_AlbSurf));
LWE_sub_SP_eval_AlbSurf = nan(size(LW_sub_SP_eval_AlbSurf));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) = LW_sub_CLM_eval_AlbSurf(:,j)./LW_atm_eval;
    LWE_sub_SP_eval_AlbSurf(:,j) = LW_sub_SP_eval_AlbSurf(:,j)./LW_atm_eval;
end

% diffuse fraction
LWE_sub_CLM_eval_Fdiff = nan(size(LW_sub_CLM_eval_Fdiff));
LWE_sub_SP_eval_Fdiff = nan(size(LW_sub_SP_eval_Fdiff));
for j=1:11
    LWE_sub_CLM_eval_Fdiff(:,j) = LW_sub_CLM_eval_Fdiff(:,j)./LW_atm_eval;
    LWE_sub_SP_eval_Fdiff(:,j) = LW_sub_SP_eval_Fdiff(:,j)./LW_atm_eval;
end

% soil water content
LWE_sub_CLM_eval_SWC = nan(size(LW_sub_CLM_eval_SWC));
LWE_sub_SP_eval_SWC = nan(size(LW_sub_SP_eval_SWC));
for j=1:11
    LWE_sub_CLM_eval_SWC(:,j) = LW_sub_CLM_eval_SWC(:,j)./LW_atm_eval;
    LWE_sub_SP_eval_SWC(:,j) = LW_sub_SP_eval_SWC(:,j)./LW_atm_eval;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_Alp = hist_obs.Values;
hist_clm_Alp_AlbSurf = nan(length(hist_obs_Alp),11);
hist_sp_Alp_AlbSurf = nan(length(hist_obs_Alp),11);
hist_clm_Alp_Fdiff = nan(length(hist_obs_Alp),11);
hist_sp_Alp_Fdiff = nan(length(hist_obs_Alp),11);
hist_clm_Alp_SWC = nan(length(hist_obs_Alp),11);
hist_sp_Alp_SWC = nan(length(hist_obs_Alp),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Alp_AlbSurf(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_AlbSurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_sp_Alp_AlbSurf(:,j) = hist_sp.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Fdiff(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Alp_Fdiff(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_Fdiff(:,j),spectrum_LWenh,'Normalization','probability');
    hist_sp_Alp_Fdiff(:,j) = hist_sp.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_SWC(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Alp_SWC(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_SWC(:,j),spectrum_LWenh,'Normalization','probability');
    hist_sp_Alp_SWC(:,j) = hist_sp.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
hold on
text(0.8+0.05*1.2,0.95*0.3,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Alp,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Alp_AlbSurf,'Color',EcstaticEmerald)
plot(spectrum_LWenh_xaxis,hist_sp_Alp_AlbSurf,'Color',CandidCoral)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(1,3,2)
hold on
text(0.8+0.05*1.2,0.95*0.3,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Alp,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Alp_Fdiff,'Color',EcstaticEmerald)
plot(spectrum_LWenh_xaxis,hist_sp_Alp_Fdiff,'Color',CandidCoral)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(1,3,3)
hold on
text(0.8+0.05*1.2,0.95*0.3,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Alp,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Alp_SWC,'Color',EcstaticEmerald)
plot(spectrum_LWenh_xaxis,hist_sp_Alp_SWC,'Color',CandidCoral)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Alptal.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 0 33.5 10];
fig.PaperSize = [28.5 9.5];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Alptal.pdf')