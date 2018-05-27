clear
close all

load LWsub_Seehornwald_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_AlbSurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Seehornwald_Fdiff.mat
LW_in_bc_CLM_Fdiff = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_Fdiff = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Seehornwald_Zsnow.mat
LW_in_bc_CLM_Zsnow = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_Zsnow = LW_in_bc_SP;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Seehornwald.mat','LWR_Atmosphere_Seehornwald_2008',...
    'LWR_Atmosphere_Seehornwald_2009','LWR_Atmosphere_Seehornwald_2010',...
    'LWR_Atmosphere_Seehornwald_2011','LWR_Atmosphere_Seehornwald_2012')
LWR_Atmosphere_Seehornwald = vertcat(LWR_Atmosphere_Seehornwald_2008,...
    LWR_Atmosphere_Seehornwald_2009,LWR_Atmosphere_Seehornwald_2010,...
    LWR_Atmosphere_Seehornwald_2011,LWR_Atmosphere_Seehornwald_2012);

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
CandidCoral = [1 0.44 0.32];
MagicMaroon = [0.65 0.32 0.35];

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
EP_2008 = 1921:4728;
EP_2009 = 10705:12864;
EP_2010 = 19465:22080;
EP_2011 = 28225:30312;
EP_2012 = 36985:39768;

%--------------  find & skip gaps and uncertain time steps  --------------%
% due to met forcing
gap_2010 = 1200:1224;   % 20 Feb 0:00 - 21 Feb 0:00
gap_2011 = 408:456;     % 18 Jan 0:00 - 20 Jan 0:00
% due to evaluation data
gap_LWsub_2008 = 159;   % 7 Jan 15:00
gap_LWsub_2010 = 134:137;   % 6 Jan 14:00 - 17:00
gap_LWsub_2010_2 = 1200;    % 20 Feb 0:00
gap_LWsub_2011 = 346:349;       % 15 Jan 10:00 - 13:00
gap_LWsub_2011_2 = 427:442;     % 18 Jan 19:00 - 19 Jan 10:00


%------------------------------  year 2008  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_1h(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_val_eval_all = LW_sub_val_eval;
% surface albedo
LW_sub_CLM_eval_AlbSurf = vertcat(LW_in_bc_CLM_AlbSurf(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_CLM_AlbSurf(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_SP_eval_AlbSurf = vertcat(LW_in_bc_SP_AlbSurf(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_SP_AlbSurf(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_CLM_eval_AlbSurf_all = LW_sub_CLM_eval_AlbSurf;
LW_sub_SP_eval_AlbSurf_all = LW_sub_SP_eval_AlbSurf;
% diffuse fraction
LW_sub_CLM_eval_Fdiff = vertcat(LW_in_bc_CLM_Fdiff(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_CLM_Fdiff(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_SP_eval_Fdiff = vertcat(LW_in_bc_SP_Fdiff(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_SP_Fdiff(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_CLM_eval_Fdiff_all = LW_sub_CLM_eval_Fdiff;
LW_sub_SP_eval_Fdiff_all = LW_sub_SP_eval_Fdiff;
% snow depth
LW_sub_CLM_eval_Zsnow = vertcat(LW_in_bc_CLM_Zsnow(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_CLM_Zsnow(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_SP_eval_Zsnow = vertcat(LW_in_bc_SP_Zsnow(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15,:),...
    LW_in_bc_SP_Zsnow(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end),:));
LW_sub_CLM_eval_Zsnow_all = LW_sub_CLM_eval_Zsnow;
LW_sub_SP_eval_Zsnow_all = LW_sub_SP_eval_Zsnow;


%------------------------------  year 2009  ------------------------------%
LW_sub_val_eval = LW_in_bc_1h(EP_2009(1):EP_2009(end),:);
temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
% surface albedo
LW_sub_CLM_eval_AlbSurf = LW_in_bc_CLM_AlbSurf(EP_2009(1):EP_2009(end),:);
LW_sub_SP_eval_AlbSurf = LW_in_bc_SP_AlbSurf(EP_2009(1):EP_2009(end),:);
temp = LW_sub_CLM_eval_AlbSurf_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_AlbSurf_all = vertcat(temp,LW_sub_CLM_eval_AlbSurf);
temp = LW_sub_SP_eval_AlbSurf_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_AlbSurf_all = vertcat(temp,LW_sub_SP_eval_AlbSurf);
% diffuse fraction
LW_sub_CLM_eval_Fdiff = LW_in_bc_CLM_Fdiff(EP_2009(1):EP_2009(end),:);
LW_sub_SP_eval_Fdiff = LW_in_bc_SP_Fdiff(EP_2009(1):EP_2009(end),:);
temp = LW_sub_CLM_eval_Fdiff_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Fdiff_all = vertcat(temp,LW_sub_CLM_eval_Fdiff);
temp = LW_sub_SP_eval_Fdiff_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Fdiff_all = vertcat(temp,LW_sub_SP_eval_Fdiff);
% snow depth
LW_sub_CLM_eval_Zsnow = LW_in_bc_CLM_Zsnow(EP_2009(1):EP_2009(end),:);
LW_sub_SP_eval_Zsnow = LW_in_bc_SP_Zsnow(EP_2009(1):EP_2009(end),:);
temp = LW_sub_CLM_eval_Zsnow_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Zsnow_all = vertcat(temp,LW_sub_CLM_eval_Zsnow);
temp = LW_sub_SP_eval_Zsnow_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Zsnow_all = vertcat(temp,LW_sub_SP_eval_Zsnow);
    
%------------------------------  year 2010  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_1h(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_1h(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
% surface albedo
LW_sub_CLM_eval_AlbSurf = vertcat(LW_in_bc_CLM_AlbSurf(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_CLM_AlbSurf(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_CLM_AlbSurf(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
LW_sub_SP_eval_AlbSurf = vertcat(LW_in_bc_SP_AlbSurf(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_SP_AlbSurf(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_SP_AlbSurf(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
temp = LW_sub_CLM_eval_AlbSurf_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_AlbSurf_all = vertcat(temp,LW_sub_CLM_eval_AlbSurf);
temp = LW_sub_SP_eval_AlbSurf_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_AlbSurf_all = vertcat(temp,LW_sub_SP_eval_AlbSurf);
% diffuse fraction
LW_sub_CLM_eval_Fdiff = vertcat(LW_in_bc_CLM_Fdiff(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_CLM_Fdiff(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_CLM_Fdiff(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
LW_sub_SP_eval_Fdiff = vertcat(LW_in_bc_SP_Fdiff(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_SP_Fdiff(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_SP_Fdiff(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
temp = LW_sub_CLM_eval_Fdiff_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Fdiff_all = vertcat(temp,LW_sub_CLM_eval_Fdiff);
temp = LW_sub_SP_eval_Fdiff_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Fdiff_all = vertcat(temp,LW_sub_SP_eval_Fdiff);
% snow depth
LW_sub_CLM_eval_Zsnow = vertcat(LW_in_bc_CLM_Zsnow(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_CLM_Zsnow(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_CLM_Zsnow(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
LW_sub_SP_eval_Zsnow = vertcat(LW_in_bc_SP_Zsnow(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14,:),...
    LW_in_bc_SP_Zsnow(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24,:),...
    LW_in_bc_SP_Zsnow(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end),:));
temp = LW_sub_CLM_eval_Zsnow_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Zsnow_all = vertcat(temp,LW_sub_CLM_eval_Zsnow);
temp = LW_sub_SP_eval_Zsnow_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Zsnow_all = vertcat(temp,LW_sub_SP_eval_Zsnow);
    
    
%------------------------------  year 2011  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_1h(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_1h(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
% surface albedo
LW_sub_CLM_eval_AlbSurf = vertcat(LW_in_bc_CLM_AlbSurf(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_CLM_AlbSurf(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_CLM_AlbSurf(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
LW_sub_SP_eval_AlbSurf = vertcat(LW_in_bc_SP_AlbSurf(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_SP_AlbSurf(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_SP_AlbSurf(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
temp = LW_sub_CLM_eval_AlbSurf_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_AlbSurf_all = vertcat(temp,LW_sub_CLM_eval_AlbSurf);
temp = LW_sub_SP_eval_AlbSurf_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_AlbSurf_all = vertcat(temp,LW_sub_SP_eval_AlbSurf);
% diffuse fraction
LW_sub_CLM_eval_Fdiff = vertcat(LW_in_bc_CLM_Fdiff(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_CLM_Fdiff(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_CLM_Fdiff(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
LW_sub_SP_eval_Fdiff = vertcat(LW_in_bc_SP_Fdiff(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_SP_Fdiff(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_SP_Fdiff(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
temp = LW_sub_CLM_eval_Fdiff_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Fdiff_all = vertcat(temp,LW_sub_CLM_eval_Fdiff);
temp = LW_sub_SP_eval_Fdiff_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Fdiff_all = vertcat(temp,LW_sub_SP_eval_Fdiff);
% snow depth
LW_sub_CLM_eval_Zsnow = vertcat(LW_in_bc_CLM_Zsnow(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_CLM_Zsnow(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_CLM_Zsnow(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
LW_sub_SP_eval_Zsnow = vertcat(LW_in_bc_SP_Zsnow(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10,:),...
    LW_in_bc_SP_Zsnow(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24,:),...
    LW_in_bc_SP_Zsnow(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end),:));
temp = LW_sub_CLM_eval_Zsnow_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Zsnow_all = vertcat(temp,LW_sub_CLM_eval_Zsnow);
temp = LW_sub_SP_eval_Zsnow_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Zsnow_all = vertcat(temp,LW_sub_SP_eval_Zsnow);

%------------------------------  year 2012  ------------------------------%
LW_sub_val_eval = LW_in_bc_1h(EP_2012(1):EP_2012(end),:);
temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
% surface albedo
LW_sub_CLM_eval_AlbSurf = LW_in_bc_CLM_AlbSurf(EP_2012(1):EP_2012(end),:);
LW_sub_SP_eval_AlbSurf = LW_in_bc_SP_AlbSurf(EP_2012(1):EP_2012(end),:);
temp = LW_sub_CLM_eval_AlbSurf_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_AlbSurf_all = vertcat(temp,LW_sub_CLM_eval_AlbSurf);
temp = LW_sub_SP_eval_AlbSurf_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_AlbSurf_all = vertcat(temp,LW_sub_SP_eval_AlbSurf);
% diffuse fraction
LW_sub_CLM_eval_Fdiff = LW_in_bc_CLM_Fdiff(EP_2012(1):EP_2012(end),:);
LW_sub_SP_eval_Fdiff = LW_in_bc_SP_Fdiff(EP_2012(1):EP_2012(end),:);
temp = LW_sub_CLM_eval_Fdiff_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Fdiff_all = vertcat(temp,LW_sub_CLM_eval_Fdiff);
temp = LW_sub_SP_eval_Fdiff_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Fdiff_all = vertcat(temp,LW_sub_SP_eval_Fdiff);
% snow depth
LW_sub_CLM_eval_Zsnow = LW_in_bc_CLM_Zsnow(EP_2012(1):EP_2012(end),:);
LW_sub_SP_eval_Zsnow = LW_in_bc_SP_Zsnow(EP_2012(1):EP_2012(end),:);
temp = LW_sub_CLM_eval_Zsnow_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_Zsnow_all = vertcat(temp,LW_sub_CLM_eval_Zsnow);
temp = LW_sub_SP_eval_Zsnow_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_Zsnow_all = vertcat(temp,LW_sub_SP_eval_Zsnow);
    
    
%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%

% CLM
RMSE_CLM = nan(11,3); MBD_CLM = nan(11,3);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLM_eval_AlbSurf_all(:,j),LW_sub_val_eval_all);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval_all),LW_sub_CLM_eval_AlbSurf_all(:,j),LW_sub_val_eval_all);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLM_eval_Fdiff_all(:,j),LW_sub_val_eval_all);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval_all),LW_sub_CLM_eval_Fdiff_all(:,j),LW_sub_val_eval_all);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLM_eval_Zsnow_all(:,j),LW_sub_val_eval_all);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval_all),LW_sub_CLM_eval_Zsnow_all(:,j),LW_sub_val_eval_all);
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
    RMSE_SP(j,1) = RMSE(length(LW_sub_val_eval_all),LW_sub_SP_eval_AlbSurf_all(:,j),LW_sub_val_eval_all);
    MBD_SP(j,1) = MBD(length(LW_sub_val_eval_all),LW_sub_SP_eval_AlbSurf_all(:,j),LW_sub_val_eval_all);
    RMSE_SP(j,2) = RMSE(length(LW_sub_val_eval_all),LW_sub_SP_eval_Fdiff_all(:,j),LW_sub_val_eval_all);
    MBD_SP(j,2) = MBD(length(LW_sub_val_eval_all),LW_sub_SP_eval_Fdiff_all(:,j),LW_sub_val_eval_all);
    RMSE_SP(j,3) = RMSE(length(LW_sub_val_eval_all),LW_sub_SP_eval_Zsnow_all(:,j),LW_sub_val_eval_all);
    MBD_SP(j,3) = MBD(length(LW_sub_val_eval_all),LW_sub_SP_eval_Zsnow_all(:,j),LW_sub_val_eval_all);
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
LWE_sub_val_eval = LW_sub_val_eval_all./LWR_Atmosphere_Seehornwald;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_AlbSurf_all));
LWE_sub_SP_eval_AlbSurf = nan(size(LW_sub_SP_eval_AlbSurf_all));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) = LW_sub_CLM_eval_AlbSurf_all(:,j)./LWR_Atmosphere_Seehornwald;
    LWE_sub_SP_eval_AlbSurf(:,j) = LW_sub_SP_eval_AlbSurf_all(:,j)./LWR_Atmosphere_Seehornwald;
end

% diffuse fraction
LWE_sub_CLM_eval_Fdiff = nan(size(LW_sub_CLM_eval_Fdiff_all));
LWE_sub_SP_eval_Fdiff = nan(size(LW_sub_SP_eval_Fdiff_all));
for j=1:11
    LWE_sub_CLM_eval_Fdiff(:,j) = LW_sub_CLM_eval_Fdiff_all(:,j)./LWR_Atmosphere_Seehornwald;
    LWE_sub_SP_eval_Fdiff(:,j) = LW_sub_SP_eval_Fdiff_all(:,j)./LWR_Atmosphere_Seehornwald;
end

% snow depth
LWE_sub_CLM_eval_Zsnow = nan(size(LW_sub_CLM_eval_Zsnow_all));
LWE_sub_SP_eval_Zsnow = nan(size(LW_sub_SP_eval_Zsnow_all));
for j=1:11
    LWE_sub_CLM_eval_Zsnow(:,j) = LW_sub_CLM_eval_Zsnow_all(:,j)./LWR_Atmosphere_Seehornwald;
    LWE_sub_SP_eval_Zsnow(:,j) = LW_sub_SP_eval_Zsnow_all(:,j)./LWR_Atmosphere_Seehornwald;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_SHW = hist_obs.Values;
hist_clm_SHW_AlbSurf = nan(length(hist_obs_SHW),11);
hist_sp_SHW_AlbSurf = nan(length(hist_obs_SHW),11);
hist_clm_SHW_Fdiff = nan(length(hist_obs_SHW),11);
hist_sp_SHW_Fdiff = nan(length(hist_obs_SHW),11);
hist_clm_SHW_Zsnow = nan(length(hist_obs_SHW),11);
hist_sp_SHW_Zsnow = nan(length(hist_obs_SHW),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_SHW_AlbSurf(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_AlbSurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_sp_SHW_AlbSurf(:,j) = hist_sp.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Fdiff(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_SHW_Fdiff(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_Fdiff(:,j),spectrum_LWenh,'Normalization','probability');
    hist_sp_SHW_Fdiff(:,j) = hist_sp.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Zsnow(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_SHW_Zsnow(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_Zsnow(:,j),spectrum_LWenh,'Normalization','probability');
    hist_sp_SHW_Zsnow(:,j) = hist_sp.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
hold on
text(0.8+0.05*1.2,0.95*0.15,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_SHW,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_SHW_AlbSurf,'Color',MagicMaroon)
plot(spectrum_LWenh_xaxis,hist_sp_SHW_AlbSurf,'Color',CandidCoral)
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
text(0.8+0.05*1.2,0.95*0.15,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_SHW,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_SHW_Fdiff,'Color',MagicMaroon)
plot(spectrum_LWenh_xaxis,hist_sp_SHW_Fdiff,'Color',CandidCoral)
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
text(0.8+0.05*1.2,0.95*0.15,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_SHW,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_SHW_Zsnow,'Color',MagicMaroon)
plot(spectrum_LWenh_xaxis,hist_sp_SHW_Zsnow,'Color',CandidCoral)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Seehornwald.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 0 33.5 10];
fig.PaperSize = [28.5 9.5];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Seehornwald.pdf')