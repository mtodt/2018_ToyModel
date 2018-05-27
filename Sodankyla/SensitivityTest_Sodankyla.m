clear
close all

load LWsub_Sodankyla_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_AlbSurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Sodankyla_Tsurf.mat
LW_in_bc_CLM_Tsurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_Tsurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Sodankyla_Zsnow.mat
LW_in_bc_CLM_Zsnow = LW_in_bc_CLM;
clear LW_in_bc_CLM
LW_in_bc_SP_Zsnow = LW_in_bc_SP;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Sodankyla_HM.mat','LWR_Atmosphere_Sodankyla')
LWR_Atmosphere_Sodankyla = LWR_Atmosphere_Sodankyla';

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
GlacialGrey = [0.4 0.7 1];
CandidCoral = [1 0.44 0.32];

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
% C: select 10 March 1:00 to 28 March 0:00
EP_C_1 = 12:443;
% C: add 29 March 1:00 to 17 April 0:00
EP_C_2 = 468:923;

LW_sub_val_eval = vertcat(LW_in_bc_C_1h(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_C_1h(EP_C_2(1):EP_C_2(end),:));
LW_sub_val_eval_avg = nan(length(LW_sub_val_eval(:,1)),1);
for l=1:length(LW_sub_val_eval_avg)
    LW_sub_val_eval_avg(l) = mean(LW_sub_val_eval(l,:));
end

% surface albedo
LW_sub_CLM_avg_AlbSurf = nan(length(LW_in_bc_CLM_AlbSurf(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_CLM_avg_AlbSurf)
        LW_sub_CLM_avg_AlbSurf(l,j) = mean(LW_in_bc_CLM_AlbSurf(l,:,j));
    end
end
LW_sub_CLM_eval_avg_AlbSurf = vertcat(LW_sub_CLM_avg_AlbSurf(EP_C_1(1):EP_C_1(end),:),...
    LW_sub_CLM_avg_AlbSurf(EP_C_2(1):EP_C_2(end),:));
LW_sub_SP_avg_AlbSurf = nan(length(LW_in_bc_SP_AlbSurf(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_SP_avg_AlbSurf)
        LW_sub_SP_avg_AlbSurf(l,j) = mean(LW_in_bc_SP_AlbSurf(l,:,j));
    end
end
LW_sub_SP_eval_avg_AlbSurf = vertcat(LW_sub_SP_avg_AlbSurf(EP_C_1(1):EP_C_1(end),:),...
    LW_sub_SP_avg_AlbSurf(EP_C_2(1):EP_C_2(end),:));

% surface temperature
LW_sub_CLM_avg_Tsurf = nan(length(LW_in_bc_CLM_Tsurf(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_CLM_avg_Tsurf)
        LW_sub_CLM_avg_Tsurf(l,j) = mean(LW_in_bc_CLM_Tsurf(l,:,j));
    end
end
LW_sub_CLM_eval_avg_Tsurf = vertcat(LW_sub_CLM_avg_Tsurf(EP_C_1(1):EP_C_1(end),:),...
    LW_sub_CLM_avg_Tsurf(EP_C_2(1):EP_C_2(end),:));
LW_sub_SP_avg_Tsurf = nan(length(LW_in_bc_SP_Tsurf(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_SP_avg_Tsurf)
        LW_sub_SP_avg_Tsurf(l,j) = mean(LW_in_bc_SP_Tsurf(l,:,j));
    end
end
LW_sub_SP_eval_avg_Tsurf = vertcat(LW_sub_SP_avg_Tsurf(EP_C_1(1):EP_C_1(end),:),...
    LW_sub_SP_avg_Tsurf(EP_C_2(1):EP_C_2(end),:));

% snow depth
LW_sub_CLM_avg_Zsnow = nan(length(LW_in_bc_CLM_Zsnow(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_CLM_avg_Zsnow)
        LW_sub_CLM_avg_Zsnow(l,j) = mean(LW_in_bc_CLM_Zsnow(l,:,j));
    end
end
LW_sub_CLM_eval_avg_Zsnow = vertcat(LW_sub_CLM_avg_Zsnow(EP_C_1(1):EP_C_1(end),:),...
    LW_sub_CLM_avg_Zsnow(EP_C_2(1):EP_C_2(end),:));
LW_sub_SP_avg_Zsnow = nan(length(LW_in_bc_SP_Zsnow(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_SP_avg_Zsnow)
        LW_sub_SP_avg_Zsnow(l,j) = mean(LW_in_bc_SP_Zsnow(l,:,j));
    end
end
LW_sub_SP_eval_avg_Zsnow = vertcat(LW_sub_SP_avg_Zsnow(EP_C_1(1):EP_C_1(end),:),...
    LW_sub_SP_avg_Zsnow(EP_C_2(1):EP_C_2(end),:));


%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%
% CLM
RMSE_CLM = nan(11,3); MBD_CLM = nan(11,3);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval_avg),...
        LW_sub_CLM_eval_avg_AlbSurf(:,j),LW_sub_val_eval_avg);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval_avg),...
        LW_sub_CLM_eval_avg_AlbSurf(:,j),LW_sub_val_eval_avg);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval_avg),...
        LW_sub_CLM_eval_avg_Tsurf(:,j),LW_sub_val_eval_avg);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval_avg),...
        LW_sub_CLM_eval_avg_Tsurf(:,j),LW_sub_val_eval_avg);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval_avg),...
        LW_sub_CLM_eval_avg_Zsnow(:,j),LW_sub_val_eval_avg);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval_avg),...
        LW_sub_CLM_eval_avg_Zsnow(:,j),LW_sub_val_eval_avg);
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
    RMSE_SP(j,1) = RMSE(length(LW_sub_val_eval_avg),...
        LW_sub_SP_eval_avg_AlbSurf(:,j),LW_sub_val_eval_avg);
    MBD_SP(j,1) = MBD(length(LW_sub_val_eval_avg),...
        LW_sub_SP_eval_avg_AlbSurf(:,j),LW_sub_val_eval_avg);
    RMSE_SP(j,2) = RMSE(length(LW_sub_val_eval_avg),...
        LW_sub_SP_eval_avg_Tsurf(:,j),LW_sub_val_eval_avg);
    MBD_SP(j,2) = MBD(length(LW_sub_val_eval_avg),...
        LW_sub_SP_eval_avg_Tsurf(:,j),LW_sub_val_eval_avg);
    RMSE_SP(j,3) = RMSE(length(LW_sub_val_eval_avg),...
        LW_sub_SP_eval_avg_Zsnow(:,j),LW_sub_val_eval_avg);
    MBD_SP(j,3) = MBD(length(LW_sub_val_eval_avg),...
        LW_sub_SP_eval_avg_Zsnow(:,j),LW_sub_val_eval_avg);
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
LWE_sub_val_eval = LW_sub_val_eval_avg./LWR_Atmosphere_Sodankyla;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_avg_AlbSurf));
LWE_sub_SP_eval_AlbSurf = nan(size(LW_sub_SP_eval_avg_AlbSurf));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) ...
        = LW_sub_CLM_eval_avg_AlbSurf(:,j)./LWR_Atmosphere_Sodankyla;
    LWE_sub_SP_eval_AlbSurf(:,j) ...
        = LW_sub_SP_eval_avg_AlbSurf(:,j)./LWR_Atmosphere_Sodankyla;
end

% surface temperature
LWE_sub_CLM_eval_Tsurf = nan(size(LW_sub_CLM_eval_avg_Tsurf));
LWE_sub_SP_eval_Tsurf = nan(size(LW_sub_SP_eval_avg_Tsurf));
for j=1:11
    LWE_sub_CLM_eval_Tsurf(:,j) ...
        = LW_sub_CLM_eval_avg_Tsurf(:,j)./LWR_Atmosphere_Sodankyla;
    LWE_sub_SP_eval_Tsurf(:,j) ...
        = LW_sub_SP_eval_avg_Tsurf(:,j)./LWR_Atmosphere_Sodankyla;
end

% snow depth
LWE_sub_CLM_eval_Zsnow = nan(size(LW_sub_CLM_eval_avg_Zsnow));
LWE_sub_SP_eval_Zsnow = nan(size(LW_sub_SP_eval_avg_Zsnow));
for j=1:11
    LWE_sub_CLM_eval_Zsnow(:,j) ...
        = LW_sub_CLM_eval_avg_Zsnow(:,j)./LWR_Atmosphere_Sodankyla;
    LWE_sub_SP_eval_Zsnow(:,j) ...
        = LW_sub_SP_eval_avg_Zsnow(:,j)./LWR_Atmosphere_Sodankyla;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_Sod = hist_obs.Values;
hist_clm_Sod_AlbSurf = nan(length(hist_obs_Sod),11);
hist_sp_Sod_AlbSurf = nan(length(hist_obs_Sod),11);
hist_clm_Sod_Tsurf = nan(length(hist_obs_Sod),11);
hist_sp_Sod_Tsurf = nan(length(hist_obs_Sod),11);
hist_clm_Sod_Zsnow = nan(length(hist_obs_Sod),11);
hist_sp_Sod_Zsnow = nan(length(hist_obs_Sod),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Sod_AlbSurf(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_AlbSurf(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_sp_Sod_AlbSurf(:,j) = hist_sp.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Tsurf(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Sod_Tsurf(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_Tsurf(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_sp_Sod_Tsurf(:,j) = hist_sp.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Zsnow(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Sod_Zsnow(:,j) = hist_clm.Values;
    hist_sp = histogram(LWE_sub_SP_eval_Zsnow(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_sp_Sod_Zsnow(:,j) = hist_sp.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
hold on
text(0.8+0.05*1.2,0.95*0.4,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Sod,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Sod_AlbSurf,'Color',GlacialGrey)
plot(spectrum_LWenh_xaxis,hist_sp_Sod_AlbSurf,'Color',CandidCoral)
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
text(0.8+0.05*1.2,0.95*0.4,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Sod,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Sod_Tsurf,'Color',GlacialGrey)
plot(spectrum_LWenh_xaxis,hist_sp_Sod_Tsurf,'Color',CandidCoral)
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
text(0.8+0.05*1.2,0.95*0.4,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Sod,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Sod_Zsnow,'Color',GlacialGrey)
plot(spectrum_LWenh_xaxis,hist_sp_Sod_Zsnow,'Color',CandidCoral)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Sodankyla.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 0 33.5 10];
fig.PaperSize = [28.5 9.5];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Sodankyla.pdf')