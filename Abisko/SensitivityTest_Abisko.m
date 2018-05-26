clear
close all

% deciduous -> no second layer -> no SNOWPACK-2LHM simulation
load LWsub_Abisko_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_AlbSurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Abisko_SWC.mat
LW_in_bc_CLM_SWC = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_SWC = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Abisko_Tsurf.mat
LW_in_bc_CLM_Tsurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_Tsurf = LW_in_bc_SP;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Abisko_HM.mat','LWR_Atmosphere_Abisko')
LWR_Atmosphere_Abisko = LWR_Atmosphere_Abisko';


%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
GorgeousGold = [1 0.85 0];


%-------------------------------------------------------------------------%
%----------------  select evaluation periods based on QC  ----------------%
EP1 = 11:34;     % evaluation period 1 - 1 day
EP2 = 107:130;   % evaluation period 2 - 1 day
EP3 = 179:226;   % evaluation period 3 - 2 days
EP4 = 395:442;   % evaluation period 4 - 2 days
EP5 = 515:586;   % evaluation period 5 - 3 days

LW_sub_val_eval = vertcat(LW_in_bc_C_1h(EP1(1):EP1(end),:),...
    LW_in_bc_C_1h(EP2(1):EP2(end),:),...
    LW_in_bc_C_1h(EP3(1):EP3(end),:),...
    LW_in_bc_C_1h(EP4(1):EP4(end),:),...
    LW_in_bc_C_1h(EP5(1):EP5(end),:));
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
LW_sub_CLM_eval_avg_AlbSurf = vertcat(LW_sub_CLM_avg_AlbSurf(EP1(1):EP1(end),:),...
    LW_sub_CLM_avg_AlbSurf(EP2(1):EP2(end),:),...
    LW_sub_CLM_avg_AlbSurf(EP3(1):EP3(end),:),...
    LW_sub_CLM_avg_AlbSurf(EP4(1):EP4(end),:),...
    LW_sub_CLM_avg_AlbSurf(EP5(1):EP5(end),:));

% soil water content
LW_sub_CLM_avg_SWC = nan(length(LW_in_bc_CLM_SWC(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_CLM_avg_SWC)
        LW_sub_CLM_avg_SWC(l,j) = mean(LW_in_bc_CLM_SWC(l,:,j));
    end
end
LW_sub_CLM_eval_avg_SWC = vertcat(LW_sub_CLM_avg_SWC(EP1(1):EP1(end),:),...
    LW_sub_CLM_avg_SWC(EP2(1):EP2(end),:),...
    LW_sub_CLM_avg_SWC(EP3(1):EP3(end),:),...
    LW_sub_CLM_avg_SWC(EP4(1):EP4(end),:),...
    LW_sub_CLM_avg_SWC(EP5(1):EP5(end),:));

% surface temperature
LW_sub_CLM_avg_Tsurf = nan(length(LW_in_bc_CLM_Tsurf(:,1,1)),11);
for j=1:11
    for l=1:length(LW_sub_CLM_avg_Tsurf)
        LW_sub_CLM_avg_Tsurf(l,j) = mean(LW_in_bc_CLM_Tsurf(l,:,j));
    end
end
LW_sub_CLM_eval_avg_Tsurf = vertcat(LW_sub_CLM_avg_Tsurf(EP1(1):EP1(end),:),...
    LW_sub_CLM_avg_Tsurf(EP2(1):EP2(end),:),...
    LW_sub_CLM_avg_Tsurf(EP3(1):EP3(end),:),...
    LW_sub_CLM_avg_Tsurf(EP4(1):EP4(end),:),...
    LW_sub_CLM_avg_Tsurf(EP5(1):EP5(end),:));


%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%
% CLM
RMSE_CLM = nan(11,3); MBD_CLM = nan(11,3);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg_AlbSurf(:,j),LW_sub_val_eval_avg);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg_AlbSurf(:,j),LW_sub_val_eval_avg);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg_SWC(:,j),LW_sub_val_eval_avg);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg_SWC(:,j),LW_sub_val_eval_avg);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg_Tsurf(:,j),LW_sub_val_eval_avg);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg_Tsurf(:,j),LW_sub_val_eval_avg);
end
RMSE_CLM_change = nan(size(RMSE_CLM));
for i=1:3
    RMSE_CLM_change(:,i) = RMSE_CLM(:,i) - RMSE_CLM(6,i);
end
MBD_CLM_change = nan(size(MBD_CLM));
for i=1:3
    MBD_CLM_change(:,i) = MBD_CLM(:,i) - MBD_CLM(6,i);
end

%------------------------  longwave  enhancement  ------------------------%
LWE_sub_val_eval = LW_sub_val_eval_avg./LWR_Atmosphere_Abisko;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_avg_AlbSurf));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) = LW_sub_CLM_eval_avg_AlbSurf(:,j)./LWR_Atmosphere_Abisko;
end

% soil water content
LWE_sub_CLM_eval_SWC = nan(size(LW_sub_CLM_eval_avg_SWC));
for j=1:11
    LWE_sub_CLM_eval_SWC(:,j) = LW_sub_CLM_eval_avg_SWC(:,j)./LWR_Atmosphere_Abisko;
end

% surface temperature
LWE_sub_CLM_eval_Tsurf = nan(size(LW_sub_CLM_eval_avg_Tsurf));
for j=1:11
    LWE_sub_CLM_eval_Tsurf(:,j) = LW_sub_CLM_eval_avg_Tsurf(:,j)./LWR_Atmosphere_Abisko;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_Abi = hist_obs.Values;
hist_clm_Abi_AlbSurf = nan(length(hist_obs_Abi),11);
hist_clm_Abi_SWC = nan(length(hist_obs_Abi),11);
hist_clm_Abi_Tsurf = nan(length(hist_obs_Abi),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Abi_AlbSurf(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_SWC(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Abi_SWC(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Tsurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Abi_Tsurf(:,j) = hist_clm.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
hold on
text(0.8+0.05*1.2,0.95*0.4,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Abi,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Abi_AlbSurf,'Color',GorgeousGold)
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
plot(spectrum_LWenh_xaxis,hist_obs_Abi,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Abi_SWC,'Color',GorgeousGold)
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
plot(spectrum_LWenh_xaxis,hist_obs_Abi,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Abi_Tsurf,'Color',GorgeousGold)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Abisko.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 0 33.5 10];
fig.PaperSize = [28.5 9.5];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Abisko.pdf')