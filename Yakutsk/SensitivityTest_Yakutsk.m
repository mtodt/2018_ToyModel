clear
close all

% deciduous -> no second layer -> no SNOWPACK-2LHM simulation
load LWsub_Yakutsk_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM_98;
clear LW_in_bc_CLM
% LW_in_bc_SP_AlbSurf = LW_in_bc_SP_98;
clear LW_in_bc_SP
load LWsub_Yakutsk_Fdiff.mat
LW_in_bc_CLM_Fdiff = LW_in_bc_CLM_98;
clear LW_in_bc_CLM
% LW_in_bc_SP_Fdiff = LW_in_bc_SP_98;
clear LW_in_bc_SP
load LWsub_Yakutsk_Zsnow.mat
LW_in_bc_CLM_Zsnow = LW_in_bc_CLM_98;
clear LW_in_bc_CLM
% LW_in_bc_SP_Zsnow = LW_in_bc_SP_98;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Yakutsk.mat','LWR_Atmosphere_Yakutsk_Night')
LWR_Atmosphere_Yakutsk_Night = LWR_Atmosphere_Yakutsk_Night';

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
ViolentViolet = [0.6 0 0.7];

%-------------------------------------------------------------------------%
%----------------------  calculate validation data  ----------------------%
% limit unrealistic sub-canopy SWR
SW_out_bc_98_1h_filled_corrected = nan(size(SW_out_bc_98_1h_filled));
for t=1:length(Rnet_bc_98_1h_filled)
    SW_out_bc_98_1h_filled_corrected(t) ...
        = min(SW_out_bc_98_1h_filled(t),SW_in_bc_98_1h_filled(t));
end
SW_net_bc_98_1h_filled_corrected = nan(size(Rnet_bc_98_1h_filled));
for l=1:length(SW_net_bc_98_1h_filled_corrected)
    SW_net_bc_98_1h_filled_corrected(l) ...
        = SW_in_bc_98_1h_filled(l) - SW_out_bc_98_1h_filled_corrected(l);
end

% calculate net sub-canopy LWR
LW_net_bc_98 = Rnet_bc_98_1h_filled - SW_net_bc_98_1h_filled_corrected;

% subtract outgoing LWR calculated from surface temperature
boltz = 5.67*10^(-8);
em_soil = 0.96; em_snow = 0.97;
em_gr_98 = nan(size(Rnet_bc_98_1h_filled));
em_gr_98(1:3195) = em_snow; em_gr_98(3196:end) = em_soil;
LW_in_bc_98 = nan(size(Rnet_bc_98_1h_filled));
for t=1:length(Rnet_bc_98_1h_filled)
    LW_in_bc_98(t) = LW_net_bc_98(t) + em_gr_98(t)*boltz*T_surf_98_1h_filled(t)^4;
end

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
EP_98_1 = 1057:1728;
EP_98_2 = 1777:3192;

LW_sub_val_eval = vertcat(LW_in_bc_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_98(EP_98_2(1):EP_98_2(end)));

% surface albedo
LW_sub_CLM_eval_AlbSurf = vertcat(LW_in_bc_CLM_AlbSurf(EP_98_1(1):EP_98_1(end),:),...
    LW_in_bc_CLM_AlbSurf(EP_98_2(1):EP_98_2(end),:));

% diffuse fraction
LW_sub_CLM_eval_Fdiff = vertcat(LW_in_bc_CLM_Fdiff(EP_98_1(1):EP_98_1(end),:),...
    LW_in_bc_CLM_Fdiff(EP_98_2(1):EP_98_2(end),:));

% snow depth
LW_sub_CLM_eval_Zsnow = vertcat(LW_in_bc_CLM_Zsnow(EP_98_1(1):EP_98_1(end),:),...
    LW_in_bc_CLM_Zsnow(EP_98_2(1):EP_98_2(end),:));

SWR_Incoming_Yakutsk = vertcat(SW_in_ac_98_1h_filled(EP_98_1(1):EP_98_1(end),:),...
    SW_in_ac_98_1h_filled(EP_98_2(1):EP_98_2(end),:));

i=0;
for l=1:length(SWR_Incoming_Yakutsk)
    if SWR_Incoming_Yakutsk(l) == 0
        i=i+1;
        LW_sub_val_eval_Night(i,:) = LW_sub_val_eval(l,:);
        LW_sub_CLM_eval_AlbSurf_Night(i,:) = LW_sub_CLM_eval_AlbSurf(l,:);
        LW_sub_CLM_eval_Fdiff_Night(i,:) = LW_sub_CLM_eval_Fdiff(l,:);
        LW_sub_CLM_eval_Zsnow_Night(i,:) = LW_sub_CLM_eval_Zsnow(l,:);
    end
end

%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%

% CLM
RMSE_CLM = nan(11,3); MBD_CLM = nan(11,3);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval_Night),...
        LW_sub_CLM_eval_AlbSurf_Night(:,j),LW_sub_val_eval_Night);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval_Night),...
        LW_sub_CLM_eval_AlbSurf_Night(:,j),LW_sub_val_eval_Night);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval_Night),...
        LW_sub_CLM_eval_Fdiff_Night(:,j),LW_sub_val_eval_Night);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval_Night),...
        LW_sub_CLM_eval_Fdiff_Night(:,j),LW_sub_val_eval_Night);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval_Night),...
        LW_sub_CLM_eval_Zsnow_Night(:,j),LW_sub_val_eval_Night);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval_Night),...
        LW_sub_CLM_eval_Zsnow_Night(:,j),LW_sub_val_eval_Night);
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
LWE_sub_val_eval = LW_sub_val_eval_Night./LWR_Atmosphere_Yakutsk_Night;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_AlbSurf_Night));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) ...
        = LW_sub_CLM_eval_AlbSurf_Night(:,j)./LWR_Atmosphere_Yakutsk_Night;
end

% diffuse fraction
LWE_sub_CLM_eval_Fdiff = nan(size(LW_sub_CLM_eval_Fdiff_Night));
for j=1:11
    LWE_sub_CLM_eval_Fdiff(:,j) ...
        = LW_sub_CLM_eval_Fdiff_Night(:,j)./LWR_Atmosphere_Yakutsk_Night;
end

% snow depth
LWE_sub_CLM_eval_Zsnow = nan(size(LW_sub_CLM_eval_Zsnow_Night));
for j=1:11
    LWE_sub_CLM_eval_Zsnow(:,j) ...
        = LW_sub_CLM_eval_Zsnow_Night(:,j)./LWR_Atmosphere_Yakutsk_Night;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_Yak = hist_obs.Values;
hist_clm_Yak_AlbSurf = nan(length(hist_obs_Yak),11);
hist_clm_Yak_Fdiff = nan(length(hist_obs_Yak),11);
hist_clm_Yak_Zsnow = nan(length(hist_obs_Yak),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Yak_AlbSurf(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Fdiff(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Yak_Fdiff(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Zsnow(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Yak_Zsnow(:,j) = hist_clm.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
hold on
text(0.8+0.05*1.2,0.95*0.2,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Yak,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Yak_AlbSurf,'Color',ViolentViolet)
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
text(0.8+0.05*1.2,0.95*0.2,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Yak,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Yak_Fdiff,'Color',ViolentViolet)
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
text(0.8+0.05*1.2,0.95*0.2,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Yak,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Yak_Zsnow,'Color',ViolentViolet)
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Yakutsk.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 0 33.5 10];
fig.PaperSize = [28.5 9.5];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Yakutsk.pdf')