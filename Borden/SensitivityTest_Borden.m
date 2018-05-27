clear
close all

% deciduous -> no second layer -> no SNOWPACK-2LHM simulation
load LWsub_Borden_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM_comb;
clear LW_in_bc_CLM
% LW_in_bc_SP_AlbSurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Borden_Fdiff.mat
LW_in_bc_CLM_Fdiff = LW_in_bc_CLM_comb;
clear LW_in_bc_CLM
% LW_in_bc_SP_Fdiff = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Borden_fPFT.mat
LW_in_bc_CLM_fPFT = LW_in_bc_CLM_comb;
clear LW_in_bc_CLM
% LW_in_bc_SP_fPFT = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Borden_Zsnow.mat
LW_in_bc_CLM_Zsnow = LW_in_bc_CLM_comb;
clear LW_in_bc_CLM
% LW_in_bc_SP_Zsnow = LW_in_bc_SP;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Borden.mat','LWR_Atmosphere_Borden')
LWR_Atmosphere_Borden = LWR_Atmosphere_Borden';
LW_in_bc_1h = LW_in_bc_1h';

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
RussianRed = [1 0.032 0.064];

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
EP_start = 25;
EP_end = 2232;

gap_LW_1 = 649:696;     % 28 Jan 1:00 to 30 Jan 0:00
gap_LW_2 = 985:1008;    % 11 Feb 1:00 to 12 Feb 0:00
gap_LW_3 = 1057:1080;   % 14 Feb 1:00 to 15 Feb 0:00
gap_LW_4 = 1273:1296;   % 23 Feb 1:00 to 24 Feb 0:00
gap_LW_7 = 1585:1608;   % 8 Mar 1:00 to 9 Mar 0:00
gap_LW_8 = 2425:2448;   % 12 Apr 1:00 to 13 Apr 0:00
gap_snow_1 = 1321:1416;   % 25 Feb 1:00 to 28 Feb 0:00
gap_snow_2 = 1969:1992;   % 24 Mar 1:00 to 25 Mar 0:00
gap_snow_3 = 2137:2160;   % 31 Mar 1:00 to 1 Apr 0:00
gap_met = 721:768;

LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_start:gap_LW_1(1)-1),...
    LW_in_bc_1h(gap_LW_1(end)+1:gap_met(1)-1),...
    LW_in_bc_1h(gap_met(end)+1:gap_LW_2(1)-1),...
    LW_in_bc_1h(gap_LW_2(end)+1:gap_LW_3(1)-1),...
    LW_in_bc_1h(gap_LW_3(end)+1:gap_LW_4(1)-1),...
    LW_in_bc_1h(gap_LW_4(end)+1:gap_snow_1(1)-1),...
    LW_in_bc_1h(gap_snow_1(end)+1:gap_LW_7(1)-1),...
    LW_in_bc_1h(gap_LW_7(end)+1:gap_snow_2(1)-1),...
    LW_in_bc_1h(gap_snow_2(end)+1:gap_snow_3(1)-1),...
    LW_in_bc_1h(gap_snow_3(end)+1:gap_LW_8(1)-1),...
    LW_in_bc_1h(gap_LW_8(end)+1:EP_end));

% surface albedo
LW_sub_CLM_eval_AlbSurf = vertcat(LW_in_bc_CLM_AlbSurf(EP_start:gap_LW_1(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_LW_1(end)+1:gap_met(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_met(end)+1:gap_LW_2(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_LW_2(end)+1:gap_LW_3(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_LW_3(end)+1:gap_LW_4(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_LW_4(end)+1:gap_snow_1(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_snow_1(end)+1:gap_LW_7(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_LW_7(end)+1:gap_snow_2(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_snow_2(end)+1:gap_snow_3(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_snow_3(end)+1:gap_LW_8(1)-1,:),...
    LW_in_bc_CLM_AlbSurf(gap_LW_8(end)+1:EP_end,:));

% diffuse fraction
LW_sub_CLM_eval_Fdiff = vertcat(LW_in_bc_CLM_Fdiff(EP_start:gap_LW_1(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_LW_1(end)+1:gap_met(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_met(end)+1:gap_LW_2(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_LW_2(end)+1:gap_LW_3(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_LW_3(end)+1:gap_LW_4(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_LW_4(end)+1:gap_snow_1(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_snow_1(end)+1:gap_LW_7(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_LW_7(end)+1:gap_snow_2(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_snow_2(end)+1:gap_snow_3(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_snow_3(end)+1:gap_LW_8(1)-1,:),...
    LW_in_bc_CLM_Fdiff(gap_LW_8(end)+1:EP_end,:));

% PFT fraction
LW_sub_CLM_eval_fPFT = vertcat(LW_in_bc_CLM_fPFT(EP_start:gap_LW_1(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_LW_1(end)+1:gap_met(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_met(end)+1:gap_LW_2(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_LW_2(end)+1:gap_LW_3(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_LW_3(end)+1:gap_LW_4(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_LW_4(end)+1:gap_snow_1(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_snow_1(end)+1:gap_LW_7(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_LW_7(end)+1:gap_snow_2(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_snow_2(end)+1:gap_snow_3(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_snow_3(end)+1:gap_LW_8(1)-1,:),...
    LW_in_bc_CLM_fPFT(gap_LW_8(end)+1:EP_end,:));

% snow depth
LW_sub_CLM_eval_Zsnow = vertcat(LW_in_bc_CLM_Zsnow(EP_start:gap_LW_1(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_LW_1(end)+1:gap_met(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_met(end)+1:gap_LW_2(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_LW_2(end)+1:gap_LW_3(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_LW_3(end)+1:gap_LW_4(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_LW_4(end)+1:gap_snow_1(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_snow_1(end)+1:gap_LW_7(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_LW_7(end)+1:gap_snow_2(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_snow_2(end)+1:gap_snow_3(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_snow_3(end)+1:gap_LW_8(1)-1,:),...
    LW_in_bc_CLM_Zsnow(gap_LW_8(end)+1:EP_end,:));

%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%
% CLM
RMSE_CLM = nan(11,4); MBD_CLM = nan(11,4);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_AlbSurf(:,j),LW_sub_val_eval);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_AlbSurf(:,j),LW_sub_val_eval);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_Fdiff(:,j),LW_sub_val_eval);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_Fdiff(:,j),LW_sub_val_eval);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_fPFT(:,j),LW_sub_val_eval);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_fPFT(:,j),LW_sub_val_eval);
    RMSE_CLM(j,4) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_Zsnow(:,j),LW_sub_val_eval);
    MBD_CLM(j,4) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_Zsnow(:,j),LW_sub_val_eval);
end
RMSE_CLM_change = nan(size(RMSE_CLM));
for i=1:4
    RMSE_CLM_change(:,i) = RMSE_CLM(:,i) - RMSE_CLM(6,i);
end
MBD_CLM_change = nan(size(MBD_CLM));
for i=1:4
    MBD_CLM_change(:,i) = MBD_CLM(:,i) - MBD_CLM(6,i);
end

%------------------------  longwave  enhancement  ------------------------%
LWE_sub_val_eval = LW_sub_val_eval./LWR_Atmosphere_Borden;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_AlbSurf));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) = LW_sub_CLM_eval_AlbSurf(:,j)./LWR_Atmosphere_Borden;
end

% diffuse fraction
LWE_sub_CLM_eval_Fdiff = nan(size(LW_sub_CLM_eval_Fdiff));
for j=1:11
    LWE_sub_CLM_eval_Fdiff(:,j) = LW_sub_CLM_eval_Fdiff(:,j)./LWR_Atmosphere_Borden;
end

% PFT fraction
LWE_sub_CLM_eval_fPFT = nan(size(LW_sub_CLM_eval_fPFT));
for j=1:11
    LWE_sub_CLM_eval_fPFT(:,j) = LW_sub_CLM_eval_fPFT(:,j)./LWR_Atmosphere_Borden;
end

% snow depth
LWE_sub_CLM_eval_Zsnow = nan(size(LW_sub_CLM_eval_Zsnow));
for j=1:11
    LWE_sub_CLM_eval_Zsnow(:,j) = LW_sub_CLM_eval_Zsnow(:,j)./LWR_Atmosphere_Borden;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_Bor = hist_obs.Values;
hist_clm_Bor_AlbSurf = nan(length(hist_obs_Bor),11);
hist_clm_Bor_Fdiff = nan(length(hist_obs_Bor),11);
hist_clm_Bor_fPFT = nan(length(hist_obs_Bor),11);
hist_clm_Bor_Zsnow = nan(length(hist_obs_Bor),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Bor_AlbSurf(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Fdiff(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Bor_Fdiff(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_fPFT(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Bor_fPFT(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Zsnow(:,j),spectrum_LWenh,'Normalization','probability');
    hist_clm_Bor_Zsnow(:,j) = hist_clm.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(2,2,1)
hold on
text(0.8+0.05*1.2,0.95*0.4,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Bor,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Bor_AlbSurf,'Color',RussianRed)
hold off
xlim([0.8 2])
ylim([0 0.4])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,2,2)
hold on
text(0.8+0.05*1.2,0.95*0.4,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Bor,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Bor_Fdiff,'Color',RussianRed)
hold off
xlim([0.8 2])
ylim([0 0.4])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,2,3)
hold on
text(0.8+0.05*1.2,0.95*0.4,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Bor,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Bor_fPFT,'Color',RussianRed)
hold off
xlim([0.8 2])
ylim([0 0.4])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,2,4)
hold on
text(0.8+0.05*1.2,0.95*0.4,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Bor,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Bor_Zsnow,'Color',RussianRed)
hold off
xlim([0.8 2])
ylim([0 0.4])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Borden.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.5 -0.5 20.5 20.5];
fig.PaperSize = [19 19];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Borden.pdf')