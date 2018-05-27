clear
close all

% deciduous -> no second layer -> no SNOWPACK-2LHM simulation
load LWsub_Cherskiy_AlbSurf.mat
LW_in_bc_CLM_AlbSurf = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_AlbSurf = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Cherskiy_Fdiff.mat
LW_in_bc_CLM_Fdiff = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_Fdiff = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Cherskiy_Ffrozen.mat
LW_in_bc_CLM_Ffrozen = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_Ffrozen = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Cherskiy_SAI.mat
LW_in_bc_CLM_SAI = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_SAI = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Cherskiy_SWC.mat
LW_in_bc_CLM_SWC = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_SWC = LW_in_bc_SP;
clear LW_in_bc_SP
load LWsub_Cherskiy_Zsnow.mat
LW_in_bc_CLM_Zsnow = LW_in_bc_CLM;
clear LW_in_bc_CLM
% LW_in_bc_SP_Zsnow = LW_in_bc_SP;
clear LW_in_bc_SP
% load atmospheric LWR for LWE
load('MLRdata_Cherskiy_HM.mat','LWR_Atmosphere_Cherskiy')
LWR_Atmosphere_Cherskiy = LWR_Atmosphere_Cherskiy';
LW_in_bc_1h = LW_in_bc_1h';

%-------------------------------------------------------------------------%
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
LightGrey = [0.65 0.65 0.65];
BastilleBleu = [0.12 0.34 0.56];

%-------------------------------------------------------------------------%
%----------------  select evaluation periods based on QC  ----------------%
EP1 = 12:1163;
EP2 = 1188:length(LW_in_bc_1h);

LW_sub_val_eval = vertcat(LW_in_bc_1h(EP1(1):EP1(end)),...
    LW_in_bc_1h(EP2(1):EP2(end)));

% surface albedo
LW_sub_CLM_eval_AlbSurf = vertcat(LW_in_bc_CLM_AlbSurf(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_AlbSurf(EP2(1):EP2(end),:));

% diffuse fraction
LW_sub_CLM_eval_Fdiff = vertcat(LW_in_bc_CLM_Fdiff(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_Fdiff(EP2(1):EP2(end),:));

% fraction of frozen soil
LW_sub_CLM_eval_Ffrozen = vertcat(LW_in_bc_CLM_Ffrozen(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_Ffrozen(EP2(1):EP2(end),:));

% SAI
LW_sub_CLM_eval_SAI = vertcat(LW_in_bc_CLM_SAI(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_SAI(EP2(1):EP2(end),:));

% soil water content
LW_sub_CLM_eval_SWC = vertcat(LW_in_bc_CLM_SWC(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_SWC(EP2(1):EP2(end),:));

% snow depth
LW_sub_CLM_eval_Zsnow = vertcat(LW_in_bc_CLM_Zsnow(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_Zsnow(EP2(1):EP2(end),:));


%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%

% CLM
RMSE_CLM = nan(11,6); MBD_CLM = nan(11,6);
for j=1:11
    RMSE_CLM(j,1) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_AlbSurf(:,j),...
        LW_sub_val_eval);
    MBD_CLM(j,1) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_AlbSurf(:,j),...
        LW_sub_val_eval);
    RMSE_CLM(j,2) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_Fdiff(:,j),...
        LW_sub_val_eval);
    MBD_CLM(j,2) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_Fdiff(:,j),...
        LW_sub_val_eval);
    RMSE_CLM(j,3) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_Ffrozen(:,j),...
        LW_sub_val_eval);
    MBD_CLM(j,3) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_Ffrozen(:,j),...
        LW_sub_val_eval);
    RMSE_CLM(j,4) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_SAI(:,j),...
        LW_sub_val_eval);
    MBD_CLM(j,4) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_SAI(:,j),...
        LW_sub_val_eval);
    RMSE_CLM(j,5) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_SWC(:,j),...
        LW_sub_val_eval);
    MBD_CLM(j,5) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_SWC(:,j),...
        LW_sub_val_eval);
    RMSE_CLM(j,6) = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval_Zsnow(:,j),...
        LW_sub_val_eval);
    MBD_CLM(j,6) = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval_Zsnow(:,j),...
        LW_sub_val_eval);
end
RMSE_CLM_change = nan(size(RMSE_CLM));
for i=1:6
    RMSE_CLM_change(:,i) = RMSE_CLM(:,i) - RMSE_CLM(6,i);
end
MBD_CLM_change = nan(size(MBD_CLM));
for i=1:6
    MBD_CLM_change(:,i) = MBD_CLM(:,i) - MBD_CLM(6,i);
end

%------------------------  longwave  enhancement  ------------------------%
LWE_sub_val_eval = LW_sub_val_eval./LWR_Atmosphere_Cherskiy;

% surface albedo
LWE_sub_CLM_eval_AlbSurf = nan(size(LW_sub_CLM_eval_AlbSurf));
for j=1:11
    LWE_sub_CLM_eval_AlbSurf(:,j) = LW_sub_CLM_eval_AlbSurf(:,j)...
        ./LWR_Atmosphere_Cherskiy;
end

% diffuse fraction
LWE_sub_CLM_eval_Fdiff = nan(size(LW_sub_CLM_eval_Fdiff));
for j=1:11
    LWE_sub_CLM_eval_Fdiff(:,j) = LW_sub_CLM_eval_Fdiff(:,j)...
        ./LWR_Atmosphere_Cherskiy;
end

% fraction of frozen soil
LWE_sub_CLM_eval_Ffrozen = nan(size(LW_sub_CLM_eval_Ffrozen));
for j=1:11
    LWE_sub_CLM_eval_Ffrozen(:,j) = LW_sub_CLM_eval_Ffrozen(:,j)...
        ./LWR_Atmosphere_Cherskiy;
end

% SAI
LWE_sub_CLM_eval_SAI = nan(size(LW_sub_CLM_eval_SAI));
for j=1:11
    LWE_sub_CLM_eval_SAI(:,j) = LW_sub_CLM_eval_SAI(:,j)...
        ./LWR_Atmosphere_Cherskiy;
end

% soil water content
LWE_sub_CLM_eval_SWC = nan(size(LW_sub_CLM_eval_SWC));
for j=1:11
    LWE_sub_CLM_eval_SWC(:,j) = LW_sub_CLM_eval_SWC(:,j)...
        ./LWR_Atmosphere_Cherskiy;
end

% snow depth
LWE_sub_CLM_eval_Zsnow = nan(size(LW_sub_CLM_eval_Zsnow));
for j=1:11
    LWE_sub_CLM_eval_Zsnow(:,j) = LW_sub_CLM_eval_Zsnow(:,j)...
        ./LWR_Atmosphere_Cherskiy;
end

spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_sub_val_eval,spectrum_LWenh,'Normalization','probability');
    hist_obs_Che = hist_obs.Values;
hist_clm_Che_AlbSurf = nan(length(hist_obs_Che),11);
hist_clm_Che_Fdiff = nan(length(hist_obs_Che),11);
hist_clm_Che_Ffrozen = nan(length(hist_obs_Che),11);
hist_clm_Che_SAI = nan(length(hist_obs_Che),11);
hist_clm_Che_SWC = nan(length(hist_obs_Che),11);
hist_clm_Che_Zsnow = nan(length(hist_obs_Che),11);
for j=1:11
    hist_clm = histogram(LWE_sub_CLM_eval_AlbSurf(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Che_AlbSurf(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Fdiff(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Che_Fdiff(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Ffrozen(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Che_Ffrozen(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_SAI(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Che_SAI(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_SWC(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Che_SWC(:,j) = hist_clm.Values;
    hist_clm = histogram(LWE_sub_CLM_eval_Zsnow(:,j),spectrum_LWenh,...
        'Normalization','probability');
    hist_clm_Che_Zsnow(:,j) = hist_clm.Values;
end

fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(2,3,1)
hold on
text(0.8+0.05*1.2,0.95*0.3,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che_AlbSurf,'Color',BastilleBleu)
hold off
xlim([0.8 2])
ylim([0 0.3])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,3,2)
hold on
text(0.8+0.05*1.2,0.95*0.3,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che_Fdiff,'Color',BastilleBleu)
hold off
xlim([0.8 2])
ylim([0 0.3])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,3,3)
hold on
text(0.8+0.05*1.2,0.95*0.3,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che_Ffrozen,'Color',BastilleBleu)
hold off
xlim([0.8 2])
ylim([0 0.3])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,3,4)
hold on
text(0.8+0.05*1.2,0.95*0.3,'d','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che_SAI,'Color',BastilleBleu)
hold off
xlim([0.8 2])
ylim([0 0.3])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,3,5)
hold on
text(0.8+0.05*1.2,0.95*0.3,'e','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che_SWC,'Color',BastilleBleu)
hold off
xlim([0.8 2])
ylim([0 0.3])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
subplot(2,3,6)
hold on
text(0.8+0.05*1.2,0.95*0.3,'f','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che_Zsnow,'Color',BastilleBleu)
hold off
xlim([0.8 2])
ylim([0 0.3])
xlabel('LW enhancement','FontSize',13,'FontWeight','bold')
ylabel('Frequency','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
pbaspect([1 1 1])
print(fig,'-dpng','-r600','PDF_Sensitivity_Cherskiy.png')
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 -0.5 33.5 20.5];
fig.PaperSize = [28.5 19];
print(fig,'-dpdf','-r600','PDF_Sensitivity_Cherskiy.pdf')