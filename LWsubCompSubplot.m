tic

clear
close all

%-----------------------------  import data  -----------------------------%
%{
Create subplot for model comparison.
1) Run ToyModelAtAlptal.m without HM parameterization for CLM4.5.
2) Save sub-canopy longwave radiation:

EP_2004 = start_2004+18:1728;   % evaluation period
gap_2004_1 = 601:624;           % gap
gap_2004_2 = 1009:1032;         % gap
gap_2004_3 = 1201:1224;         % gap
gap_2004_4 = 1633:1728;         % gap
LW_sub_val_eval_Alp = vertcat(LW_in_bc_1h_all(EP_2004(1):gap_2004_1(1)-1,1),...
LW_in_bc_1h_all(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
LW_in_bc_1h_all(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
LW_in_bc_1h_all(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
LW_in_bc_1h_all(gap_2004_4(end)+1:EP_2004(end),1));
LW_sub_CLM_eval_Alp = vertcat(LW_in_bc_CLM(EP_2004(1):gap_2004_1(1)-1,1),...
LW_in_bc_CLM(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
LW_in_bc_CLM(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
LW_in_bc_CLM(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
LW_in_bc_CLM(gap_2004_4(end)+1:EP_2004(end),1));
LW_sub_SP_eval_Alp = vertcat(LW_in_bc_SP(EP_2004(1):gap_2004_1(1)-1,1),...
LW_in_bc_SP(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
LW_in_bc_SP(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
LW_in_bc_SP(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
LW_in_bc_SP(gap_2004_4(end)+1:EP_2004(end),1));
save('ModelComp_LWsub_Alptal.mat','LW_sub_val_eval_Alp','LW_sub_CLM_eval_Alp','LW_sub_SP_eval_Alp')

3) Run ToyModelAtSeehornwald.m without HM parameterization for CLM4.5.
4) Save sub-canopy longwave radiation:

EP_2008 = 1921:4728;        % evaluation period
EP_2009 = 10705:12864;      % evaluation period
EP_2010 = 19465:22080;      % evaluation period
EP_2011 = 28225:30312;      % evaluation period
EP_2012 = 36985:39768;      % evaluation period
gap_2010 = 1200:1224;       % gap due to met forcing
gap_2011 = 408:456;         % gap due to met forcing
gap_LWsub_2008 = 159;       % gap due to evaluation data
gap_LWsub_2010 = 134:137;   % gap due to evaluation data
gap_LWsub_2010_2 = 1200;    % gap due to evaluation data
gap_LWsub_2011 = 346:349;   % gap due to evaluation data
gap_LWsub_2011_2 = 427:442; % gap due to evaluation data
LW_sub_val_eval_SHW08 = vertcat(LW_in_bc_1h(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_1h(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_CLM_eval_SHW08 = vertcat(LW_in_bc_CLM(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_CLM(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_SP_eval_SHW08 = vertcat(LW_in_bc_SP(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_SP(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_val_eval_SHW09 = LW_in_bc_1h(EP_2009(1):EP_2009(end));
LW_sub_CLM_eval_SHW09 = LW_in_bc_CLM(EP_2009(1):EP_2009(end));
LW_sub_SP_eval_SHW09 = LW_in_bc_SP(EP_2009(1):EP_2009(end));
LW_sub_val_eval_SHW10 = vertcat(LW_in_bc_1h(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_1h(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_1h(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_CLM_eval_SHW10 = vertcat(LW_in_bc_CLM(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_CLM(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_CLM(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_SP_eval_SHW10 = vertcat(LW_in_bc_SP(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_SP(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_SP(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_val_eval_SHW11 = vertcat(LW_in_bc_1h(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_1h(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_1h(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_CLM_eval_SHW11 = vertcat(LW_in_bc_CLM(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
LW_in_bc_CLM(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
LW_in_bc_CLM(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_SP_eval_SHW11 = vertcat(LW_in_bc_SP(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_SP(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_SP(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_val_eval_SHW12 = LW_in_bc_1h(EP_2012(1):EP_2012(end));
LW_sub_CLM_eval_SHW12 = LW_in_bc_CLM(EP_2012(1):EP_2012(end));
LW_sub_SP_eval_SHW12 = LW_in_bc_SP(EP_2012(1):EP_2012(end));
save('ModelComp_LWsub_Seehornwald.mat','LW_sub_val_eval_SHW08','LW_sub_val_eval_SHW09',...
    'LW_sub_val_eval_SHW10','LW_sub_val_eval_SHW11','LW_sub_val_eval_SHW12',...
    'LW_sub_CLM_eval_SHW08','LW_sub_CLM_eval_SHW09','LW_sub_CLM_eval_SHW10',...
    'LW_sub_CLM_eval_SHW11','LW_sub_CLM_eval_SHW12','LW_sub_SP_eval_SHW08',...
    'LW_sub_SP_eval_SHW09','LW_sub_SP_eval_SHW10','LW_sub_SP_eval_SHW11','LW_sub_SP_eval_SHW12')

5) Run ToyModelAtSodankyla.m for site 'C' without HM parameterization for CLM4.5.
6) Save sub-canopy longwave radiation:

EP_C_1 = 12:443;
EP_C_2 = 468:923;
LW_sub_val_eval = vertcat(LW_in_bc_C_1h(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_C_1h(EP_C_2(1):EP_C_2(end),:));
LW_sub_val_eval_Sod = nan(length(LW_sub_val_eval(:,1)),1);
for l=1:length(LW_sub_val_eval_Sod)
    LW_sub_val_eval_Sod(l) = mean(LW_sub_val_eval(l,:));
end
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_CLM(EP_C_2(1):EP_C_2(end),:));
LW_sub_CLM_eval_Sod = nan(length(LW_sub_CLM_eval(:,1)),1);
for l=1:length(LW_sub_CLM_eval_Sod)
    LW_sub_CLM_eval_Sod(l) = mean(LW_sub_CLM_eval(l,:));
end
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_SP(EP_C_2(1):EP_C_2(end),:));
LW_sub_SP_eval_Sod = nan(length(LW_sub_SP_eval(:,1)),1);
for l=1:length(LW_sub_SP_eval_Sod)
    LW_sub_SP_eval_Sod(l) = mean(LW_sub_SP_eval(l,:));
end
save('ModelComp_LWsub_Sodankyla.mat','LW_sub_val_eval_Sod','LW_sub_CLM_eval_Sod','LW_sub_SP_eval_Sod')
%}
load ModelComp_LWsub_Alptal.mat
load ModelComp_LWsub_Seehornwald.mat
load ModelComp_LWsub_Sodankyla.mat

%---------------------------  plot essentials  ---------------------------%
GlacialGrey = [0.4 0.7 1];
MagicMaroon = [0.65 0.32 0.35];
EcstaticEmerald = [0.17 0.52 0.5];
CandidCoral = [1 0.44 0.32];

%-------------------------  calculate RMSE & MB  -------------------------%
% Alptal
RMSE_CLM_Alp = RMSE(length(LW_sub_val_eval_Alp),LW_sub_CLM_eval_Alp,LW_sub_val_eval_Alp);
MB_CLM_Alp = MBD(length(LW_sub_val_eval_Alp),LW_sub_CLM_eval_Alp,LW_sub_val_eval_Alp);
RMSE_SP_Alp = RMSE(length(LW_sub_val_eval_Alp),LW_sub_SP_eval_Alp,LW_sub_val_eval_Alp);
MB_SP_Alp = MBD(length(LW_sub_val_eval_Alp),LW_sub_SP_eval_Alp,LW_sub_val_eval_Alp);

% Seehornwald
RMSE_CLM_SHW = RMSE(length(LW_sub_val_eval_SHW),LW_sub_CLM_eval_SHW,LW_sub_val_eval_SHW);
MB_CLM_SHW = MBD(length(LW_sub_val_eval_SHW),LW_sub_CLM_eval_SHW,LW_sub_val_eval_SHW);
RMSE_SP_SHW = RMSE(length(LW_sub_val_eval_SHW),LW_sub_SP_eval_SHW,LW_sub_val_eval_SHW);
MB_SP_SHW = MBD(length(LW_sub_val_eval_SHW),LW_sub_SP_eval_SHW,LW_sub_val_eval_SHW);

% Sodankyla
RMSE_CLM_Sod = RMSE(length(LW_sub_val_eval_Sod),LW_sub_CLM_eval_Sod,LW_sub_val_eval_Sod);
MB_CLM_Sod = MBD(length(LW_sub_val_eval_Sod),LW_sub_CLM_eval_Sod,LW_sub_val_eval_Sod);
RMSE_SP_Sod = RMSE(length(LW_sub_val_eval_Sod),LW_sub_SP_eval_Sod,LW_sub_val_eval_Sod);
MB_SP_Sod = MBD(length(LW_sub_val_eval_Sod),LW_sub_SP_eval_Sod,LW_sub_val_eval_Sod);

%-----------------------------  create plot  -----------------------------%
% sub-canopy LW comparison
fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
subplot(1,3,1)
hold on
text(160,435,'a','FontSize',13,'FontWeight','bold')
plot([150 450],[150 450],'k')
plot(LW_sub_val_eval_Alp,LW_sub_CLM_eval_Alp,'Color',EcstaticEmerald,...
    'Marker','.','LineStyle','none')
plot(LW_sub_val_eval_Alp,LW_sub_SP_eval_Alp,'Color',CandidCoral,...
    'Marker','.','LineStyle','none')
text(440,230,'MB [W m^{-2}]:','FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
str1 = num2str(MB_CLM_Alp,3);
str2 = num2str(MB_SP_Alp,3);
text(440,205,['\color[rgb]{0.17 0.52 0.5}',str1,'\color{black}, ',...
    '\color[rgb]{1 0.44 0.32}',str2],'FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
text(440,180,'RMSE [W m^{-2}]:','FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
str1 = num2str(RMSE_CLM_Alp,4);
str2 = num2str(RMSE_SP_Alp,3);
text(440,155,['\color[rgb]{0.17 0.52 0.5}',str1,'\color{black}, ',...
    '\color[rgb]{1 0.44 0.32}',str2],'FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
xlim([150 450])
ylim([150 450])
xlabel('observed sub-canopy LWR [W m^{-2}]','FontSize',13,'FontWeight','bold')
ylabel('simulated sub-canopy LWR [W m^{-2}]','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
box on
pbaspect([1 1 1])
subplot(1,3,2)
hold on
text(160,435,'b','FontSize',13,'FontWeight','bold')
plot([150 450],[150 450],'k')
plot(LW_sub_val_eval_SHW,LW_sub_CLM_eval_SHW,'Color',MagicMaroon,...
    'Marker','.','LineStyle','none')
plot(LW_sub_val_eval_SHW,LW_sub_SP_eval_SHW,'Color',CandidCoral,...
    'Marker','.','LineStyle','none')
text(440,230,'MB [W m^{-2}]:','FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
str1 = num2str(MB_CLM_SHW,3);
str2 = num2str(MB_SP_SHW,3);
text(440,205,['\color[rgb]{0.65 0.32 0.35}',str1,'\color{black}, ',...
    '\color[rgb]{1 0.44 0.32}',str2],'FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
text(440,180,'RMSE [W m^{-2}]:','FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
str1 = num2str(RMSE_CLM_SHW,4);
str2 = num2str(RMSE_SP_SHW,4);
text(440,155,['\color[rgb]{0.65 0.32 0.35}',str1,'\color{black}, ',...
    '\color[rgb]{1 0.44 0.32}',str2],'FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
xlim([150 450])
ylim([150 450])
xlabel('observed sub-canopy LWR [W m^{-2}]','FontSize',13,'FontWeight','bold')
ylabel('simulated sub-canopy LWR [W m^{-2}]','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
box on
pbaspect([1 1 1])
subplot(1,3,3)
hold on
text(160,435,'c','FontSize',13,'FontWeight','bold')
plot([150 450],[150 450],'k')
plot(LW_sub_val_eval_Sod,LW_sub_SP_eval_Sod,'Color',CandidCoral,...
    'Marker','.','LineStyle','none')
plot(LW_sub_val_eval_Sod,LW_sub_CLM_eval_Sod,'Color',GlacialGrey,...
    'Marker','.','LineStyle','none')
text(440,230,'MB [W m^{-2}]:','FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
str1 = num2str(MB_CLM_Sod,2);
str2 = num2str(MB_SP_Sod,4);
text(440,205,['\color[rgb]{0.4 0.7 1}',str1,'\color{black}, ',...
    '\color[rgb]{1 0.44 0.32}',str2],'FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
text(440,180,'RMSE [W m^{-2}]:','FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
str1 = num2str(RMSE_CLM_Sod,4);
str2 = num2str(RMSE_SP_Sod,4);
text(440,155,['\color[rgb]{0.4 0.7 1}',str1,'\color{black}, ',...
    '\color[rgb]{1 0.44 0.32}',str2],'FontSize',13,'FontWeight','bold',...
    'HorizontalAlignment', 'right','VerticalAlignment', 'bottom')
xlim([150 450])
ylim([150 450])
xlabel('observed sub-canopy LWR [W m^{-2}]','FontSize',13,'FontWeight','bold')
ylabel('simulated sub-canopy LWR [W m^{-2}]','FontSize',13,'FontWeight','bold')
set(gca,'FontSize',13,'FontWeight','bold','LineWidth',2)
box on
pbaspect([1 1 1])
fig.PaperUnits = 'inches';
fig.PaperPosition = [-2.5 0 33.5 10];
fig.PaperSize = [28.5 9.5];
print(fig,'-dpng','-r600','Overview_LWsub_Scatter.png')
print(fig,'-dpdf','-r600','Overview_LWsub_Scatter.pdf')

toc