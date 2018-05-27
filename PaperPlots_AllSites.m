tic

clear
close all
%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end
GlacialGrey = [0.4 0.7 1];
MagicMaroon = [0.65 0.32 0.35];
EcstaticEmerald = [0.17 0.52 0.5];
GorgeousGold = [1 0.85 0];
BastilleBleu = [0.12 0.34 0.56];
RussianRed = [1 0.032 0.064];
ViolentViolet = [0.6 0 0.7];


%-------------------------------------------------------------------------%
%-----------------------------  import data  -----------------------------%
%-------------------------------------------------------------------------%
load MLRdata_Abisko_HM.mat
load MLRdata_Abisko.mat
load MLRdata_Alptal_HM.mat
load MLRdata_Alptal.mat
load MLRdata_Borden_HM.mat
load MLRdata_Borden.mat
load MLRdata_Cherskiy_HM.mat
load MLRdata_Cherskiy.mat
load MLRdata_Seehornwald_HM.mat
load MLRdata_Seehornwald.mat
load MLRdata_Sodankyla_HM.mat
load MLRdata_Sodankyla.mat
load MLRdata_Yakutsk_HM.mat
load MLRdata_Yakutsk.mat

%-------------------------------------------------------------------------%
%---------------------  subplots of Toy Model sites  ---------------------%
%-------------------------------------------------------------------------%
%{
As Alptal and Seehornwald feature rail measurements and thus a spatial 
average, we use radiometer averages for Abisko and Sodankyla.
%}
    % Abisko
LWR_Subcanopy_Observations_Abisko_avg = nan(size(LWR_Atmosphere_Abisko));
LWR_Subcanopy_CLM_Abisko_avg = nan(size(LWR_Atmosphere_Abisko));
LWR_Subcanopy_CLM_Abisko_avg_HM = nan(size(LWR_Atmosphere_Abisko));
for l=1:length(LWR_Atmosphere_Abisko)
    LWR_Subcanopy_Observations_Abisko_avg(l) = ...
        mean(LWR_Subcanopy_Observations_Abisko(l,:));
    LWR_Subcanopy_CLM_Abisko_avg(l) = ...
        mean(LWR_Subcanopy_CLM_Abisko(l,:));
    LWR_Subcanopy_CLM_Abisko_avg_HM(l) = ...
        mean(LWR_Subcanopy_CLM_Abisko_HM(l,:));
end
    % Sodankyla
LWR_Subcanopy_Observations_Sodankyla_avg = nan(size(LWR_Atmosphere_Sodankyla));
LWR_Subcanopy_CLM_Sodankyla_avg = nan(size(LWR_Atmosphere_Sodankyla));
LWR_Subcanopy_CLM_Sodankyla_avg_HM = nan(size(LWR_Atmosphere_Sodankyla));
for l=1:length(LWR_Atmosphere_Sodankyla)
    LWR_Subcanopy_Observations_Sodankyla_avg(l) = ...
        mean(LWR_Subcanopy_Observations_Sodankyla(l,:));
    LWR_Subcanopy_CLM_Sodankyla_avg(l) = ...
        mean(LWR_Subcanopy_CLM_Sodankyla(l,:));
    LWR_Subcanopy_CLM_Sodankyla_avg_HM(l) = ...
        mean(LWR_Subcanopy_CLM_Sodankyla_HM(l,:));
end

%{
For Alptal and Seehornwald, seasons are aggregated.
%}
    % Alptal
SWR_Incoming_Alptal = vertcat(SWR_Incoming_Alptal_2004,...
    SWR_Incoming_Alptal_2005,SWR_Incoming_Alptal_2006,...
    SWR_Incoming_Alptal_2007);
Cosine_Solar_Angle_Alptal = vertcat(Cosine_Solar_Angle_Alptal_2004,...
    Cosine_Solar_Angle_Alptal_2005,Cosine_Solar_Angle_Alptal_2006,...
    Cosine_Solar_Angle_Alptal_2007);
Sky_Emissivity_Alptal = vertcat(Sky_Emissivity_Alptal_2004,...
    Sky_Emissivity_Alptal_2005,Sky_Emissivity_Alptal_2006,...
    Sky_Emissivity_Alptal_2007);
LWR_Atmosphere_Alptal = vertcat(LWR_Atmosphere_Alptal_2004,...
    LWR_Atmosphere_Alptal_2005,LWR_Atmosphere_Alptal_2006,...
    LWR_Atmosphere_Alptal_2007);
LWR_Subcanopy_Observations_Alptal = vertcat(LWR_Subcanopy_Observations_Alptal_2004,...
    LWR_Subcanopy_Observations_Alptal_2005,LWR_Subcanopy_Observations_Alptal_2006,...
    LWR_Subcanopy_Observations_Alptal_2007);
LWR_Subcanopy_CLM_Alptal = vertcat(LWR_Subcanopy_CLM_Alptal_2004,...
    LWR_Subcanopy_CLM_Alptal_2005,LWR_Subcanopy_CLM_Alptal_2006,...
    LWR_Subcanopy_CLM_Alptal_2007);
LWR_Subcanopy_CLM_Alptal_HM = vertcat(LWR_Subcanopy_CLM_Alptal_2004_HM,...
    LWR_Subcanopy_CLM_Alptal_2005_HM,LWR_Subcanopy_CLM_Alptal_2006_HM,...
    LWR_Subcanopy_CLM_Alptal_2007_HM);
    % Seehornwald
SWR_Incoming_Seehornwald = vertcat(SWR_Incoming_Seehornwald_2008,...
    SWR_Incoming_Seehornwald_2009,SWR_Incoming_Seehornwald_2010,...
    SWR_Incoming_Seehornwald_2011,SWR_Incoming_Seehornwald_2012);
Cosine_Solar_Angle_Seehornwald = vertcat(Cosine_Solar_Angle_Seehornwald_2008,...
    Cosine_Solar_Angle_Seehornwald_2009,Cosine_Solar_Angle_Seehornwald_2010,...
    Cosine_Solar_Angle_Seehornwald_2011,Cosine_Solar_Angle_Seehornwald_2012);
Sky_Emissivity_Seehornwald = vertcat(Sky_Emissivity_Seehornwald_2008,...
    Sky_Emissivity_Seehornwald_2009,Sky_Emissivity_Seehornwald_2010,...
    Sky_Emissivity_Seehornwald_2011,Sky_Emissivity_Seehornwald_2012);
LWR_Atmosphere_Seehornwald = vertcat(LWR_Atmosphere_Seehornwald_2008,...
    LWR_Atmosphere_Seehornwald_2009,LWR_Atmosphere_Seehornwald_2010,...
    LWR_Atmosphere_Seehornwald_2011,LWR_Atmosphere_Seehornwald_2012);
LWR_Subcanopy_Observations_Seehornwald = vertcat(LWR_Subcanopy_Observations_Seehornwald_2008,...
    LWR_Subcanopy_Observations_Seehornwald_2009,LWR_Subcanopy_Observations_Seehornwald_2010,...
    LWR_Subcanopy_Observations_Seehornwald_2011,LWR_Subcanopy_Observations_Seehornwald_2012);
LWR_Subcanopy_CLM_Seehornwald = vertcat(LWR_Subcanopy_CLM_Seehornwald_2008,...
    LWR_Subcanopy_CLM_Seehornwald_2009,LWR_Subcanopy_CLM_Seehornwald_2010,...
    LWR_Subcanopy_CLM_Seehornwald_2011,LWR_Subcanopy_CLM_Seehornwald_2012);
LWR_Subcanopy_CLM_Seehornwald_HM = vertcat(LWR_Subcanopy_CLM_Seehornwald_2008_HM,...
    LWR_Subcanopy_CLM_Seehornwald_2009_HM,LWR_Subcanopy_CLM_Seehornwald_2010_HM,...
    LWR_Subcanopy_CLM_Seehornwald_2011_HM,LWR_Subcanopy_CLM_Seehornwald_2012_HM);


%-------------------------  scattered LWR error  -------------------------%
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
% Alptal
sp1 = subplot(4,2,1);
hold on
text(0.4+0.05*0.8,25-0.05*50,'a','FontSize',12,'FontWeight','bold')
bla = (LWR_Subcanopy_CLM_Alptal - LWR_Subcanopy_Observations_Alptal)...
    ./LWR_Subcanopy_Observations_Alptal*100;
scatter(Sky_Emissivity_Alptal,bla,13,SWR_Incoming_Alptal)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
% xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Seehornwald
sp2 = subplot(4,2,2);
hold on
text(0.4+0.05*0.8,25-0.05*50,'b','FontSize',12,'FontWeight','bold')
bla = (LWR_Subcanopy_CLM_Seehornwald - LWR_Subcanopy_Observations_Seehornwald)...
    ./LWR_Subcanopy_Observations_Seehornwald*100;
scatter(Sky_Emissivity_Seehornwald,bla,13,SWR_Incoming_Seehornwald)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
% xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
% ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Sodankyla
sp3 = subplot(4,2,3);
hold on
text(0.4+0.05*0.8,25-0.05*50,'c','FontSize',12,'FontWeight','bold')
bla = (LWR_Subcanopy_CLM_Sodankyla_avg - LWR_Subcanopy_Observations_Sodankyla_avg)...
    ./LWR_Subcanopy_Observations_Sodankyla_avg*100;
scatter(Sky_Emissivity_Sodankyla,bla,13,SWR_Incoming_Sodankyla)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
% xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Cherskiy
sp4 = subplot(4,2,4);
hold on
text(0.4+0.05*0.8,25-0.05*50,'d','FontSize',12,'FontWeight','bold')
bla = (LWR_Subcanopy_CLM_Cherskiy - LWR_Subcanopy_Observations_Cherskiy)...
    ./LWR_Subcanopy_Observations_Cherskiy*100;
scatter(Sky_Emissivity_Cherskiy,bla,13,SWR_Incoming_Cherskiy)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
% xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
% ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Abisko
sp5 = subplot(4,2,5);
hold on
text(0.4+0.05*0.8,25-0.05*50,'e','FontSize',12,'FontWeight','bold')
bla = (LWR_Subcanopy_CLM_Abisko_avg - LWR_Subcanopy_Observations_Abisko_avg)...
    ./LWR_Subcanopy_Observations_Abisko_avg*100;
scatter(Sky_Emissivity_Abisko,bla,13,SWR_Incoming_Abisko)
hold off
colormap(hot(100))
% cb1=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
% xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
% ylabel(cb1,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Yakutsk
sp6 = subplot(4,2,6);
hold on
text(0.4+0.05*0.8,25-0.05*50,'f','FontSize',12,'FontWeight','bold')
% bla = (LWR_Subcanopy_CLM_Yakutsk - LWR_Subcanopy_Observations_Yakutsk)...
%     ./LWR_Subcanopy_Observations_Yakutsk*100;
bla = (LWR_Subcanopy_CLM_Yakutsk_Night - LWR_Subcanopy_Observations_Yakutsk_Night)...
    ./LWR_Subcanopy_Observations_Yakutsk_Night*100;
scatter(Sky_Emissivity_Yakutsk_Night,bla,13,SWR_Incoming_Yakutsk_Night)
hold off
colormap(hot(100))
% cb2=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
% ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
% ylabel(cb2,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Borden
sp7 = subplot(4,2,7);
hold on
text(0.4+0.05*0.8,25-0.05*50,'g','FontSize',12,'FontWeight','bold')
bla = (LWR_Subcanopy_CLM_Borden - LWR_Subcanopy_Observations_Borden)...
    ./LWR_Subcanopy_Observations_Borden*100;
scatter(Sky_Emissivity_Borden,bla,13,SWR_Incoming_Borden)
hold off
colormap(hot(100))
cb3=colorbar('SouthOutside');
caxis([0 1100])
ylim([-25 25])
xlim([0.4 1.2])
set(gca,'XTick',0.4:0.2:1.2)
xlabel('effective emissivity of the sky','FontSize',12,'FontWeight','bold')
ylabel('sub-canopy LWR error [%]','FontSize',12,'FontWeight','bold')
ylabel(cb3,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% positioning & print
set(sp1,'Position',[0.15,0.76,0.3,0.17])
set(sp2,'Position',[0.55,0.76,0.3,0.17])
set(sp3,'Position',[0.15,0.54,0.3,0.17])
set(sp4,'Position',[0.55,0.54,0.3,0.17])
set(sp5,'Position',[0.15,0.32,0.3,0.17])
set(sp6,'Position',[0.55,0.32,0.3,0.17])
set(sp7,'Position',[0.15,0.1,0.3,0.17])
set(cb3,'Position',[0.55,0.1,0.3,0.03])
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.5 -0.5 11.5 16];
fig.PaperSize = [10 15];
print(fig,'-dpng','-r600','Overview_LWsub_ScatterEmskySWac.png')
set(gcf, 'Renderer', 'opengl')
print(fig,'-dpdf','-r600','Overview_LWsub_ScatterEmskySWac.pdf')


%----------------------  scattered LWE comparisons  ----------------------%
% calculate LWE
LWE_Subcanopy_Observations_Abisko_avg ...
    = LWR_Subcanopy_Observations_Abisko_avg./LWR_Atmosphere_Abisko;
LWE_Subcanopy_Observations_Alptal ...
    = LWR_Subcanopy_Observations_Alptal./LWR_Atmosphere_Alptal;
LWE_Subcanopy_Observations_Borden ...
    = LWR_Subcanopy_Observations_Borden./LWR_Atmosphere_Borden;
LWE_Subcanopy_Observations_Cherskiy ...
    = LWR_Subcanopy_Observations_Cherskiy./LWR_Atmosphere_Cherskiy;
LWE_Subcanopy_Observations_Seehornwald ...
    = LWR_Subcanopy_Observations_Seehornwald./LWR_Atmosphere_Seehornwald;
LWE_Subcanopy_Observations_Sodankyla_avg ...
    = LWR_Subcanopy_Observations_Sodankyla_avg./LWR_Atmosphere_Sodankyla;
LWE_Subcanopy_Observations_Yakutsk ...
    = LWR_Subcanopy_Observations_Yakutsk_Night./LWR_Atmosphere_Yakutsk_Night;
LWE_Subcanopy_CLM_Abisko_avg ...
    = LWR_Subcanopy_CLM_Abisko_avg./LWR_Atmosphere_Abisko;
LWE_Subcanopy_CLM_Alptal ...
    = LWR_Subcanopy_CLM_Alptal./LWR_Atmosphere_Alptal;
LWE_Subcanopy_CLM_Borden ...
    = LWR_Subcanopy_CLM_Borden./LWR_Atmosphere_Borden;
LWE_Subcanopy_CLM_Cherskiy ...
    = LWR_Subcanopy_CLM_Cherskiy./LWR_Atmosphere_Cherskiy;
LWE_Subcanopy_CLM_Seehornwald ...
    = LWR_Subcanopy_CLM_Seehornwald./LWR_Atmosphere_Seehornwald;
LWE_Subcanopy_CLM_Sodankyla_avg ...
    = LWR_Subcanopy_CLM_Sodankyla_avg./LWR_Atmosphere_Sodankyla;
LWE_Subcanopy_CLM_Yakutsk ...
    = LWR_Subcanopy_CLM_Yakutsk_Night./LWR_Atmosphere_Yakutsk_Night;
LWE_Subcanopy_CLM_Abisko_avg_HM ...
    = LWR_Subcanopy_CLM_Abisko_avg_HM./LWR_Atmosphere_Abisko;
LWE_Subcanopy_CLM_Alptal_HM ...
    = LWR_Subcanopy_CLM_Alptal_HM./LWR_Atmosphere_Alptal;
LWE_Subcanopy_CLM_Borden_HM ...
    = LWR_Subcanopy_CLM_Borden_HM./LWR_Atmosphere_Borden;
LWE_Subcanopy_CLM_Cherskiy_HM ...
    = LWR_Subcanopy_CLM_Cherskiy_HM./LWR_Atmosphere_Cherskiy;
LWE_Subcanopy_CLM_Seehornwald_HM ...
    = LWR_Subcanopy_CLM_Seehornwald_HM./LWR_Atmosphere_Seehornwald;
LWE_Subcanopy_CLM_Sodankyla_avg_HM ...
    = LWR_Subcanopy_CLM_Sodankyla_avg_HM./LWR_Atmosphere_Sodankyla;
LWE_Subcanopy_CLM_Yakutsk_HM ...
    = LWR_Subcanopy_CLM_Yakutsk_HM_Night./LWR_Atmosphere_Yakutsk_Night;

% create graph
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
% Alptal
sp1 = subplot(4,2,1);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'a','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Alptal,...
    LWE_Subcanopy_CLM_Alptal,13,SWR_Incoming_Alptal)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
% xlabel('observed LWE','FontSize',12,'FontWeight','bold')
ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Seehornwald
sp2 = subplot(4,2,2);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'b','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Seehornwald,LWE_Subcanopy_CLM_Seehornwald,...
    13,SWR_Incoming_Seehornwald)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
% xlabel('observed LWE','FontSize',12,'FontWeight','bold')
% ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Sodankyla
sp3 = subplot(4,2,3);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'c','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Sodankyla_avg,LWE_Subcanopy_CLM_Sodankyla_avg,...
    13,SWR_Incoming_Sodankyla)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
% xlabel('observed LWE','FontSize',12,'FontWeight','bold')
ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Cherskiy
sp4 = subplot(4,2,4);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'d','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Cherskiy,LWE_Subcanopy_CLM_Cherskiy,...
    13,SWR_Incoming_Cherskiy)
hold off
colormap(hot(100))
% c=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
% xlabel('observed LWE','FontSize',12,'FontWeight','bold')
% ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
% ylabel(c,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Abisko
sp5 = subplot(4,2,5);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'e','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Abisko_avg,LWE_Subcanopy_CLM_Abisko_avg,...
    13,SWR_Incoming_Abisko)
hold off
colormap(hot(100))
% cb1=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
% xlabel('observed LWE','FontSize',12,'FontWeight','bold')
ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
% ylabel(cb1,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Yakutsk
sp6 = subplot(4,2,6);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'f','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Yakutsk,LWE_Subcanopy_CLM_Yakutsk,...
    13,SWR_Incoming_Yakutsk_Night)
hold off
colormap(hot(100))
% cb2=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
xlabel('observed LWE','FontSize',12,'FontWeight','bold')
% ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
% ylabel(cb2,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% Borden
sp7 = subplot(4,2,7);
hold on
plot([0.8 2.2],[0.8 2.2],'k')
text(0.8+0.05*1.4,2.2-0.05*1.4,'g','FontSize',12,'FontWeight','bold')
scatter(LWE_Subcanopy_Observations_Borden,LWE_Subcanopy_CLM_Borden,...
    13,SWR_Incoming_Borden)
hold off
colormap(hot(100))
cb3=colorbar('SouthOutside');
caxis([0 1100])
ylim([0.8 2.2])
xlim([0.8 2.2])
set(gca,'XTick',0.8:0.2:2.2)
set(gca,'YTick',0.8:0.2:2.2)
xlabel('observed LWE','FontSize',12,'FontWeight','bold')
ylabel('simulated LWE','FontSize',12,'FontWeight','bold')
ylabel(cb3,'incoming SWR [W m^{-2}]','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
box on
% positioning & print
set(sp1,'Position',[0.15,0.76,0.3,0.17])
set(sp2,'Position',[0.55,0.76,0.3,0.17])
set(sp3,'Position',[0.15,0.54,0.3,0.17])
set(sp4,'Position',[0.55,0.54,0.3,0.17])
set(sp5,'Position',[0.15,0.32,0.3,0.17])
set(sp6,'Position',[0.55,0.32,0.3,0.17])
set(sp7,'Position',[0.15,0.1,0.3,0.17])
set(cb3,'Position',[0.55,0.1,0.3,0.03])
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.5 -0.5 11.5 16];
fig.PaperSize = [10 15];
print(fig,'-dpng','-r600','Overview_LWE_ScatterEmskySWac.png')
set(gcf, 'Renderer', 'opengl')
print(fig,'-dpdf','-r600','Overview_LWE_ScatterEmskySWac.pdf')


%-----------------------------  PDFs of LWE  -----------------------------%
% calculate PDFs
spectrum_LWenh = 0.8:0.05:2;
spectrum_LWenh_xaxis = 0.825:0.05:1.975;
hist_obs = histogram(LWE_Subcanopy_Observations_Abisko_avg,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_Abi = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Abisko_avg,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_Abi = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Abisko_avg_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_Abi = hist_clm_hm.Values;
hist_obs = histogram(LWE_Subcanopy_Observations_Alptal,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_Alp = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Alptal,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_Alp = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Alptal_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_Alp = hist_clm_hm.Values;
hist_obs = histogram(LWE_Subcanopy_Observations_Borden,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_Bor = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Borden,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_Bor = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Borden_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_Bor = hist_clm_hm.Values;
hist_obs = histogram(LWE_Subcanopy_Observations_Cherskiy,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_Che = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Cherskiy,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_Che = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Cherskiy_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_Che = hist_clm_hm.Values;
hist_obs = histogram(LWE_Subcanopy_Observations_Seehornwald,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_SHW = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Seehornwald,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_SHW = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Seehornwald_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_SHW = hist_clm_hm.Values;
hist_obs = histogram(LWE_Subcanopy_Observations_Sodankyla_avg,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_Sod = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Sodankyla_avg,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_Sod = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Sodankyla_avg_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_Sod = hist_clm_hm.Values;
hist_obs = histogram(LWE_Subcanopy_Observations_Yakutsk,spectrum_LWenh,...
    'Normalization','probability');
    hist_obs_Yak = hist_obs.Values;
hist_clm = histogram(LWE_Subcanopy_CLM_Yakutsk,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_Yak = hist_clm.Values;
hist_clm_hm = histogram(LWE_Subcanopy_CLM_Yakutsk_HM,spectrum_LWenh,...
    'Normalization','probability');
    hist_clm_hm_Yak = hist_clm_hm.Values;


% create graph
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
% Alptal
sp1 = subplot(4,2,1);
hold on
text(0.8+0.05*1.2,0.95*0.3,'a','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Alp,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Alp,'Color',EcstaticEmerald,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_Alp,'Color',EcstaticEmerald,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
% xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% Seehornwald
sp2 = subplot(4,2,2);
hold on
text(0.8+0.05*1.2,0.95*0.15,'b','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_SHW,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_SHW,'Color',MagicMaroon,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_SHW,'Color',MagicMaroon,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
% xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
% ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% Sodankyla
sp3 = subplot(4,2,3);
hold on
text(0.8+0.05*1.2,0.95*0.3,'c','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Sod,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Sod,'Color',GlacialGrey,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_Sod,'Color',GlacialGrey,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
% xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% Cherskiy
sp4 = subplot(4,2,4);
hold on
text(0.8+0.05*1.2,0.95*0.3,'d','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Che,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Che,'Color',BastilleBleu,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_Che,'Color',BastilleBleu,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
% xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
% ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% Abisko
sp5 = subplot(4,2,5);
hold on
text(0.8+0.05*1.2,0.95*0.4,'e','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Abi,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Abi,'Color',GorgeousGold,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_Abi,'Color',GorgeousGold,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
% xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% Yakutsk
sp6 = subplot(4,2,6);
hold on
text(0.8+0.05*1.2,0.95*0.2,'f','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Yak,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Yak,'Color',ViolentViolet,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_Yak,'Color',ViolentViolet,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
% ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% Borden
sp7 = subplot(4,2,7);
hold on
text(0.8+0.05*1.2,0.95*0.4,'g','FontSize',12,'FontWeight','bold')
plot(spectrum_LWenh_xaxis,hist_obs_Bor,'k','LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_Bor,'Color',RussianRed,'LineWidth',2)
plot(spectrum_LWenh_xaxis,hist_clm_hm_Bor,'Color',RussianRed,...
    'LineWidth',2,'LineStyle',':')
hold off
xlim([0.8 2])
xlabel('LW enhancement','FontSize',12,'FontWeight','bold')
ylabel('Frequency','FontSize',12,'FontWeight','bold')
set(gca,'FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.8:0.2:2)
box on
% positioning & print
set(sp1,'Position',[0.15,0.76,0.3,0.17])
set(sp2,'Position',[0.55,0.76,0.3,0.17])
set(sp3,'Position',[0.15,0.54,0.3,0.17])
set(sp4,'Position',[0.55,0.54,0.3,0.17])
set(sp5,'Position',[0.15,0.32,0.3,0.17])
set(sp6,'Position',[0.55,0.32,0.3,0.17])
set(sp7,'Position',[0.15,0.1,0.3,0.17])
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.5 -0.5 11.5 16];
fig.PaperSize = [10 15];
print(fig,'-dpng','-r600','Overview_PDFsLWE.png')
set(gcf, 'Renderer', 'opengl')
print(fig,'-dpdf','-r600','Overview_PDFsLWE.pdf')


%-------------------------------------------------------------------------%
%------------------  ranges of vegetation temperatures  ------------------%
%-------------------------------------------------------------------------%
%{
Instead of showing multiple diurnal cycles, we create a bar plot of
temperature ranges. As measurements of vegetation temperature aren't always
available we infer it from sub-canopy longwave radiation using CLM's own
calculations.
%}
%------------------  infer real vegetation temperature  ------------------%
% Abisko
T_veg_Abi = nan(length(LWR_Atmosphere_Abisko),2);
T_veg_Abi(:,1) = VegTempInfer(PAI_Abisko,LWR_Atmosphere_Abisko,...
    LWR_Subcanopy_Observations_Abisko_avg);
T_veg_Abi(:,2) = VegTempInfer(PAI_Abisko,LWR_Atmosphere_Abisko,...
    LWR_Subcanopy_CLM_Abisko_avg);

% Alptal
T_veg_Alp04 = nan(length(LWR_Atmosphere_Alptal_2004),2);
T_veg_Alp04(:,1) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2004,...
    LWR_Subcanopy_Observations_Alptal_2004);
T_veg_Alp04(:,2) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2004,...
    LWR_Subcanopy_CLM_Alptal_2004);
T_veg_Alp05 = nan(length(LWR_Atmosphere_Alptal_2005),2);
T_veg_Alp05(:,1) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2005,...
    LWR_Subcanopy_Observations_Alptal_2005);
T_veg_Alp05(:,2) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2005,...
    LWR_Subcanopy_CLM_Alptal_2005);
T_veg_Alp06 = nan(length(LWR_Atmosphere_Alptal_2006),2);
T_veg_Alp06(:,1) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2006,...
    LWR_Subcanopy_Observations_Alptal_2006);
T_veg_Alp06(:,2) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2006,...
    LWR_Subcanopy_CLM_Alptal_2006);
T_veg_Alp07 = nan(length(LWR_Atmosphere_Alptal_2007),2);
T_veg_Alp07(:,1) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2007,...
    LWR_Subcanopy_Observations_Alptal_2007);
T_veg_Alp07(:,2) = VegTempInfer(PAI_Alptal,LWR_Atmosphere_Alptal_2007,...
    LWR_Subcanopy_CLM_Alptal_2007);

% Cherskiy
T_veg_Che = nan(length(LWR_Atmosphere_Cherskiy),2);
T_veg_Che(:,1) = VegTempInfer(PAI_Cherskiy,LWR_Atmosphere_Cherskiy,...
    LWR_Subcanopy_Observations_Cherskiy);
T_veg_Che(:,2) = VegTempInfer(PAI_Cherskiy,LWR_Atmosphere_Cherskiy,...
    LWR_Subcanopy_CLM_Cherskiy);

% Seehornwald
T_veg_SHW08 = nan(length(LWR_Atmosphere_Seehornwald_2008),2);
T_veg_SHW08(:,1) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2008,...
    LWR_Subcanopy_Observations_Seehornwald_2008);
T_veg_SHW08(:,2) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2008,...
    LWR_Subcanopy_CLM_Seehornwald_2008);
T_veg_SHW09 = nan(length(LWR_Atmosphere_Seehornwald_2009),2);
T_veg_SHW09(:,1) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2009,...
    LWR_Subcanopy_Observations_Seehornwald_2009);
T_veg_SHW09(:,2) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2009,...
    LWR_Subcanopy_CLM_Seehornwald_2009);
T_veg_SHW10 = nan(length(LWR_Atmosphere_Seehornwald_2010),2);
T_veg_SHW10(:,1) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2010,...
    LWR_Subcanopy_Observations_Seehornwald_2010);
T_veg_SHW10(:,2) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2010,...
    LWR_Subcanopy_CLM_Seehornwald_2010);
T_veg_SHW11 = nan(length(LWR_Atmosphere_Seehornwald_2011),2);
T_veg_SHW11(:,1) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2011,...
    LWR_Subcanopy_Observations_Seehornwald_2011);
T_veg_SHW11(:,2) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2011,...
    LWR_Subcanopy_CLM_Seehornwald_2011);
T_veg_SHW12 = nan(length(LWR_Atmosphere_Seehornwald_2012),2);
T_veg_SHW12(:,1) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2012,...
    LWR_Subcanopy_Observations_Seehornwald_2012);
T_veg_SHW12(:,2) = VegTempInfer(PAI_Seehornwald,LWR_Atmosphere_Seehornwald_2012,...
    LWR_Subcanopy_CLM_Seehornwald_2012);

% Sodankyla
T_veg_Sod = nan(length(LWR_Atmosphere_Sodankyla),2);
T_veg_Sod(:,1) = VegTempInfer(PAI_Sodankyla,LWR_Atmosphere_Sodankyla,...
    LWR_Subcanopy_Observations_Sodankyla_avg);
T_veg_Sod(:,2) = VegTempInfer(PAI_Sodankyla,LWR_Atmosphere_Sodankyla,...
    LWR_Subcanopy_CLM_Sodankyla_avg);


%---------------  calculate vegetation temperature ranges  ---------------%
%{
There basically are two approaches:
1) calculate the average diurnal cycle and afterwards its range; and
2) calculate diurnal minima and maxima for each day and average each of
them afterwards.
%}
% diurnal cycles
T_veg_CLM_Abi_dc = nan(24,1); T_veg_est_Abi_dc = nan(24,1);
T_veg_CLM_Alp04_dc = nan(24,1); T_veg_est_Alp04_dc = nan(24,1);
T_veg_CLM_Alp05_dc = nan(24,1); T_veg_est_Alp05_dc = nan(24,1);
T_veg_CLM_Alp06_dc = nan(24,1); T_veg_est_Alp06_dc = nan(24,1);
T_veg_CLM_Alp07_dc = nan(24,1); T_veg_est_Alp07_dc = nan(24,1);
T_veg_CLM_Che_dc = nan(24,1); T_veg_est_Che_dc = nan(24,1);
T_veg_CLM_SHW08_dc = nan(24,1); T_veg_est_SHW08_dc = nan(24,1);
T_veg_CLM_SHW09_dc = nan(24,1); T_veg_est_SHW09_dc = nan(24,1);
T_veg_CLM_SHW10_dc = nan(24,1); T_veg_est_SHW10_dc = nan(24,1);
T_veg_CLM_SHW11_dc = nan(24,1); T_veg_est_SHW11_dc = nan(24,1);
T_veg_CLM_SHW12_dc = nan(24,1); T_veg_est_SHW12_dc = nan(24,1);
T_veg_CLM_Sod_dc = nan(24,1); T_veg_est_Sod_dc = nan(24,1);
for h=1:24
% Abisko
    T_veg_est_Abi_dc(h) = mean(T_veg_Abi(h:24:end-24+h,1));
    T_veg_CLM_Abi_dc(h) = mean(T_veg_Abi(h:24:end-24+h,2));
% Alptal 2004
    T_veg_est_Alp04_dc(h) = mean(T_veg_Alp04(h:24:end-24+h,1));
    T_veg_CLM_Alp04_dc(h) = mean(T_veg_Alp04(h:24:end-24+h,2));
% Alptal 2005
    T_veg_est_Alp05_dc(h) = mean(T_veg_Alp05(h:24:end-24+h,1));
    T_veg_CLM_Alp05_dc(h) = mean(T_veg_Alp05(h:24:end-24+h,2));
% Alptal 2006
    T_veg_est_Alp06_dc(h) = mean(T_veg_Alp06(h:24:end-24+h,1));
    T_veg_CLM_Alp06_dc(h) = mean(T_veg_Alp06(h:24:end-24+h,2));
% Alptal 2007
    T_veg_est_Alp07_dc(h) = mean(T_veg_Alp07(h:24:end-24+h,1));
    T_veg_CLM_Alp07_dc(h) = mean(T_veg_Alp07(h:24:end-24+h,2));
% Cherskiy
    T_veg_est_Che_dc(h) = mean(T_veg_Che(h:24:end-24+h,1));
    T_veg_CLM_Che_dc(h) = mean(T_veg_Che(h:24:end-24+h,2));
% Seehornwald 2008
    T_veg_est_SHW08_dc(h) = mean(T_veg_SHW08(h:24:end-24+h,1));
    T_veg_CLM_SHW08_dc(h) = mean(T_veg_SHW08(h:24:end-24+h,2));
% Seehornwald 2009
    T_veg_est_SHW09_dc(h) = mean(T_veg_SHW09(h:24:end-24+h,1));
    T_veg_CLM_SHW09_dc(h) = mean(T_veg_SHW09(h:24:end-24+h,2));
% Seehornwald 2010
    T_veg_est_SHW10_dc(h) = mean(T_veg_SHW10(h:24:end-24+h,1));
    T_veg_CLM_SHW10_dc(h) = mean(T_veg_SHW10(h:24:end-24+h,2));
% Seehornwald 2011
    T_veg_est_SHW11_dc(h) = mean(T_veg_SHW11(h:24:end-24+h,1));
    T_veg_CLM_SHW11_dc(h) = mean(T_veg_SHW11(h:24:end-24+h,2));
% Seehornwald 2012
    T_veg_est_SHW12_dc(h) = mean(T_veg_SHW12(h:24:end-24+h,1));
    T_veg_CLM_SHW12_dc(h) = mean(T_veg_SHW12(h:24:end-24+h,2));
% Sodankylä
    T_veg_est_Sod_dc(h) = mean(T_veg_Sod(h:24:end-24+h,1));
    T_veg_CLM_Sod_dc(h) = mean(T_veg_Sod(h:24:end-24+h,2));
end

fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
plot([1 1],[min(T_veg_CLM_Alp04_dc) max(T_veg_CLM_Alp04_dc)]-273.15,...
    'Color',EcstaticEmerald,'LineWidth',23)
plot([2 2],[min(T_veg_est_Alp04_dc) max(T_veg_est_Alp04_dc)]-273.15,...
    'k','LineWidth',23)
plot([4 4],[min(T_veg_CLM_Alp05_dc) max(T_veg_CLM_Alp05_dc)]-273.15,...
    'Color',EcstaticEmerald,'LineWidth',23)
plot([5 5],[min(T_veg_est_Alp05_dc) max(T_veg_est_Alp05_dc)]-273.15,...
    'k','LineWidth',23)
plot([7 7],[min(T_veg_CLM_Alp06_dc) max(T_veg_CLM_Alp06_dc)]-273.15,...
    'Color',EcstaticEmerald,'LineWidth',23)
plot([8 8],[min(T_veg_est_Alp06_dc) max(T_veg_est_Alp06_dc)]-273.15,...
    'k','LineWidth',23)
plot([10 10],[min(T_veg_CLM_Alp07_dc) max(T_veg_CLM_Alp07_dc)]-273.15,...
    'Color',EcstaticEmerald,'LineWidth',23)
plot([11 11],[min(T_veg_est_Alp07_dc) max(T_veg_est_Alp07_dc)]-273.15,...
    'k','LineWidth',23)
plot([15 15],[min(T_veg_CLM_SHW08_dc) max(T_veg_CLM_SHW08_dc)]-273.15,...
    'Color',MagicMaroon,'LineWidth',23)
plot([16 16],[min(T_veg_est_SHW08_dc) max(T_veg_est_SHW08_dc)]-273.15,...
    'k','LineWidth',23)
plot([18 18],[min(T_veg_CLM_SHW09_dc) max(T_veg_CLM_SHW09_dc)]-273.15,...
    'Color',MagicMaroon,'LineWidth',23)
plot([19 19],[min(T_veg_est_SHW09_dc) max(T_veg_est_SHW09_dc)]-273.15,...
    'k','LineWidth',23)
plot([21 21],[min(T_veg_CLM_SHW10_dc) max(T_veg_CLM_SHW10_dc)]-273.15,...
    'Color',MagicMaroon,'LineWidth',23)
plot([22 22],[min(T_veg_est_SHW10_dc) max(T_veg_est_SHW10_dc)]-273.15,...
    'k','LineWidth',23)
plot([24 24],[min(T_veg_CLM_SHW11_dc) max(T_veg_CLM_SHW11_dc)]-273.15,...
    'Color',MagicMaroon,'LineWidth',23)
plot([25 25],[min(T_veg_est_SHW11_dc) max(T_veg_est_SHW11_dc)]-273.15,...
    'k','LineWidth',23)
plot([27 27],[min(T_veg_CLM_SHW12_dc) max(T_veg_CLM_SHW12_dc)]-273.15,...
    'Color',MagicMaroon,'LineWidth',23)
plot([28 28],[min(T_veg_est_SHW12_dc) max(T_veg_est_SHW12_dc)]-273.15,...
    'k','LineWidth',23)
plot([32 32],[min(T_veg_CLM_Sod_dc) max(T_veg_CLM_Sod_dc)]-273.15,...
    'Color',GlacialGrey,'LineWidth',23)
plot([33 33],[min(T_veg_est_Sod_dc) max(T_veg_est_Sod_dc)]-273.15,...
    'k','LineWidth',23)
plot([37 37],[min(T_veg_CLM_Che_dc) max(T_veg_CLM_Che_dc)]-273.15,...
    'Color',BastilleBleu,'LineWidth',23)
plot([38 38],[min(T_veg_est_Che_dc) max(T_veg_est_Che_dc)]-273.15,...
    'k','LineWidth',23)
plot([42 42],[min(T_veg_CLM_Abi_dc) max(T_veg_CLM_Abi_dc)]-273.15,...
    'Color',GorgeousGold,'LineWidth',23)
plot([43 43],[min(T_veg_est_Abi_dc) max(T_veg_est_Abi_dc)]-273.15,...
    'k','LineWidth',23)
xlim([0 44])
box on
ylabel('vegetation temperature [°C]')
set(gca,'XTick',[6 21.5 32.5 37.5 42.5],'XTickLabel',...
    {'Alptal','Seehornwald','Sodankylä','Cherskiy','Abisko'})
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
rotateXLabels(gca,30)
fig.PaperUnits = 'inches';
fig.PaperPosition = [-0.5 0 10.5 6];
fig.PaperSize = [9.5 5.9];
print(fig,'-dpng','-r600','SiteComp_VegTemp.png')
print(fig,'-dpdf','-r600','SiteComp_VegTemp.pdf')



toc