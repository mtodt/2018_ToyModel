tic

clear
close all
%{
Use MLR input to create met plots - evaluation period already cut out.
%}

%-----------------------------  import data  -----------------------------%
load MLRdata_Abisko.mat
load MLRdata_Alptal.mat
load MLRdata_Borden.mat
load MLRdata_Cherskiy.mat
load MLRdata_Seehornwald.mat
load MLRdata_Sodankyla.mat
load MLRdata_Yakutsk.mat
SWmax_Yak = 767.9;  % need value for figure

%----------------------  calculate air temperature  ----------------------%
boltz = 5.67*10^(-8);
    % Abisko
T_air_Abi = nan(size(LWR_Atmosphere_Abisko));
for l=1:length(LWR_Atmosphere_Abisko)
    T_air_Abi(l) = nthroot(LWR_Atmosphere_Abisko(l)/...
        (Sky_Emissivity_Abisko(l)*boltz),4)-273.15;
end
    % Alptal
T_air_Alp04 = nan(size(LWR_Atmosphere_Alptal_2004));
for l=1:length(LWR_Atmosphere_Alptal_2004)
    T_air_Alp04(l) = nthroot(LWR_Atmosphere_Alptal_2004(l)/...
        (Sky_Emissivity_Alptal_2004(l)*boltz),4)-273.15;
end
T_air_Alp05 = nan(size(LWR_Atmosphere_Alptal_2005));
for l=1:length(LWR_Atmosphere_Alptal_2005)
    T_air_Alp05(l) = nthroot(LWR_Atmosphere_Alptal_2005(l)/...
        (Sky_Emissivity_Alptal_2005(l)*boltz),4)-273.15;
end
T_air_Alp06 = nan(size(LWR_Atmosphere_Alptal_2006));
for l=1:length(LWR_Atmosphere_Alptal_2006)
    T_air_Alp06(l) = nthroot(LWR_Atmosphere_Alptal_2006(l)/...
        (Sky_Emissivity_Alptal_2006(l)*boltz),4)-273.15;
end
T_air_Alp07 = nan(size(LWR_Atmosphere_Alptal_2007));
for l=1:length(LWR_Atmosphere_Alptal_2007)
    T_air_Alp07(l) = nthroot(LWR_Atmosphere_Alptal_2007(l)/...
        (Sky_Emissivity_Alptal_2007(l)*boltz),4)-273.15;
end
    % Borden
T_air_Bor = nan(size(LWR_Atmosphere_Borden));
for l=1:length(LWR_Atmosphere_Borden)
    T_air_Bor(l) = nthroot(LWR_Atmosphere_Borden(l)/...
        (Sky_Emissivity_Borden(l)*boltz),4)-273.15;
end
    % Cherskiy
T_air_Che = nan(size(LWR_Atmosphere_Cherskiy));
for l=1:length(LWR_Atmosphere_Cherskiy)
    T_air_Che(l) = nthroot(LWR_Atmosphere_Cherskiy(l)/...
        (Sky_Emissivity_Cherskiy(l)*boltz),4)-273.15;
end
    % Seehornwald
T_air_SHW08 = nan(size(LWR_Atmosphere_Seehornwald_2008));
for l=1:length(LWR_Atmosphere_Seehornwald_2008)
    T_air_SHW08(l) = nthroot(LWR_Atmosphere_Seehornwald_2008(l)/...
        (Sky_Emissivity_Seehornwald_2008(l)*boltz),4)-273.15;
end
T_air_SHW09 = nan(size(LWR_Atmosphere_Seehornwald_2009));
for l=1:length(LWR_Atmosphere_Seehornwald_2009)
    T_air_SHW09(l) = nthroot(LWR_Atmosphere_Seehornwald_2009(l)/...
        (Sky_Emissivity_Seehornwald_2009(l)*boltz),4)-273.15;
end
T_air_SHW10 = nan(size(LWR_Atmosphere_Seehornwald_2010));
for l=1:length(LWR_Atmosphere_Seehornwald_2010)
    T_air_SHW10(l) = nthroot(LWR_Atmosphere_Seehornwald_2010(l)/...
        (Sky_Emissivity_Seehornwald_2010(l)*boltz),4)-273.15;
end
T_air_SHW11 = nan(size(LWR_Atmosphere_Seehornwald_2011));
for l=1:length(LWR_Atmosphere_Seehornwald_2011)
    T_air_SHW11(l) = nthroot(LWR_Atmosphere_Seehornwald_2011(l)/...
        (Sky_Emissivity_Seehornwald_2011(l)*boltz),4)-273.15;
end
T_air_SHW12 = nan(size(LWR_Atmosphere_Seehornwald_2012));
for l=1:length(LWR_Atmosphere_Seehornwald_2012)
    T_air_SHW12(l) = nthroot(LWR_Atmosphere_Seehornwald_2012(l)/...
        (Sky_Emissivity_Seehornwald_2012(l)*boltz),4)-273.15;
end
    % Sodankyla
T_air_Sod = nan(size(LWR_Atmosphere_Sodankyla));
for l=1:length(LWR_Atmosphere_Sodankyla)
    T_air_Sod(l) = nthroot(LWR_Atmosphere_Sodankyla(l)/...
        (Sky_Emissivity_Sodankyla(l)*boltz),4)-273.15;
end
    % Yakutsk
T_air_Yak = nan(size(LWR_Atmosphere_Yakutsk_Night));
for l=1:length(LWR_Atmosphere_Yakutsk_Night)
    T_air_Yak(l) = nthroot(LWR_Atmosphere_Yakutsk_Night(l)/...
        (Sky_Emissivity_Yakutsk_Night(l)*boltz),4)-273.15;
end

%---------------------------  plot essentials  ---------------------------%
if exist('fig_num','var') == 0
    fig_num = 1;
end

screen_size = get(0, 'ScreenSize');

GlacialGrey = [0.4 0.7 1];
MagicMaroon = [0.65 0.32 0.35];
EcstaticEmerald = [0.17 0.52 0.5];
GorgeousGold = [1 0.85 0];
BastilleBleu = [0.12 0.34 0.56];
RussianRed = [1 0.032 0.064];
ViolentViolet = [0.6 0 0.7];

%---------------------------  produce figures  ---------------------------%
% temperature vs insolation
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
% range line
plot([min(T_air_Abi) max(T_air_Abi)],[max(SWR_Incoming_Abisko) ...
    max(SWR_Incoming_Abisko)],'Color',GorgeousGold,'LineWidth',2)
plot([min(T_air_Alp04) max(T_air_Alp04)],[max(SWR_Incoming_Alptal_2004)...
    max(SWR_Incoming_Alptal_2004)],'Color',EcstaticEmerald,'LineWidth',2)
plot([min(T_air_Bor) max(T_air_Bor)],[max(SWR_Incoming_Borden) ...
    max(SWR_Incoming_Borden)],'Color',RussianRed,'LineWidth',2)
plot([min(T_air_Che) max(T_air_Che)],[max(SWR_Incoming_Cherskiy) ...
    max(SWR_Incoming_Cherskiy)],'Color',BastilleBleu,'LineWidth',2)
plot([min(T_air_SHW08) max(T_air_SHW08)],[max(SWR_Incoming_Seehornwald_2008) ...
    max(SWR_Incoming_Seehornwald_2008)],'Color',MagicMaroon,'LineWidth',2)
plot([min(T_air_Sod) max(T_air_Sod)],[max(SWR_Incoming_Sodankyla) ...
    max(SWR_Incoming_Sodankyla)],'Color',GlacialGrey,'LineWidth',2)
plot([min(T_air_Yak) max(T_air_Yak)],[SWmax_Yak SWmax_Yak],...
    'Color',ViolentViolet,'LineWidth',2)
plot([min(T_air_Alp05) max(T_air_Alp05)],[max(SWR_Incoming_Alptal_2005)...
    max(SWR_Incoming_Alptal_2005)],'Color',EcstaticEmerald,'LineWidth',2)
plot([min(T_air_Alp06) max(T_air_Alp06)],[max(SWR_Incoming_Alptal_2006)...
    max(SWR_Incoming_Alptal_2006)],'Color',EcstaticEmerald,'LineWidth',2)
plot([min(T_air_Alp07) max(T_air_Alp07)],[max(SWR_Incoming_Alptal_2007)...
    max(SWR_Incoming_Alptal_2007)],'Color',EcstaticEmerald,'LineWidth',2)
plot([min(T_air_SHW09) max(T_air_SHW09)],[max(SWR_Incoming_Seehornwald_2009) ...
    max(SWR_Incoming_Seehornwald_2009)],'Color',MagicMaroon,'LineWidth',2)
plot([min(T_air_SHW10) max(T_air_SHW10)],[max(SWR_Incoming_Seehornwald_2010) ...
    max(SWR_Incoming_Seehornwald_2010)],'Color',MagicMaroon,'LineWidth',2)
plot([min(T_air_SHW11) max(T_air_SHW11)],[max(SWR_Incoming_Seehornwald_2011) ...
    max(SWR_Incoming_Seehornwald_2011)],'Color',MagicMaroon,'LineWidth',2)
plot([min(T_air_SHW12) max(T_air_SHW12)],[max(SWR_Incoming_Seehornwald_2012) ...
    max(SWR_Incoming_Seehornwald_2012)],'Color',MagicMaroon,'LineWidth',2)
% mean dot
plot(mean(T_air_Abi),max(SWR_Incoming_Abisko),'MarkerEdgeColor',GorgeousGold,...
    'MarkerFaceColor',GorgeousGold,'Marker','o','MarkerSize',13)
plot(mean(T_air_Alp04),max(SWR_Incoming_Alptal_2004),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','o','MarkerSize',13)
plot(mean(T_air_Alp05),max(SWR_Incoming_Alptal_2005),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','o','MarkerSize',13)
plot(mean(T_air_Alp06),max(SWR_Incoming_Alptal_2006),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','o','MarkerSize',13)
plot(mean(T_air_Alp07),max(SWR_Incoming_Alptal_2007),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','o','MarkerSize',13)
plot(mean(T_air_Bor),max(SWR_Incoming_Borden),'MarkerEdgeColor',RussianRed,...
    'MarkerFaceColor',RussianRed,'Marker','o','MarkerSize',13)
plot(mean(T_air_Che),max(SWR_Incoming_Cherskiy),'MarkerEdgeColor',BastilleBleu,...
    'MarkerFaceColor',BastilleBleu,'Marker','o','MarkerSize',13)
plot(mean(T_air_SHW08),max(SWR_Incoming_Seehornwald_2008),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','o','MarkerSize',13)
plot(mean(T_air_SHW09),max(SWR_Incoming_Seehornwald_2009),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','o','MarkerSize',13)
plot(mean(T_air_SHW10),max(SWR_Incoming_Seehornwald_2010),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','o','MarkerSize',13)
plot(mean(T_air_SHW11),max(SWR_Incoming_Seehornwald_2011),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','o','MarkerSize',13)
plot(mean(T_air_SHW12),max(SWR_Incoming_Seehornwald_2012),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','o','MarkerSize',13)
plot(mean(T_air_Sod),max(SWR_Incoming_Sodankyla),'MarkerEdgeColor',GlacialGrey,...
    'MarkerFaceColor',GlacialGrey,'Marker','o','MarkerSize',13)
plot(mean(T_air_Yak),SWmax_Yak,'MarkerEdgeColor',ViolentViolet,...
    'MarkerFaceColor',ViolentViolet,'Marker','o','MarkerSize',13)
% min & max triangles
plot(min(T_air_Abi),max(SWR_Incoming_Abisko),'MarkerEdgeColor',GorgeousGold,...
    'MarkerFaceColor',GorgeousGold,'Marker','>','MarkerSize',13)
plot(min(T_air_Alp04),max(SWR_Incoming_Alptal_2004),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','>','MarkerSize',13)
plot(min(T_air_Alp05),max(SWR_Incoming_Alptal_2005),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','>','MarkerSize',13)
plot(min(T_air_Alp06),max(SWR_Incoming_Alptal_2006),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','>','MarkerSize',13)
plot(min(T_air_Alp07),max(SWR_Incoming_Alptal_2007),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','>','MarkerSize',13)
plot(min(T_air_Bor),max(SWR_Incoming_Borden),'MarkerEdgeColor',RussianRed,...
    'MarkerFaceColor',RussianRed,'Marker','>','MarkerSize',13)
plot(min(T_air_Che),max(SWR_Incoming_Cherskiy),'MarkerEdgeColor',BastilleBleu,...
    'MarkerFaceColor',BastilleBleu,'Marker','>','MarkerSize',13)
plot(min(T_air_SHW08),max(SWR_Incoming_Seehornwald_2008),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','>','MarkerSize',13)
plot(min(T_air_SHW09),max(SWR_Incoming_Seehornwald_2009),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','>','MarkerSize',13)
plot(min(T_air_SHW10),max(SWR_Incoming_Seehornwald_2010),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','>','MarkerSize',13)
plot(min(T_air_SHW11),max(SWR_Incoming_Seehornwald_2011),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','>','MarkerSize',13)
plot(min(T_air_SHW12),max(SWR_Incoming_Seehornwald_2012),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','>','MarkerSize',13)
plot(min(T_air_Sod),max(SWR_Incoming_Sodankyla),'MarkerEdgeColor',GlacialGrey,...
    'MarkerFaceColor',GlacialGrey,'Marker','>','MarkerSize',13)
plot(min(T_air_Yak),SWmax_Yak,'MarkerEdgeColor',ViolentViolet,...
    'MarkerFaceColor',ViolentViolet,'Marker','>','MarkerSize',13)
plot(max(T_air_Abi),max(SWR_Incoming_Abisko),'MarkerEdgeColor',GorgeousGold,...
    'MarkerFaceColor',GorgeousGold,'Marker','<','MarkerSize',13)
plot(max(T_air_Alp04),max(SWR_Incoming_Alptal_2004),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','<','MarkerSize',13)
plot(max(T_air_Alp05),max(SWR_Incoming_Alptal_2005),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','<','MarkerSize',13)
plot(max(T_air_Alp06),max(SWR_Incoming_Alptal_2006),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','<','MarkerSize',13)
plot(max(T_air_Alp07),max(SWR_Incoming_Alptal_2007),'MarkerEdgeColor',EcstaticEmerald,...
    'MarkerFaceColor',EcstaticEmerald,'Marker','<','MarkerSize',13)
plot(max(T_air_Bor),max(SWR_Incoming_Borden),'MarkerEdgeColor',RussianRed,...
    'MarkerFaceColor',RussianRed,'Marker','<','MarkerSize',13)
plot(max(T_air_Che),max(SWR_Incoming_Cherskiy),'MarkerEdgeColor',BastilleBleu,...
    'MarkerFaceColor',BastilleBleu,'Marker','<','MarkerSize',13)
plot(max(T_air_SHW08),max(SWR_Incoming_Seehornwald_2008),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','<','MarkerSize',13)
plot(max(T_air_SHW09),max(SWR_Incoming_Seehornwald_2009),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','<','MarkerSize',13)
plot(max(T_air_SHW10),max(SWR_Incoming_Seehornwald_2010),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','<','MarkerSize',13)
plot(max(T_air_SHW11),max(SWR_Incoming_Seehornwald_2011),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','<','MarkerSize',13)
plot(max(T_air_SHW12),max(SWR_Incoming_Seehornwald_2012),'MarkerEdgeColor',MagicMaroon,...
    'MarkerFaceColor',MagicMaroon,'Marker','<','MarkerSize',13)
plot(max(T_air_Sod),max(SWR_Incoming_Sodankyla),'MarkerEdgeColor',GlacialGrey,...
    'MarkerFaceColor',GlacialGrey,'Marker','<','MarkerSize',13)
plot(max(T_air_Yak),SWmax_Yak,'MarkerEdgeColor',ViolentViolet,...
    'MarkerFaceColor',ViolentViolet,'Marker','<','MarkerSize',13)
hold off
xlabel('Temperature [\circC]','FontSize',17,'FontWeight','bold')
ylabel('Insolation [W m^{-2}]','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
legend('location','northoutside','orientation','horizontal',...
    'Abisko','Alptal','Borden','Cherskiy','Seehornwald','Sodankylä','Yakutsk')
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 13 7];
fig.PaperSize = [12.7 6.6];
print(fig,'-dpng','-r600','TempSW_Paper.png')
print(fig,'-dpdf','-r600','TempSW_Paper.pdf')
close

% PDF effective emissivity of the sky
spectrum_emeff = 0.4:0.025:1.2;
spectrum_emeff_xaxis = 0.4125:0.025:1.1875;

hist = histogram(Sky_Emissivity_Abisko,spectrum_emeff,'Normalization','probability');
    hist_Abi = hist.Values;
hist = histogram(Sky_Emissivity_Alptal_2004,spectrum_emeff,'Normalization','probability');
    hist_Alp04 = hist.Values;
hist = histogram(Sky_Emissivity_Alptal_2005,spectrum_emeff,'Normalization','probability');
    hist_Alp05 = hist.Values;
hist = histogram(Sky_Emissivity_Alptal_2006,spectrum_emeff,'Normalization','probability');
    hist_Alp06 = hist.Values;
hist = histogram(Sky_Emissivity_Alptal_2007,spectrum_emeff,'Normalization','probability');
    hist_Alp07 = hist.Values;
hist = histogram(vertcat(Sky_Emissivity_Alptal_2004,Sky_Emissivity_Alptal_2005,...
    Sky_Emissivity_Alptal_2006,Sky_Emissivity_Alptal_2007),spectrum_emeff,...
    'Normalization','probability');
    hist_Alp = hist.Values;  
hist = histogram(Sky_Emissivity_Borden,spectrum_emeff,'Normalization','probability');
    hist_Bor = hist.Values;
hist = histogram(Sky_Emissivity_Cherskiy,spectrum_emeff,'Normalization','probability');
    hist_Che = hist.Values;
hist = histogram(Sky_Emissivity_Seehornwald_2008,spectrum_emeff,'Normalization','probability');
    hist_SHW08 = hist.Values;
hist = histogram(Sky_Emissivity_Seehornwald_2009,spectrum_emeff,'Normalization','probability');
    hist_SHW09 = hist.Values;
hist = histogram(Sky_Emissivity_Seehornwald_2010,spectrum_emeff,'Normalization','probability');
    hist_SHW10 = hist.Values;
hist = histogram(Sky_Emissivity_Seehornwald_2011,spectrum_emeff,'Normalization','probability');
    hist_SHW11 = hist.Values;
hist = histogram(Sky_Emissivity_Seehornwald_2012,spectrum_emeff,'Normalization','probability');
    hist_SHW12 = hist.Values;
hist = histogram(vertcat(Sky_Emissivity_Seehornwald_2008,Sky_Emissivity_Seehornwald_2009,...
    Sky_Emissivity_Seehornwald_2010,Sky_Emissivity_Seehornwald_2011,...
    Sky_Emissivity_Seehornwald_2012),spectrum_emeff,'Normalization','probability');
    hist_SHW = hist.Values;
hist = histogram(Sky_Emissivity_Sodankyla,spectrum_emeff,'Normalization','probability');
    hist_Sod = hist.Values;
hist = histogram(Sky_Emissivity_Yakutsk_Night,spectrum_emeff,'Normalization','probability');
    hist_Yak = hist.Values;
close

% one line per site, individual seasons pale in background
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
plot(spectrum_emeff_xaxis,hist_Alp04,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp05,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp06,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp07,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW08,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW09,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW10,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW11,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW12,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Abi,'Color',GorgeousGold,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Alp,'Color',EcstaticEmerald,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Bor,'Color',RussianRed,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Che,'Color',BastilleBleu,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_SHW,'Color',MagicMaroon,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Sod,'Color',GlacialGrey,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Yak,'Color',ViolentViolet,'LineWidth',4)
hold off
xlim([0.4 1.2])
%ylim([0 0.17])
xlabel('effective emissivity of the sky','FontSize',17,'FontWeight','bold')
ylabel('Frequency','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.4:0.1:1.2)
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 13 7];
fig.PaperSize = [12.7 6.6];
print(fig,'-dpng','-r600','EmEffPDF_Paper.png')
print(fig,'-dpdf','-r600','EmEffPDF_Paper.pdf')
    % as subplots
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
subplot(2,1,1)
hold on
plot(spectrum_emeff_xaxis,hist_Alp,'Color',EcstaticEmerald,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_SHW,'Color',MagicMaroon,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Sod,'Color',GlacialGrey,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Alp04,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp05,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp06,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp07,'Color',EcstaticEmerald,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW08,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW09,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW10,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW11,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_SHW12,'Color',MagicMaroon,'LineWidth',0.25)
plot(spectrum_emeff_xaxis,hist_Alp,'Color',EcstaticEmerald,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_SHW,'Color',MagicMaroon,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Sod,'Color',GlacialGrey,'LineWidth',4)
hold off
xlim([0.4 1.2])
ylim([0 0.2])
xlabel('effective emissivity of the sky','FontSize',17,'FontWeight','bold')
ylabel('Frequency','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.4:0.1:1.2)
legend('location','northoutside','orientation','horizontal',...
    'Alptal','Seehornwald','Sodankylä')
box on
subplot(2,1,2)
hold on
plot(spectrum_emeff_xaxis,hist_Abi,'Color',GorgeousGold,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Bor,'Color',RussianRed,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Che,'Color',BastilleBleu,'LineWidth',4)
plot(spectrum_emeff_xaxis,hist_Yak,'Color',ViolentViolet,'LineWidth',4)
hold off
xlim([0.4 1.2])
ylim([0 0.2])
xlabel('effective emissivity of the sky','FontSize',17,'FontWeight','bold')
ylabel('Frequency','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0.4:0.1:1.2)
legend('location','northoutside','orientation','horizontal',...
    'Abisko','Borden','Cherskiy','Yakutsk')
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 13 14];
fig.PaperSize = [12.7 13.2];
print(fig,'-dpng','-r600','EmEffPDF_Paper_split.png')
print(fig,'-dpdf','-r600','EmEffPDF_Paper_split.pdf')