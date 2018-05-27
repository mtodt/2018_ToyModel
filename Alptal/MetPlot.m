function MetPlot(SiteName,SiteLabel,Season,time,timeticks,datelabel,...
    LW_atm,T_air,SW_in)
%{
Required input to create meteorological overvier figure used for initial
Toy Model paper. Met graph was outsourced due to its complexity and since
it doesn't change at a site if Toy Model parameters get varied.
General & Time:
    - SiteName  -> site for file name (string)
    - SiteLabel -> site for x-axis label (string)
    - Season    -> year (string)
    - time      -> time in days since 1 January 0 (array)
    - timeticks -> ticks on time axis for labels (array)
    - datelabel -> desired labels on time axis (string array)
For example, used for initial Alptal figure:
timeticks = [time(4009) time(4418) time(4765) ... ];
datelabel = {'15 Jan' '01 Feb' '15 Feb' ...};
Meteorological forcing:
    - LW_atm    -> atmospheric longwave radiation forcing (array)
    - T_air     -> air temperature forcing (array)
    - SW_in     -> incoming shortwave radiation forcing (array)
%}
    
%---------------------------  Plot Essentials  ---------------------------%
close all
if exist('fig_num','var') == 0
    fig_num = 1;
end
% colours
GargantuanGreen = [0.1 0.5 0.2];
matlabblue = [0 0.447 0.741];
matlabred = [0.85 0.325 0.098];
matlabyellow = [0.929 0.694 0.125];

%---------------------------  Post-Processing  ---------------------------%
% effective emissivity of the sky
em_sky = nan(size(T_air));
for t=1:length(time)
    em_sky(t) = LW_atm(t)/(5.67*10^(-8)*T_air(t)^4);
end
T_air = T_air-273.15;   % convert from K for figures
% plot axes' limits
    % met
SWmax = max(SW_in);
SWmax = SWmax/100; SWmax = ceil(SWmax); SWmax = SWmax*100;
if mod(SWmax/100,2) == 1
    SWmax = SWmax+100;
end
LWmax = max(LW_atm);
LWmax = LWmax/50; LWmax = ceil(LWmax); LWmax = LWmax*50;
LWmin = min(LW_atm);
LWmin = LWmin/50; LWmin = floor(LWmin); LWmin = LWmin*50;
Tairmax = max(T_air);
Tairmax = Tairmax/5; Tairmax = ceil(Tairmax); Tairmax = Tairmax*5;
Tairmin = min(T_air);
Tairmin = Tairmin/5; Tairmin = floor(Tairmin); Tairmin = Tairmin*5;

%--------------------------------  Graph  --------------------------------%
x = time;
y1 = LW_atm;
y2 = T_air;
y3 = SW_in;
y4 = nan(size(time));
y5 = nan(size(time));
y6 = em_sky;
fig=figure(fig_num);fig_num = fig_num+1;
set(gcf,'Position',get(0,'ScreenSize'))
hold on
[ax,h1,h2] = plotyy(x,y1,x,y2);
ax(3) = axes('yaxislocation','left','Color','none','XColor','k','YColor',matlabyellow,'box','off');
h3 = line(x,y3,'Parent',ax(3),'Color',matlabyellow);
ax(4) = axes('yaxislocation','right','Color','none','XColor','k','YColor',GargantuanGreen,'box','off');
h4 = line(x,y6,'Parent',ax(4),'Color',GargantuanGreen);
set(ax(1),'XTick',timeticks,'XTickLabel','');%set(ax(1),'XTick',[4418 5114 5858 6578]); datetick(ax(1),'x','dd mmm')
set(ax(2),'XTick',timeticks,'XTickLabel','');%set(ax(2),'XTick',[4418 5114 5858 6578]); datetick(ax(2),'x','dd mmm')
set(ax(3),'XTick',timeticks,'XTickLabel','');%set(ax(3),'XTick',[4418 5114 5858 6578]); datetick(ax(3),'x','dd mmm')
set(ax(4),'XTick',timeticks,'XTickLabel',datelabel);%set(ax(4),'XTick',[4418 5114 5858 6578]); datetick(ax(4),'x','dd mmm')
set(ax(1),'Box','off')
set(ax(2),'Box','off')
set(ax(3),'Box','off')
set(ax(4),'Box','off')
xlabel([SiteLabel ' observation period ' Season],'FontSize',17,'FontWeight','bold')
ylabel(ax(1),{'atmospheric longwave','radiation [W m^{-2}]'},'FontSize',17,'FontWeight','bold')
ylabel(ax(2),{'above-canopy','air temperature [°C]'},'FontSize',17,'FontWeight','bold')
ylabel(ax(3),{'incoming shortwave','radiation [W m^{-2}]'},'FontSize',17,'FontWeight','bold')
ylabel(ax(4),{'effective emissivity','of the sky'},'FontSize',17,'FontWeight','bold')
set(ax(1),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(2),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(3),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(4),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
xlim(ax(1),[time(1) time(end)]); ylim(ax(1),[LWmin-325 LWmax+75]);
xlim(ax(2),[time(1) time(end)]); ylim(ax(2),[Tairmin-(Tairmax-Tairmin) Tairmax+1.1*(Tairmax-Tairmin)]);
xlim(ax(3),[time(1) time(end)]); ylim(ax(3),[0 2.5*SWmax]);
xlim(ax(4),[time(1) time(end)]); ylim(ax(4),[-3 1.2]);
set(ax(1),'YTick',LWmin+50:50:LWmax-50)
if Tairmax-Tairmin <= 35
    set(ax(2),'YTick',Tairmin+5:5:Tairmax-5)
elseif Tairmax-Tairmin > 35
    if mod(Tairmax-Tairmin,10) == 0
        set(ax(2),'YTick',Tairmin+5:10:Tairmax-5)
    else
        set(ax(2),'YTick',Tairmin+10:10:Tairmax-5)
    end
end
set(ax(3),'YTick',0:200:SWmax-200)
set(ax(4),'YTick',0.6:0.2:1.2)
ylab1 = get(ax(1),'YLabel');% set(ylab1,'Position',get(ylab1,'Position') - [6 -(LWmax-LWmin+400)/5 0])
pos1 = get(ylab1,'Position');
set(ylab1,'Position',[pos1(1)-0.01*(time(end)-time(1)) LWmin+0.5*(LWmax-LWmin) pos1(3)])
ylab2 = get(ax(2),'YLabel');% set(ylab2,'Position',get(ylab2,'Position') + [4 -0.05*(Tairmax-Tairmin) 0])
pos2 = get(ylab2,'Position');
set(ylab2,'Position',[pos2(1)+0.01*(time(end)-time(1)) Tairmin+0.5*(Tairmax-Tairmin) pos2(3)])
ylab3 = get(ax(3),'YLabel');% set(ylab3,'Position',get(ylab3,'Position') - [6 SWmax/1.25 0])
pos3 = get(ylab3,'Position');
set(ylab3,'Position',[pos3(1)-0.01*(time(end)-time(1)) 0.5*SWmax pos3(3)])
ylab4 = get(ax(4),'YLabel');% set(ylab4,'Position',get(ylab4,'Position') + [4 1.78 0])
pos4 = get(ylab4,'Position');
set(ylab4,'Position',[pos4(1)+0.01*(time(end)-time(1)) 0.9 pos4(3)])
ax(5) = axes('yaxislocation','left','Color','none','XColor','k','YColor','k');
ax(6) = axes('yaxislocation','right','Color','none','XColor','k','YColor','k');
h5 = line(x,y4,'Parent',ax(5),'Color','k');
h6 = line(x,y5,'Parent',ax(6),'Color','k');
set(ax(5),'Box','off'); set(ax(6),'Box','off');
set(ax(5),'XTick',timeticks,'XTickLabel','');%set(ax(5),'XTick',[4418 5114 5858 6578]); datetick(ax(5),'x','dd mmm');
set(ax(6),'XTick',timeticks,'XTickLabel','');%set(ax(6),'XTick',[4418 5114 5858 6578]); datetick(ax(6),'x','dd mmm');
set(ax(5),'FontSize',17,'FontWeight','bold','LineWidth',2)
set(ax(6),'FontSize',17,'FontWeight','bold','LineWidth',2)
xlim(ax(5),[time(1) time(end)]); ylim(ax(5),[0 1]);
xlim(ax(6),[time(1) time(end)]); ylim(ax(6),[0 1]);
set(ax(5),'YTick',[]); set(ax(6),'YTick',[]);
print(fig,'-dpng','-r600',[SiteName '_' Season '_MetForcEmsky.png'])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 20 10];
fig.PaperSize = [20 11];
print(fig,'-dpdf','-r600',[SiteName '_' Season '_MetForcEmsky.pdf'])
close
end