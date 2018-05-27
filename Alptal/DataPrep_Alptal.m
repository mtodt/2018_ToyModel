tic
clear
close all

%-------------------------------------------------------------------------%
%-----------------------------  import data  -----------------------------%
%-------------------------------------------------------------------------%

%-----------------------  open site met variables  -----------------------%
%{
Hourly resolution for period 1 October 2003 0:00 to 7 October 2007 23:00.
%}
bla = xlsread('MeteoData_Alptal_Open_Site.xls');
T_air_open_1h = bla(:,1);   % air temperature at 3.5m(?) height [°C]
    T_air_open_1h = T_air_open_1h + 273.15;
RH_open_1h = bla(:,2);      % relative humidity at 3.5m(?) height [%]
Precip_1h = bla(:,4);       % precipitation in [mm]
LW_in_open_1h = bla(:,6);      % LWR to fill gaps
SW_in_open_1h = bla(:,5);
%{
Unfortunately, there are no air temperature measurements above the canopy
beyond SnowMIP2, so we have to use open site measurements.
Precipitation measurements feature a gap from 21 June 2004 0:00 to 27 June
2004 23:00. Fortunately, we will exclude that period when cutting out
individual season, sincec we can't run the Toy Model continuously as
radiation measurements for the forest aren't continuous (see below). We
could use forcing data from the open site, but eh.
%}

%----------------------  forest site met variables  ----------------------%
%{
Hourly resolution for period 1 January 2004 0:00 to 31 December 2007 23:00.
%}
bla = xlsread('MeteoData_Alptal_Forested_Site.xls');
Wind_1h = bla(:,3);         % wind speed at 35m height [m s^{-1}]
T_soil_20_1h = bla(:,1);    % soil temperature at 20cm depth [°C]
    T_soil_20_1h = T_soil_20_1h + 273.15;
T_snow_10_1h = bla(:,4);    % snow temperature at 40cm height [°C]
    T_snow_10_1h = T_snow_10_1h + 273.15;
T_snow_40_1h = bla(:,5);    % snow temperature at 10cm height [°C]
    T_snow_40_1h = T_snow_40_1h + 273.15;
clear bla
bla = readtable('Alptal air temperature data.txt','Delimiter','\t','HeaderLines',6);
T_air_ac_1h = table2array(bla(:,8));    % above-canopy air temperature at 35m
T_air_ac_1h = T_air_ac_1h + 273.15;
clear bla

% manual snow measurements
bla = xlsread('SnowData_Alptal_Both_Sites.xlsx',1); % sheet 2 for open site
z_snow_forest = bla(:,7);   % sub-canopy snow depth [cm]
    z_snow_forest = z_snow_forest/100;
year = bla(:,1);
month = bla(:,2);
day = bla(:,3);
time_snow = nan(size(year));
for l=1:length(year)
    time_snow(l) = datenum(year(l),month(l),day(l),12,0,0);
end
clear year month day bla
    
% ground water level for 2004 (see SoilWaterProxy_Alptal.m)
bla = xlsread('data_Alptal_SNOWMIP2_200304.xls');
GWL_0304 = bla(:,13);
T_air_ac_0304 = bla(:,8);         % air temperature above the forest canopy [°C]
    T_air_ac_0304 = T_air_ac_0304 + 273.15;
T_open_0304 = bla(:,1);         % air temperature above the forest canopy [°C]
    T_open_0304 = T_open_0304 + 273.15;
clear bla

%--------------------  radiation data - rail & tower  --------------------%
% 2004
%{
Temporal resolution of 10 mins for period 23 January 6:30 to 19 July 10:00.
%}
bla = xlsread('AlptalRadiationData0304_calibrated_xp.xlsx');
SW_in_bc_2004 = bla(:,6);
SW_out_bc_2004 = bla(:,7);
SW_in_ac_2004 = bla(:,8);
SW_out_ac_2004 = bla(:,9);
LW_in_bc_2004 = bla(:,10);
LW_out_bc_2004 = bla(:,11);
LW_in_ac_2004 = bla(:,12);
LW_out_ac_2004 = bla(:,13);
year = bla(:,1);
month = bla(:,2);
day = bla(:,3);
hour = bla(:,4);
min = bla(:,5);
time_rad_2004 = nan(size(LW_in_bc_2004));
for l=1:length(LW_in_bc_2004)
    time_rad_2004(l) = datenum(year(l),month(l),day(l),hour(l),min(l),0);
end
clear hour month day hour min bla

% 2004/05
%{
Temporal resolution of 10 mins for period 14 October 0:10 to 10 June 6:10.
%}
bla = xlsread('AlptalRadiationData0405_calibrated_xp.xlsx');
SW_in_bc_2005 = bla(:,6);
SW_out_bc_2005 = bla(:,7);
SW_in_ac_2005 = bla(:,8);
SW_out_ac_2005 = bla(:,9);
LW_in_bc_2005 = bla(:,10);
LW_out_bc_2005 = bla(:,11);
LW_in_ac_2005 = bla(:,12);
    LW_in_ac_2005(14295:14311) = nan;   % gap 21 Jan 6:30 to 9:10
    % Can be filled by using open site LWR measurements.
LW_out_ac_2005 = bla(:,13);
year = bla(:,1);
month = bla(:,2);
day = bla(:,3);
hour = bla(:,4);
min = bla(:,5);
time_rad_2005 = nan(size(LW_in_bc_2005));
for l=1:length(LW_in_bc_2005)
    time_rad_2005(l) = datenum(year(l),month(l),day(l),hour(l),min(l),0);
end
clear hour month day hour min bla

% 2005/06
%{
Temporal resolution of 10 mins for period 15 October 14:00 to 16 May 6:10.
However, there's no SWR for 15 October, so we skip these timesteps and start
16 October 0:10.
%}
bla = xlsread('AlptalRadiationData0506_calibrated_xp.xlsx');
SW_in_bc_2006 = bla(62:end,6);
SW_out_bc_2006 = bla(62:end,7);
SW_in_ac_2006 = bla(62:end,8);
SW_out_ac_2006 = bla(62:end,9);
LW_in_bc_2006 = bla(62:end,10);
LW_out_bc_2006 = bla(62:end,11);
LW_in_ac_2006 = bla(62:end,12);
LW_out_ac_2006 = bla(62:end,13);
year = bla(62:end,1)-1900;   % typo: says "3905" and "3906"
month = bla(62:end,2);
day = bla(62:end,3);
hour = bla(62:end,4);
min = bla(62:end,5);
time_rad_2006 = nan(size(LW_in_bc_2006));
for l=1:length(LW_in_bc_2006)
    time_rad_2006(l) = datenum(year(l),month(l),day(l),hour(l),min(l),0);
end
clear hour month day hour min bla

% 2006/07
%{
Temporal resolution of 10 mins for period 27 November 15:10 to 30 April 6:10.
%}
bla = xlsread('AlptalRadiationData0607_calibrated_xp.xlsx');
SW_in_bc_2007 = bla(:,6);
SW_out_bc_2007 = bla(:,7);
SW_in_ac_2007 = bla(:,8);
    SW_in_ac_2007(2289:2417) = nan; % gap 13 Dec 12:30 to 14 Dec 9:50
    % We'll start simulation afterwards.
SW_out_ac_2007 = bla(:,9);
LW_in_bc_2007 = bla(:,10);
LW_out_bc_2007 = bla(:,11);
LW_in_ac_2007 = bla(:,12);
LW_out_ac_2007 = bla(:,13);
year = bla(:,1);
month = bla(:,2);
day = bla(:,3);
hour = bla(:,4);
min = bla(:,5);
time_rad_2007 = nan(size(LW_in_bc_2007));
for l=1:length(LW_in_bc_2007)
    time_rad_2007(l) = datenum(year(l),month(l),day(l),hour(l),min(l),0);
end
clear year month day hour min bla

%--------------------------  ground water level  --------------------------
%{
Two measurements of ground water level representing topography are
available, one for mounds on which the trees grow, and one for depressions
in between. The measurements for mounds are used to calculate a fraction of
the approximate range over the 4-year period.
Hourly measurements start at 0:00 on 1 January 2004 and are available until
23:00 on 31 May 2007, indicating that the start, not the end, of the
timestep is given.
%}
filename = 'C:\Users\w15044103\Documents\Modelling\Toy Model\Experiments\Alptal_extended\Alptal2004-2007.txt';
delimiter = '\t';
startRow = 2;
formatSpec = '%q%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
datetime_GWL = dataArray{:,1};
GWL_1h_depression = dataArray{:,2};
GWL_1h_mound = dataArray{:,3};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


%-------------------------------------------------------------------------%
%---------------------  cut out simulation periods  ----------------------%
%-------------------------------------------------------------------------%

%-----------------------------  season 2004  -----------------------------%
% 23 January 7:10/8:00 to 19 July 10:00
T_air_open_1h_2004 = T_air_open_1h(2745-1:7019-1);
RH_open_1h_2004 = RH_open_1h(2745-1:7019-1);
Precip_1h_2004 = Precip_1h(2745-1:7019-1);
Wind_1h_2004 = Wind_1h(537-1:4811-1);
T_snow_10_1h_2004 = T_snow_10_1h(537-1:4811-1);
T_snow_40_1h_2004 = T_snow_40_1h(537-1:4811-1);
T_soil_20_1h_2004 = T_soil_20_1h(537-1:4811-1);
T_air_ac_1h_2004 = T_air_ac_1h(537-1:4811-1);   % 1-hour shift necessary?
GWL_1h_depression_2004 = GWL_1h_depression(536:4810);   % 1-hour shift (s.a.)
GWL_1h_2004 = GWL_0304(4209-1:end);   % only last until 23 May
T_air_ac_1h_2004_SnowMIP = T_air_ac_0304(4209-1:end);   % only last until 23 May
T_open_1h_2004_SnowMIP = T_open_0304(4209-1:end);   % only last until 23 May
%{
There also seems to be a lead for SnowMIP2 data as air temperature peaks at
12:00 but should probably be later (as for modified open air temperature).
The whole things is weird, but insolation really shouldn't peak at 13:00.
%}

i=0;
for l=4+6:6:length(time_rad_2004)
    i=i+1;
    SW_in_ac_2004_1h(i) = mean(SW_in_ac_2004(l-5:l));
    SW_in_bc_2004_1h(i) = mean(SW_in_bc_2004(l-5:l));
    SW_out_bc_2004_1h(i) = mean(SW_out_bc_2004(l-5:l));
    LW_in_ac_2004_1h(i) = mean(LW_in_ac_2004(l-5:l));
    LW_in_bc_2004_1h(i) = mean(LW_in_bc_2004(l-5:l));
    LW_out_bc_2004_1h(i) = mean(LW_out_bc_2004(l-5:l));
    time_1h_2004(i) = time_rad_2004(l)-1/24;
end

alb_surf_2004_1h = nan(size(time_1h_2004));
for l=1:length(alb_surf_2004_1h)
    if SW_in_bc_2004_1h(l) > 1
        alb_surf_2004_1h(l) = SW_out_bc_2004_1h(l)/SW_in_bc_2004_1h(l);
    end
    if ~isnan(alb_surf_2004_1h(l))
        alb_surf_2004_1h(l) = min(alb_surf_2004_1h(l),1);
    end
end

%-----------------------------  season 2005  -----------------------------%
% 14 October 0:10/1:00 to 10 June 6:00
T_air_open_1h_2005 = T_air_open_1h(9098-1:14839-1);
RH_open_1h_2005 = RH_open_1h(9098-1:14839-1);
Precip_1h_2005 = Precip_1h(9098-1:14839-1);
LW_in_open_1h_2005 = LW_in_open_1h(9098-1:14839-1);
Wind_1h_2005 = Wind_1h(6890-1:12631-1);
T_snow_10_1h_2005 = T_snow_10_1h(6890-1:12631-1);
T_snow_40_1h_2005 = T_snow_40_1h(6890-1:12631-1);
T_soil_20_1h_2005 = T_soil_20_1h(6890-1:12631-1);
T_air_ac_1h_2005 = T_air_ac_1h(6890-1:12631-1);   % 1-hour shift necessary?
GWL_1h_depression_2005 = GWL_1h_depression(6889:12630);   % 1-hour shift (s.a.)

i=0;
for l=6:6:length(time_rad_2005)-1
    i=i+1;
    SW_in_ac_2005_1h(i) = mean(SW_in_ac_2005(l-5:l));
    SW_in_bc_2005_1h(i) = mean(SW_in_bc_2005(l-5:l));
    SW_out_bc_2005_1h(i) = mean(SW_out_bc_2005(l-5:l));
    LW_in_ac_2005_1h(i) = mean(LW_in_ac_2005(l-5:l));
    LW_in_bc_2005_1h(i) = mean(LW_in_bc_2005(l-5:l));
    LW_out_bc_2005_1h(i) = mean(LW_out_bc_2005(l-5:l));
    time_1h_2005(i) = time_rad_2005(l)-1/24;
end

alb_surf_2005_1h = nan(size(time_1h_2005));
for l=1:length(alb_surf_2005_1h)
    if SW_in_bc_2005_1h(l) > 1
        alb_surf_2005_1h(l) = SW_out_bc_2005_1h(l)/SW_in_bc_2005_1h(l);
    end
    if ~isnan(alb_surf_2005_1h(l))
        alb_surf_2005_1h(l) = min(alb_surf_2005_1h(l),1);
    end
end

% fill 4h gap from open site data
scaling = (LW_in_ac_2005_1h(2387)-LW_in_ac_2005_1h(2382))...
    /(LW_in_open_1h_2005(2387)-LW_in_open_1h_2005(2382));
for t=2383:2386
    LW_in_ac_2005_1h(t) = LW_in_ac_2005_1h(t-1) + scaling*...
        (LW_in_open_1h_2005(t)-LW_in_open_1h_2005(t-1));
end

%-----------------------------  season 2006  -----------------------------%
% 16 October 0:10/1:00 to 16 May 6:00
T_air_open_1h_2006 = T_air_open_1h(17906-1:22999-1);
RH_open_1h_2006 = RH_open_1h(17906-1:22999-1);
Precip_1h_2006 = Precip_1h(17906-1:22999-1);
Wind_1h_2006 = Wind_1h(15698-1:20791-1);
T_snow_10_1h_2006 = T_snow_10_1h(15698-1:20791-1);
T_snow_40_1h_2006 = T_snow_40_1h(15698-1:20791-1);
T_soil_20_1h_2006 = T_soil_20_1h(15698-1:20791-1);
T_air_ac_1h_2006 = T_air_ac_1h(15698-1:20791-1);   % 1-hour shift necessary?
GWL_1h_depression_2006 = GWL_1h_depression(15697:20790);   % 1-hour shift (s.a.)

i=0;
for l=6:6:length(time_rad_2006)-1
    i=i+1;
    SW_in_ac_2006_1h(i) = mean(SW_in_ac_2006(l-5:l));
    SW_in_bc_2006_1h(i) = mean(SW_in_bc_2006(l-5:l));
    SW_out_bc_2006_1h(i) = mean(SW_out_bc_2006(l-5:l));
    LW_in_ac_2006_1h(i) = mean(LW_in_ac_2006(l-5:l));
    LW_in_bc_2006_1h(i) = mean(LW_in_bc_2006(l-5:l));
    LW_out_bc_2006_1h(i) = mean(LW_out_bc_2006(l-5:l));
    time_1h_2006(i) = time_rad_2006(l)-1/24;
end

alb_surf_2006_1h = nan(size(time_1h_2006));
for l=1:length(alb_surf_2006_1h)
    if SW_in_bc_2006_1h(l) > 1
        alb_surf_2006_1h(l) = SW_out_bc_2006_1h(l)/SW_in_bc_2006_1h(l);
    end
    if ~isnan(alb_surf_2006_1h(l))
        alb_surf_2006_1h(l) = min(alb_surf_2006_1h(l),1);
    end
end

%-----------------------------  season 2007  -----------------------------%
% 27 November 15:10/16:00 to 30 April 6:00
T_air_open_1h_2007 = T_air_open_1h(27689-1:31375-1);
RH_open_1h_2007 = RH_open_1h(27689-1:31375-1);
Precip_1h_2007 = Precip_1h(27689-1:31375-1);
Wind_1h_2007 = Wind_1h(25481-1:29167-1);
T_snow_10_1h_2007 = T_snow_10_1h(25481-1:29167-1);
T_snow_40_1h_2007 = T_snow_40_1h(25481-1:29167-1);
T_soil_20_1h_2007 = T_soil_20_1h(25481-1:29167-1);
T_air_ac_1h_2007 = T_air_ac_1h(25481-1:29167-1);   % 1-hour shift necessary?
GWL_1h_depression_2007 = GWL_1h_depression(25480:29166);   % 1-hour shift (s.a.)

i=0;
for l=6:6:length(time_rad_2007)-1
    i=i+1;
    SW_in_ac_2007_1h(i) = mean(SW_in_ac_2007(l-5:l));
    SW_in_bc_2007_1h(i) = mean(SW_in_bc_2007(l-5:l));
    SW_out_bc_2007_1h(i) = mean(SW_out_bc_2007(l-5:l));
    LW_in_ac_2007_1h(i) = mean(LW_in_ac_2007(l-5:l));
    LW_in_bc_2007_1h(i) = mean(LW_in_bc_2007(l-5:l));
    LW_out_bc_2007_1h(i) = mean(LW_out_bc_2007(l-5:l));
    time_1h_2007(i) = time_rad_2007(l)-1/24;
end

alb_surf_2007_1h = nan(size(time_1h_2007));
for l=1:length(alb_surf_2007_1h)
    if SW_in_bc_2007_1h(l) > 1
        alb_surf_2007_1h(l) = SW_out_bc_2007_1h(l)/SW_in_bc_2007_1h(l);
    end
    if ~isnan(alb_surf_2007_1h(l))
        alb_surf_2007_1h(l) = min(alb_surf_2007_1h(l),1);
    end
end


%-------------------------------------------------------------------------%
%--------------------  check for snow on radiometers  --------------------%
%-------------------------------------------------------------------------%
% 2004
LW_in_bc_2004_1h(42+1:65+1) = nan;
LW_in_bc_2004_1h(450+1:473+1) = nan;
LW_in_bc_2004_1h(642+1:665+1) = nan;
LW_in_bc_2004_1h(1074+1:1169+1) = nan;  % two separate snowfall events
LW_in_bc_2004_1h(1434+1:1553+1) = nan;  % seems to have affected above canopy as well
LW_in_bc_2004_1h(1770+1:1817+1) = nan;  % might have affected above canopy as well

% 2005
LW_in_bc_2005_1h(746+1:768+1) = nan;    % yup, heavey snowfall
    % Might last until 770, but day starts at 770. But November anyway.
LW_in_bc_2005_1h(2425+1:2496+1) = nan;  % several snowfall events, better be safe
LW_in_bc_2005_1h(2641+1:2664+1) = nan;
LW_in_bc_2005_1h(2929+1:3144+1) = nan;  % several snowfall events
LW_in_bc_2005_1h(3481+1:3528+1) = nan;

% 2006
LW_in_bc_2006_1h(1231+1:1246+1) = nan;  % December anyway, so doesn't matter
LW_in_bc_2006_1h(2777+1:2834+1) = nan;  % yup, heavy snowfall

% 2007
LW_in_bc_2007_1h(856+1:877+1) = nan;    % yup, heavy snowfall
LW_in_bc_2007_1h(2689+1:2805+1) = nan;  % yup, heavy snowfall


%-------------------------------------------------------------------------%
%-----------  extrapolate snow depth from manual measurements  -----------%
%-------------------------------------------------------------------------%

% 1) set manual measurements
z_snow_forest_1h_2004 = nan(size(time_1h_2004));
for i=1:18
    bla = find(time_1h_2004==time_snow(i));
    z_snow_forest_1h_2004(bla) = z_snow_forest(i);
end
z_snow_forest_1h_2005 = nan(size(time_1h_2005));
for i=19:31
    bla = find(time_1h_2005==time_snow(i));
    z_snow_forest_1h_2005(bla) = z_snow_forest(i);
end
z_snow_forest_1h_2006 = nan(size(time_1h_2006));
for i=32:51
    bla = find(time_1h_2006==time_snow(i));
    z_snow_forest_1h_2006(bla) = z_snow_forest(i);
end
z_snow_forest_1h_2007 = nan(size(time_1h_2007));
for i=52:58
    bla = find(time_1h_2007==time_snow(i));
    z_snow_forest_1h_2007(bla) = z_snow_forest(i);
end

% 2) interpolate between measurements based on educated guesses - air
    % temperature if snow depth > 40cm
% 2004
T_pos = nan(size(T_snow_40_1h_2004));
T_pos_air = nan(size(T_air_open_1h_2004));
SnowFall = nan(size(Precip_1h_2004));
for t=1:length(z_snow_forest_1h_2004)
    T_pos(t) = max(T_snow_40_1h_2004(t)-273.15,0);
    T_pos_air(t) = max(T_air_open_1h_2004(t)-273.15,0);
    f_snow = RainSnowPartitioningAlptal(T_air_open_1h_2004(t));
    SnowFall(t) = f_snow*Precip_1h_2004(t);
end
bla = find(~isnan(z_snow_forest_1h_2004));
for i=1:length(bla)-1
    if z_snow_forest_1h_2004(bla(i)) < z_snow_forest_1h_2004(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            z_snow_forest_1h_2004(b) = z_snow_forest_1h_2004(b-1)...
                + (z_snow_forest_1h_2004(bla(i+1))-z_snow_forest_1h_2004(bla(i)))*...
                SnowFall(b)/sum(SnowFall(bla(i):bla(i+1)-1));
        end
    elseif z_snow_forest_1h_2004(bla(i)) > z_snow_forest_1h_2004(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            if sum(T_pos(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2004(b) = z_snow_forest_1h_2004(b-1)...
                    - (z_snow_forest_1h_2004(bla(i))-z_snow_forest_1h_2004(bla(i+1)))*...
                    T_pos(b)/sum(T_pos(bla(i):bla(i+1)-1));
            elseif sum(T_pos_air(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2004(b) = z_snow_forest_1h_2004(b-1)...
                    - (z_snow_forest_1h_2004(bla(i))-z_snow_forest_1h_2004(bla(i+1)))*...
                    T_pos_air(b)/sum(T_pos_air(bla(i):bla(i+1)-1));
            else
                z_snow_forest_1h_2004(b) = z_snow_forest_1h_2004(b-1)...
                    - (z_snow_forest_1h_2004(bla(i))-z_snow_forest_1h_2004(bla(i+1)))...
                    /(bla(i+1)-bla(i));
            end
        end
    elseif z_snow_forest_1h_2004(bla(i)) == z_snow_forest_1h_2004(bla(i+1))
        z_snow_forest_1h_2004(bla(i)+1:bla(i+1)-1) = z_snow_forest_1h_2004(bla(i));
    end
end

% 2005
T_pos = nan(size(T_snow_40_1h_2005));
T_pos_air = nan(size(T_air_open_1h_2005));
SnowFall = nan(size(Precip_1h_2005));
for t=1:length(z_snow_forest_1h_2005)
    T_pos(t) = max(T_snow_40_1h_2005(t)-273.15,0);
    T_pos_air(t) = max(T_air_open_1h_2005(t)-273.15,0);
    f_snow = RainSnowPartitioningAlptal(T_air_open_1h_2005(t));
    SnowFall(t) = f_snow*Precip_1h_2005(t);
end
bla = find(~isnan(z_snow_forest_1h_2005));
for i=1:length(bla)-1
    if z_snow_forest_1h_2005(bla(i)) < z_snow_forest_1h_2005(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            z_snow_forest_1h_2005(b) = z_snow_forest_1h_2005(b-1)...
                + (z_snow_forest_1h_2005(bla(i+1))-z_snow_forest_1h_2005(bla(i)))*...
                SnowFall(b)/sum(SnowFall(bla(i):bla(i+1)-1));
        end
    elseif z_snow_forest_1h_2005(bla(i)) > z_snow_forest_1h_2005(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            if sum(T_pos(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2005(b) = z_snow_forest_1h_2005(b-1)...
                    - (z_snow_forest_1h_2005(bla(i))-z_snow_forest_1h_2005(bla(i+1)))*...
                    T_pos(b)/sum(T_pos(bla(i):bla(i+1)-1));
            elseif sum(T_pos_air(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2005(b) = z_snow_forest_1h_2005(b-1)...
                    - (z_snow_forest_1h_2005(bla(i))-z_snow_forest_1h_2005(bla(i+1)))*...
                    T_pos_air(b)/sum(T_pos_air(bla(i):bla(i+1)-1));
            else
                z_snow_forest_1h_2005(b) = z_snow_forest_1h_2005(b-1)...
                    - (z_snow_forest_1h_2005(bla(i))-z_snow_forest_1h_2005(bla(i+1)))...
                    /(bla(i+1)-bla(i));
            end
        end
    elseif z_snow_forest_1h_2005(bla(i)) == z_snow_forest_1h_2005(bla(i+1))
        z_snow_forest_1h_2005(bla(i)+1:bla(i+1)-1) = z_snow_forest_1h_2005(bla(i));
    end
end

% 2006
T_pos = nan(size(T_snow_40_1h_2006));
T_pos_air = nan(size(T_air_open_1h_2006));
SnowFall = nan(size(Precip_1h_2006));
for t=1:length(z_snow_forest_1h_2006)
    T_pos(t) = max(T_snow_40_1h_2006(t)-273.15,0);
    T_pos_air(t) = max(T_air_open_1h_2006(t)-273.15,0);
    f_snow = RainSnowPartitioningAlptal(T_air_open_1h_2006(t));
    SnowFall(t) = f_snow*Precip_1h_2006(t);
end
bla = find(~isnan(z_snow_forest_1h_2006));
for i=1:length(bla)-1
    if z_snow_forest_1h_2006(bla(i)) < z_snow_forest_1h_2006(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            z_snow_forest_1h_2006(b) = z_snow_forest_1h_2006(b-1)...
                + (z_snow_forest_1h_2006(bla(i+1))-z_snow_forest_1h_2006(bla(i)))*...
                SnowFall(b)/sum(SnowFall(bla(i):bla(i+1)-1));
        end
    elseif z_snow_forest_1h_2006(bla(i)) > z_snow_forest_1h_2006(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            if sum(T_pos(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2006(b) = z_snow_forest_1h_2006(b-1)...
                    - (z_snow_forest_1h_2006(bla(i))-z_snow_forest_1h_2006(bla(i+1)))*...
                    T_pos(b)/sum(T_pos(bla(i):bla(i+1)-1));
            elseif sum(T_pos_air(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2006(b) = z_snow_forest_1h_2006(b-1)...
                    - (z_snow_forest_1h_2006(bla(i))-z_snow_forest_1h_2006(bla(i+1)))*...
                    T_pos_air(b)/sum(T_pos_air(bla(i):bla(i+1)-1));
            else
                z_snow_forest_1h_2006(b) = z_snow_forest_1h_2006(b-1)...
                    - (z_snow_forest_1h_2006(bla(i))-z_snow_forest_1h_2006(bla(i+1)))...
                    /(bla(i+1)-bla(i));
            end
        end
    elseif z_snow_forest_1h_2006(bla(i)) == z_snow_forest_1h_2006(bla(i+1))
        z_snow_forest_1h_2006(bla(i)+1:bla(i+1)-1) = z_snow_forest_1h_2006(bla(i));
    end
end

% 2007
T_pos = nan(size(T_snow_40_1h_2007));
T_pos_air = nan(size(T_air_open_1h_2007));
SnowFall = nan(size(Precip_1h_2007));
for t=1:length(z_snow_forest_1h_2007)
    T_pos(t) = max(T_snow_40_1h_2007(t)-273.15,0);
    T_pos_air(t) = max(T_air_open_1h_2007(t)-273.15,0);
    f_snow = RainSnowPartitioningAlptal(T_air_open_1h_2007(t));
    SnowFall(t) = f_snow*Precip_1h_2007(t);
end
bla = find(~isnan(z_snow_forest_1h_2007));
for i=1:length(bla)-1
    if z_snow_forest_1h_2007(bla(i)) < z_snow_forest_1h_2007(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            z_snow_forest_1h_2007(b) = z_snow_forest_1h_2007(b-1)...
                + (z_snow_forest_1h_2007(bla(i+1))-z_snow_forest_1h_2007(bla(i)))*...
                SnowFall(b)/sum(SnowFall(bla(i):bla(i+1)-1));
        end
    elseif z_snow_forest_1h_2007(bla(i)) > z_snow_forest_1h_2007(bla(i+1))
        for b=bla(i)+1:bla(i+1)-1
            if sum(T_pos(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2007(b) = z_snow_forest_1h_2007(b-1)...
                    - (z_snow_forest_1h_2007(bla(i))-z_snow_forest_1h_2007(bla(i+1)))*...
                    T_pos(b)/sum(T_pos(bla(i):bla(i+1)-1));
            elseif sum(T_pos_air(bla(i):bla(i+1)-1)) > 0
                z_snow_forest_1h_2007(b) = z_snow_forest_1h_2007(b-1)...
                    - (z_snow_forest_1h_2007(bla(i))-z_snow_forest_1h_2007(bla(i+1)))*...
                    T_pos_air(b)/sum(T_pos_air(bla(i):bla(i+1)-1));
            else
                z_snow_forest_1h_2007(b) = z_snow_forest_1h_2007(b-1)...
                    - (z_snow_forest_1h_2007(bla(i))-z_snow_forest_1h_2007(bla(i+1)))...
                    /(bla(i+1)-bla(i));
            end
        end
    elseif z_snow_forest_1h_2007(bla(i)) == z_snow_forest_1h_2007(bla(i+1))
        z_snow_forest_1h_2007(bla(i)+1:bla(i+1)-1) = z_snow_forest_1h_2007(bla(i));
    end
end


%-------------------------------------------------------------------------%
%----------------------  adjust simulation periods  ----------------------%
%-------------------------------------------------------------------------%
%{
We will check for likely meltout dates using outgoing sub-canopy LWR and
soil temperature, or restrict by manual snow measurements, to limit
evaluation period. For evaluation period start, we choose (arbitrarily)
1 January, except for 2004 when measurements only start on 23 January.
As the Toy Model doesn't need much spin-up time, we start all simulations
on 31 December (except for 2004) and run until the end of DoY 120, or early
morning for 2007, which is 29 April for 2004 and 30 April for 2005 - 2007.
Meltout happened in each year prior to that day and the last manual snow
measurements were conducted in early or mid-April.
%}
simperiod_2004 = 1:2346;    start_2004 = 121*24-length(simperiod_2004)+1;
simperiod_2005 = 1874:4777;
simperiod_2006 = 1826:4729;
simperiod_2007 = 803:length(time_1h_2007);

Precip_1h_all = nan(121*24,4);
Precip_1h_all(start_2004:end,1) ...
    = Precip_1h_2004(simperiod_2004(1):simperiod_2004(end));
Precip_1h_all(:,2) = Precip_1h_2005(simperiod_2005(1):simperiod_2005(end));
Precip_1h_all(:,3) = Precip_1h_2006(simperiod_2006(1):simperiod_2006(end));
Precip_1h_all(1:length(simperiod_2007),4) ...
    = Precip_1h_2007(simperiod_2007(1):simperiod_2007(end));

RH_open_1h_all = nan(121*24,4);
RH_open_1h_all(start_2004:end,1) ...
    = RH_open_1h_2004(simperiod_2004(1):simperiod_2004(end));
RH_open_1h_all(:,2) = RH_open_1h_2005(simperiod_2005(1):simperiod_2005(end));
RH_open_1h_all(:,3) = RH_open_1h_2006(simperiod_2006(1):simperiod_2006(end));
RH_open_1h_all(1:length(simperiod_2007),4) ...
    = RH_open_1h_2007(simperiod_2007(1):simperiod_2007(end));

Wind_1h_all = nan(121*24,4);
Wind_1h_all(start_2004:end,1) ...
    = Wind_1h_2004(simperiod_2004(1):simperiod_2004(end));
Wind_1h_all(:,2) = Wind_1h_2005(simperiod_2005(1):simperiod_2005(end));
Wind_1h_all(:,3) = Wind_1h_2006(simperiod_2006(1):simperiod_2006(end));
Wind_1h_all(1:length(simperiod_2007),4) ...
    = Wind_1h_2007(simperiod_2007(1):simperiod_2007(end));

T_air_ac_1h_all = nan(121*24,4);
T_air_ac_1h_2004_SnowMIP_cut = nan(121*24,1);
T_air_ac_1h_2004_SnowMIP_cut(start_2004:end) ...
    = T_air_ac_1h_2004_SnowMIP(simperiod_2004(1):simperiod_2004(end));
T_open_1h_2004_SnowMIP_cut = nan(121*24,1);
T_open_1h_2004_SnowMIP_cut(start_2004:end) ...
    = T_open_1h_2004_SnowMIP(simperiod_2004(1):simperiod_2004(end));
T_air_ac_1h_all(start_2004:end,1) ...
    = T_air_ac_1h_2004(simperiod_2004(1):simperiod_2004(end));
T_air_ac_1h_all(:,2) = T_air_ac_1h_2005(simperiod_2005(1):simperiod_2005(end));
T_air_ac_1h_all(:,3) = T_air_ac_1h_2006(simperiod_2006(1):simperiod_2006(end));
T_air_ac_1h_all(1:length(simperiod_2007),4) ...
    = T_air_ac_1h_2007(simperiod_2007(1):simperiod_2007(end));

T_air_open_1h_all = nan(121*24,4);
T_air_open_1h_all(start_2004:end,1) ...
    = T_air_open_1h_2004(simperiod_2004(1):simperiod_2004(end));
%    T_air_open_1h_2004_cmpac = nan(121*24,1);
%    T_air_open_1h_2004_cmpac(start_2004:end) ...
%       = T_air_ac_1h_2004(simperiod_2004(1):simperiod_2004(end));
T_air_open_1h_all(:,2) = T_air_open_1h_2005(simperiod_2005(1):simperiod_2005(end));
T_air_open_1h_all(:,3) = T_air_open_1h_2006(simperiod_2006(1):simperiod_2006(end));
T_air_open_1h_all(1:length(simperiod_2007),4) ...
    = T_air_open_1h_2007(simperiod_2007(1):simperiod_2007(end));

LW_in_ac_1h_all = nan(121*24,4);
LW_in_ac_1h_all(start_2004:end,1) ...
    = LW_in_ac_2004_1h(simperiod_2004(1):simperiod_2004(end));
LW_in_ac_1h_all(:,2) = LW_in_ac_2005_1h(simperiod_2005(1):simperiod_2005(end));
LW_in_ac_1h_all(:,3) = LW_in_ac_2006_1h(simperiod_2006(1):simperiod_2006(end));
LW_in_ac_1h_all(1:length(simperiod_2007),4) ...
    = LW_in_ac_2007_1h(simperiod_2007(1):simperiod_2007(end));

SW_in_ac_1h_all = nan(121*24,4);
SW_in_ac_1h_all(start_2004:end,1) ...
    = SW_in_ac_2004_1h(simperiod_2004(1):simperiod_2004(end));
SW_in_ac_1h_all(:,2) = SW_in_ac_2005_1h(simperiod_2005(1):simperiod_2005(end));
SW_in_ac_1h_all(:,3) = SW_in_ac_2006_1h(simperiod_2006(1):simperiod_2006(end));
SW_in_ac_1h_all(1:length(simperiod_2007),4) ...
    = SW_in_ac_2007_1h(simperiod_2007(1):simperiod_2007(end));

alb_surf_1h_all = nan(121*24,4);
alb_surf_1h_all(start_2004:end,1) ...
    = alb_surf_2004_1h(simperiod_2004(1):simperiod_2004(end));
alb_surf_1h_all(:,2) = alb_surf_2005_1h(simperiod_2005(1):simperiod_2005(end));
alb_surf_1h_all(:,3) = alb_surf_2006_1h(simperiod_2006(1):simperiod_2006(end));
alb_surf_1h_all(1:length(simperiod_2007),4) ...
    = alb_surf_2007_1h(simperiod_2007(1):simperiod_2007(end));

T_soil_20_1h_all = nan(121*24,4);
T_soil_20_1h_all(start_2004:end,1) ...
    = T_soil_20_1h_2004(simperiod_2004(1):simperiod_2004(end));
T_soil_20_1h_all(:,2) = T_soil_20_1h_2005(simperiod_2005(1):simperiod_2005(end));
T_soil_20_1h_all(:,3) = T_soil_20_1h_2006(simperiod_2006(1):simperiod_2006(end));
T_soil_20_1h_all(1:length(simperiod_2007),4) ...
    = T_soil_20_1h_2007(simperiod_2007(1):simperiod_2007(end));

T_snow_10_1h_all = nan(121*24,4);
T_snow_10_1h_all(start_2004:end,1) ...
    = T_snow_10_1h_2004(simperiod_2004(1):simperiod_2004(end));
T_snow_10_1h_all(:,2) = T_snow_10_1h_2005(simperiod_2005(1):simperiod_2005(end));
T_snow_10_1h_all(:,3) = T_snow_10_1h_2006(simperiod_2006(1):simperiod_2006(end));
T_snow_10_1h_all(1:length(simperiod_2007),4) ...
    = T_snow_10_1h_2007(simperiod_2007(1):simperiod_2007(end));

T_snow_40_1h_all = nan(121*24,4);
T_snow_40_1h_all(start_2004:end,1) ...
    = T_snow_40_1h_2004(simperiod_2004(1):simperiod_2004(end));
T_snow_40_1h_all(:,2) = T_snow_40_1h_2005(simperiod_2005(1):simperiod_2005(end));
T_snow_40_1h_all(:,3) = T_snow_40_1h_2006(simperiod_2006(1):simperiod_2006(end));
T_snow_40_1h_all(1:length(simperiod_2007),4) ...
    = T_snow_40_1h_2007(simperiod_2007(1):simperiod_2007(end));

z_snow_forest_1h_all = nan(121*24,4);
z_snow_forest_1h_all(start_2004:end,1) ...
    = z_snow_forest_1h_2004(simperiod_2004(1):simperiod_2004(end));
    % first manual snow depth measurement too late
    bla = find(~isnan(z_snow_forest_1h_all(:,1)),1);
    z_snow_forest_1h_all(1:bla-1,1) = z_snow_forest_1h_all(bla,1);
z_snow_forest_1h_all(:,2) = z_snow_forest_1h_2005(simperiod_2005(1):simperiod_2005(end));
z_snow_forest_1h_all(:,3) = z_snow_forest_1h_2006(simperiod_2006(1):simperiod_2006(end));
z_snow_forest_1h_all(1:length(simperiod_2007),4) ...
    = z_snow_forest_1h_2007(simperiod_2007(1):simperiod_2007(end));

LW_out_bc_1h_all = nan(121*24,4);
LW_out_bc_1h_all(start_2004:end,1) ...
    = LW_out_bc_2004_1h(simperiod_2004(1):simperiod_2004(end));
LW_out_bc_1h_all(:,2) = LW_out_bc_2005_1h(simperiod_2005(1):simperiod_2005(end));
LW_out_bc_1h_all(:,3) = LW_out_bc_2006_1h(simperiod_2006(1):simperiod_2006(end));
LW_out_bc_1h_all(1:length(simperiod_2007),4) ...
    = LW_out_bc_2007_1h(simperiod_2007(1):simperiod_2007(end));

LW_in_bc_1h_all = nan(121*24,4);
LW_in_bc_1h_all(start_2004:end,1) ...
    = LW_in_bc_2004_1h(simperiod_2004(1):simperiod_2004(end));
LW_in_bc_1h_all(:,2) = LW_in_bc_2005_1h(simperiod_2005(1):simperiod_2005(end));
LW_in_bc_1h_all(:,3) = LW_in_bc_2006_1h(simperiod_2006(1):simperiod_2006(end));
LW_in_bc_1h_all(1:length(simperiod_2007),4) ...
    = LW_in_bc_2007_1h(simperiod_2007(1):simperiod_2007(end));

time_1h_all = nan(121*24,4);
time_1h_all(start_2004:end,1) ...
    = time_1h_2004(simperiod_2004(1):simperiod_2004(end));
time_1h_all(:,2) = time_1h_2005(simperiod_2005(1):simperiod_2005(end));
time_1h_all(:,3) = time_1h_2006(simperiod_2006(1):simperiod_2006(end));
time_1h_all(1:length(simperiod_2007),4) ...
    = time_1h_2007(simperiod_2007(1):simperiod_2007(end));


%-------------------------------------------------------------------------%
%-------------  calculate variables  required for Toy Model  -------------%
% diffuse fraction via effective emissivity of the sky
em_sky_1h_all = nan(size(T_air_open_1h_all));
for t=1:length(em_sky_1h_all(:,1))
    for i=1:4
        em_sky_1h_all(t,i) = LW_in_ac_1h_all(t,i)/(5.67*10^(-8)...
            *T_air_open_1h_all(t,i)^4);
    end
end
em_sky_min = min([min(em_sky_1h_all(:,1)) min(em_sky_1h_all(:,2))...
    min(em_sky_1h_all(:,3)) min(em_sky_1h_all(:,4))]);
em_sky_max = max([max(em_sky_1h_all(:,1)) max(em_sky_1h_all(:,2))...
    max(em_sky_1h_all(:,3)) max(em_sky_1h_all(:,4))]);

frac_diff_1h_all = nan(size(em_sky_1h_all));
for t=1:length(frac_diff_1h_all(:,1))
    for i=1:4
        frac_diff_1h_all(t,i) = (em_sky_1h_all(t,i)-em_sky_min)...
            /(em_sky_max-em_sky_min);
    end
end

% fraction of ground covered by snow
%{
Since outgoing LWR corresponds to surface temperatures of distinctly more
than 0°C although there's still snow on the ground, snow cover fraction
gets adjusted to 0.5 (although it's probably understory vegetation and not
soil). This is done for snow depth of 15cm and less as this values seems
like a good fit based on surface albedo measurements.
%}
frac_snow_1h_all = nan(size(z_snow_forest_1h_all));
for t=1:length(frac_snow_1h_all(:,1))
    for i=1:4
        if z_snow_forest_1h_all(t,i) > 0.15
            frac_snow_1h_all(t,i) = 1;
        elseif z_snow_forest_1h_all(t,i) <= 0.15 && z_snow_forest_1h_all(t,i) > 0
            frac_snow_1h_all(t,i) = 0.5;
        elseif z_snow_forest_1h_all(t,i) == 0
            frac_snow_1h_all(t,i) = 0;
        end
    end
end

% surface temperature from outgoing sub-canopy LWR
T_surf_1h_all = nan(size(time_1h_all));
for t=1:length(time_1h_all(:,1))
    for i=1:4
        em_gr = frac_snow_1h_all(t,i)*0.97 + (1-frac_snow_1h_all(t,i))*0.96;
        T_surf_1h_all(t,i) = nthroot(LW_out_bc_1h_all(t,i)/(5.67*10^(-8)*em_gr),4);
    end
end

% create soil water proxy (see SoilWaterProxy_Alptal.m)
GWLmax = 700;
GWLmin = 0;
GWL_1h_depression_2004_scaled = (GWL_1h_depression_2004-GWLmin)/(GWLmax-GWLmin);
SWC_proxy_1h_2004 = 1-GWL_1h_depression_2004_scaled;
GWL_1h_depression_2005_scaled = (GWL_1h_depression_2005-GWLmin)/(GWLmax-GWLmin);
SWC_proxy_1h_2005 = 1-GWL_1h_depression_2005_scaled;
GWL_1h_depression_2006_scaled = (GWL_1h_depression_2006-GWLmin)/(GWLmax-GWLmin);
SWC_proxy_1h_2006 = 1-GWL_1h_depression_2006_scaled;
GWL_1h_depression_2007_scaled = (GWL_1h_depression_2007-GWLmin)/(GWLmax-GWLmin);
SWC_proxy_1h_2007 = 1-GWL_1h_depression_2007_scaled;

SWProxy_1h_all = nan(121*24,4);
SWProxy_1h_all(start_2004:end,1) ...
    = SWC_proxy_1h_2004(simperiod_2004(1):simperiod_2004(end));
SWProxy_1h_all(:,2) = SWC_proxy_1h_2005(simperiod_2005(1):simperiod_2005(end));
SWProxy_1h_all(:,3) = SWC_proxy_1h_2006(simperiod_2006(1):simperiod_2006(end));
SWProxy_1h_all(1:length(simperiod_2007),4) ...
    = SWC_proxy_1h_2007(simperiod_2007(1):simperiod_2007(end));


%-------------------------------------------------------------------------%
%---------------------------  Note and Saving  ---------------------------%
%{
I'm generally assuming a given timestamp refers to the mean (or sum) over
the previous timestep. Therefore, radiation data are averaged over the
previous 10 minutes and hourly averages are set to the last timestep used
for the average. For hourly values (open or forest), this would mean
averages over the previous hour. However, open site radiation leads
forest radiometer data by 1 hour, indicating that hourly values are set to
the beginning of the averaged period. This would also explain why
measurements start at 0:00 and end at 23:00.
Then again, insolation for forest radiometers peaks at 13:00 while open
site insolation peaks at 12:00 assuming a timestamp refers to the end of
the averaging period.
Therefore, we "correct" time and open site forcing to fit radiometer data,
but peaks set to 12:00 (as we need time for solar angle calculation).
%}
save('ForcingData_ToyModel_Alptal.mat','time_1h_all','start_2004',...
    'Precip_1h_all','RH_open_1h_all','Wind_1h_all','T_air_ac_1h_all',...
    'LW_in_ac_1h_all','frac_diff_1h_all','SW_in_ac_1h_all','alb_surf_1h_all',...
    'z_snow_forest_1h_all','frac_snow_1h_all','T_surf_1h_all','SWProxy_1h_all',...
    'T_soil_20_1h_all','LW_out_bc_1h_all','LW_in_bc_1h_all','T_air_ac_1h_2004_SnowMIP_cut')

toc