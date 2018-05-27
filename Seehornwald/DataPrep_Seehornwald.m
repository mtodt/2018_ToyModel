tic
clear
close all
%--------------------------------------------------------------------------
%------------------------------  load data  -------------------------------
%--------------------------------------------------------------------------

%----------------------  NABEL meteorological data  ----------------------%
Met = ...
    readtable('Data_Meteo_Oct07_Apr14_raw.txt','Delimiter','\t');
time_met = table2array(Met(:,1));
T_air = table2array(Met(:,2));
RH = table2array(Met(:,3));
% Rad_net = table2array(Met(:,4));
SW_in = table2array(Met(:,5));  % why negative?
P_air = table2array(Met(:,6));
Prec = table2array(Met(:,7));
% Wind_dir = table2array(Met(:,8));
Wind_speed = table2array(Met(:,9));
% Wind_speed_max = table2array(Met(:,10));

%-------------------  radiation data - rail and tower  -------------------%
Forest = ...
    readtable('Rabest_AllYears_071012_130329_600sec_LWcalib_SWcalib.txt',...
    'Delimiter','\t');
time_for = table2array(Forest(:,1));
SW_in_bc = table2array(Forest(:,2));
SW_out_bc = table2array(Forest(:,3));
LW_in_bc = table2array(Forest(:,4));
LW_out_bc = table2array(Forest(:,5));
SW_in_ac = table2array(Forest(:,7));
SW_out_ac = table2array(Forest(:,8));
LW_in_ac = table2array(Forest(:,9));
LW_out_ac = table2array(Forest(:,10));
% T_CNR1_bc = table2array(Forest(:,6));
% T_CNR1_ac = table2array(Forest(:,11));
% T_rail = table2array(Forest(:,13));
%{
Offset between start of measurement files about 11.5 days. Meteo measurements
start on 1 October 2007 at 0:10, while radiation measurements start on 12
October 2007 at 13:50. Also, meteo measurements end on 30 April 2014, while
radiation measurements end on 29 March 2013.
-> For clean start, both measurements cut at 13 October 2007 0:10.
10-minute intervals likely unnecessary - although interesting for
variations of LW enhancement - so that hourly averages are calculated,
which also cancels out changes in SVF due to movement along the rail.
Furthermore, this facilitates comparison to Alptal data.
Cut down to period:
13 October 2007 1:00 - 29 March 2013 0:00 (averaged for previous hour)
%}

time_met_for = time_met(1667:end-57218);
T_air_for = T_air(1667:end-57218);
RH_for = RH(1667:end-57218);
SW_in_for = SW_in(1667:end-57218);
P_air_for = P_air(1667:end-57218);
Prec_for = Prec(1667:end-57218);
Wind_speed_for = Wind_speed(1667:end-57218);

i=0;
time_1h = nan(47856,1);
T_air_ac_1h = nan(47856,1);
RH_ac_1h = nan(47856,1);
P_air_ac_1h = nan(47856,1);
Prec_ac_1h = nan(47856,1);
Wind_ac_1h = nan(47856,1);
SW_in_bc_1h = nan(47856,1);
SW_out_bc_1h = nan(47856,1);
LW_in_bc_1h = nan(47856,1);
LW_out_bc_1h = nan(47856,1);
SW_in_ac_1h = nan(47856,1);
SW_out_ac_1h = nan(47856,1);
LW_in_ac_1h = nan(47856,1);
%LW_out_ac_1h = nan(47856,1);
for l=68:6:length(time_for)-94
    i=i+1;
    time_1h(i) = time_for(l);
    T_air_ac_1h(i) = mean(T_air_for(l-5:l));
    RH_ac_1h(i) = mean(RH_for(l-5:l));
    P_air_ac_1h(i) = mean(P_air_for(l-5:l));
    Prec_ac_1h(i) = sum(Prec_for(l-5:l));   % total amount during timestep?
    Wind_ac_1h(i) = mean(Wind_speed_for(l-5:l));
    SW_in_bc_1h(i) = mean(SW_in_bc(l-5:l));
    SW_out_bc_1h(i) = mean(SW_out_bc(l-5:l));
    LW_in_bc_1h(i) = mean(LW_in_bc(l-5:l));
    LW_out_bc_1h(i) = mean(LW_out_bc(l-5:l));
    SW_in_ac_1h(i) = mean(SW_in_ac(l-5:l));
    SW_out_ac_1h(i) = mean(SW_out_ac(l-5:l));
    LW_in_ac_1h(i) = mean(LW_in_ac(l-5:l));
    %LW_out_ac_1h(i) = mean(LW_out_ac(l-5:l));
end
SW_in_ac_1h(22821:22822) = 0;   % weirdly -10^{-14}

%-----------------------------  ground data  -----------------------------%
% soil moisture
%{
Soil water content [%] from FluxNet measurements, which provides 3 values
but only one for the period 17 December 2006 - 31 December 2014 (the other
two stop in January 2010). Values are provided half-hourly.
%}
soil_water_content ...
    = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,32,[1,32,315552,34]);
% cut down to 13 October 2007 1:00 - 29 March 2013 0:00
SWC = soil_water_content(188978:284689,1);

% pick hourly values
i=0;
SWC_1h = nan(size(LW_in_bc_1h));
for l=2:2:length(SWC)
    i=i+1;
    SWC_1h(i) = mean(SWC(l-1:l));
end

% soil temperature
%{
Soil temperature (???) from FluxNet measurements, which provides 3 values
but only one for the entirity of period 17 December 2006 - 31 December 2014
(the other two stop in January 2010). Values are provided half-hourly.
%}
TS = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,26,[1,26,315552,28]);
% cut down to 13 October 2007 1:00 - 29 March 2013 0:00
TS_cut = TS(188978:284689,:);
% pick hourly values
i=0;
TS_cut_1h = nan(size(LW_in_bc_1h));
for l=2:2:length(TS_cut)
    i=i+1;
    for j=1:3
        TS_cut_1h(i,j) = mean(TS_cut(l-1:l,j));
    end
end
%{
Values are almost entirely >0°C for the continuous time series (except for
January/February 2010). The other two are always warmer indicating a
temperature profile, so we assume unfrozen ground continuously as we don't
know about measurement depth (maybe on FluxNet website?).
%}

% snow depth
%{
Snow depth based on two different measurements.
1) Occasional, irregular manual snow course measurements for open and
   forested site, available for seasons 2007/08 and 2008/09.
2) Continuous automatic measurements of daily resolution from SLF snow
   station (SR50). This is an open site less than 1km away from the forest.
We create a scaling ratio from manual measurements, which is then applied
to daily automatic measurements from snow station to get continuous values
for the entire simulation period. Daily values are then applied to every
hour of the day.
%}
% manual measurements
    % 2007/08
bla = xlsread('Snow_Data_Seehornwald_2008.xlsx');
SD_open_2008 = bla(:,4:14);     % 10 locations & average
SD_forest_2008 = bla(:,15:25);  % 10 locations & average
clear bla
    % 2008/09
bla = xlsread('Snow_Data_Seehornwald_2009.xlsx');
SD_open_2009 = bla(:,4:14);     % 10 locations & average
SD_forest_2009 = bla(:,15:35);  % 20 locations & average
clear bla
%{
Less measurements for 2007/08 than 2008/09 (9 field trips vs 35 field trips).
Using linear regression to determine scaling ratio yields similar slopes
(0.58... to 0.57...) when using both years or onyl 2008/09, but the residual
is larger when including 2007/08, so we only use 2008/09. This gives
y = 0.5739 * x - 0.04934, and 0.5739 is used to scale continuous measurements.
%}
scaling = 0.5739;
bla = importdata('HS_Data_Davos_1998-2016.txt','\t',15);
z_snow_open = bla.data(:,7);
clear bla
    % convert from cm to m
z_snow_open = z_snow_open/100;
    % cut down to 13 October 2007 - 28 March 2013
z_snow_open_cut = z_snow_open(3329:5322);
    % spread out daily values (6:00) as hourly values (1:00 to 0:00)
z_snow_open_cut_1h = nan(length(z_snow_open_cut)*24,1);
for l=1:length(z_snow_open_cut)
    z_snow_open_cut_1h((l-1)*24+1:l*24) = z_snow_open_cut(l);
end
    % scale to forest approximations
z_snow_forest_1h = scaling * z_snow_open_cut_1h;

% surface albedo
alb_gr_1h = nan(size(time_1h));
for t=1:length(time_1h)
    if SW_out_bc_1h(t) > 1 && SW_in_bc_1h(t) > 1
        alb_gr_1h(t) = SW_out_bc_1h(t)/SW_in_bc_1h(t);
        alb_gr_1h(t) = min(alb_gr_1h(t),1);
    end
end

% snow fraction
%{
Using surface albedo measurements (filtered for noon) there's seems to be a
threshold of snow depth for which surface albedo indicates at least patches
of soil are visible. This threshold is ca. 0.25m. Therefore, we apply snow
fraction values of 1 for snow depths of >=0.25, 0.5 (yes, this is crude)
for snow depths of <0.25 but >0 and 0 for snow depths of 0. For 2008/09,
there is a clear coinciding linear decrease of surface albedo and snow
depth so that we use a linear decrease for snow cover fraction as well.
%}
frac_snow_1h = zeros(size(z_snow_forest_1h));
for l=1:length(frac_snow_1h)
    if z_snow_forest_1h(l) >= 0.25
        frac_snow_1h(l) = 1;
    elseif z_snow_forest_1h(l) < 0.25 && z_snow_forest_1h(l) >0
        frac_snow_1h(l) = 0.5;
    end
end
    % correction based on surface albedo measurements and snow depth
frac_snow_1h(21313:21504) = 0.5;
frac_snow_1h(21505:21817) = 1;
frac_snow_1h(38929:39649) = 0.5;
for l=12865:13464
    frac_snow_1h(l) = frac_snow_1h(l-1) - 1/601;
end

% surface temperature
T_surf_1h = nan(size(time_1h));
for t=1:length(time_1h)
    em_gr = frac_snow_1h(t)*0.97 + (1-frac_snow_1h(t))*0.96;
    T_surf_1h(t) = nthroot(LW_out_bc_1h(t)/(5.67*10^(-8)*em_gr),4);
end
    % linear interpolation of sub-daily gaps
T_surf_1h_filled = T_surf_1h;
for t=2:length(T_surf_1h)
    if isnan(T_surf_1h(t)) && ~isnan(T_surf_1h(t-1))
        gap_start = t;
        gap = find(~isnan(T_surf_1h(t:end)),1)-1;
        gap_end = t+gap-1;  % -1 because gap includes t
        if gap < 24
            diff = (T_surf_1h(gap_end+1)-T_surf_1h(gap_start-1))/(gap+2);
            for g=gap_start:gap_end
                T_surf_1h_filled(g) = T_surf_1h_filled(g-1) + diff;
            end
        end
    end
end


%--------------------------------------------------------------------------
%-----------------------------  gap filling  ------------------------------
%--------------------------------------------------------------------------
% using FluxNet tower measurements (which strangely differ)
T_air_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,2,[1,2,315552,2]);
T_air_FluxNet = T_air_FluxNet(188978:284689,1);
i=0;
T_air_FluxNet_1h = nan(size(T_air_ac_1h));
for l=1:2:length(T_air_FluxNet)-1 % comp by eye
    i=i+1;
    T_air_FluxNet_1h(i) = T_air_FluxNet(l);
end

LW_in_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,7,[1,7,315552,7]);
LW_in_FluxNet = LW_in_FluxNet(188978:284689,1);
i=1; % comp by eye
LW_in_FluxNet_1h = nan(size(LW_in_ac_1h));
for l=1:2:length(LW_in_FluxNet)-3
    i=i+1;
    LW_in_FluxNet_1h(i) = LW_in_FluxNet(l);
end

SW_in_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,5,[1,5,315552,5]);
SW_in_FluxNet = SW_in_FluxNet(188978:284689,1);
i=0;
SW_in_FluxNet_1h = nan(size(SW_in_ac_1h));
for l=1:2:length(SW_in_FluxNet)-1 % comp by eye and mean diff
    i=i+1;
    SW_in_FluxNet_1h(i) = SW_in_FluxNet(l);
end

SW_in_pot_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,4,[1,4,315552,4]);
SW_in_pot_FluxNet = SW_in_pot_FluxNet(188978:284689,1);
i=0;
SW_in_pot_FluxNet_1h = nan(size(SW_in_ac_1h));
for l=2:2:length(SW_in_pot_FluxNet) % more logical to avoid 0
    i=i+1;
    SW_in_pot_FluxNet_1h(i) = mean(SW_in_pot_FluxNet(l-1:l));
end

SW_out_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
        1,22,[1,22,315552,22]);
SW_out_FluxNet = SW_out_FluxNet(188978:284689,1);
for l=1:length(SW_out_FluxNet)
    if SW_out_FluxNet(l) == -9999
        SW_out_FluxNet(l) = nan;
    end
end
i=0;
SW_out_FluxNet_1h = nan(size(SW_out_ac_1h));
for l=2:2:length(SW_out_FluxNet)
    i=i+1;
    SW_out_FluxNet_1h(i) = mean(SW_out_FluxNet(l-1:l));
end

Wind_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
        1,15,[1,15,315552,15]);
Wind_FluxNet = Wind_FluxNet(188978:284689,1);
i=0;
Wind_FluxNet_1h = nan(size(Wind_ac_1h));
for l=1:2:length(Wind_FluxNet)-1 % comp by eye
    i=i+1;
    Wind_FluxNet_1h(i) = Wind_FluxNet(l);
end

Pair_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
        1,11,[1,11,315552,11]);
Pair_FluxNet = Pair_FluxNet(188978:284689,1);
i=1; % comp by eye
Pair_FluxNet_1h = nan(size(P_air_ac_1h));
for l=3:2:length(Pair_FluxNet)-1 % comp by eye
    i=i+1;
    Pair_FluxNet_1h(i) = mean(Pair_FluxNet(l-1:l));
end
Pair_FluxNet_1h = Pair_FluxNet_1h*10;

RH_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
        1,19,[1,19,315552,19]);
RH_FluxNet = RH_FluxNet(188978:284689,1);
for l=1:length(RH_FluxNet)
    if RH_FluxNet(l) == -9999
        RH_FluxNet(l) = nan;
    end
end
i=0;
RH_FluxNet_1h = nan(size(RH_ac_1h));
for l=1:2:length(RH_FluxNet)-1 % comp by eye
    i=i+1;
    RH_FluxNet_1h(i) = RH_FluxNet(l);
end

Prec_FluxNet = csvread('FLX_CH-Dav_FLUXNET2015_SUBSET_HH_1997-2014_1-3.csv',...
    1,13,[1,13,315552,13]);
Prec_FluxNet = Prec_FluxNet(188978:284689,1);
i=0; % comp by eye - but either lag at start or lead later on, weird change
Prec_FluxNet_1h = nan(size(Prec_ac_1h));
for l=2:2:length(Prec_FluxNet)
    i=i+1;
    Prec_FluxNet_1h(i) = sum(Prec_FluxNet(l-1:l));
end

% fill gaps by interpolation from FluxNet
T_air_ac_1h_filled = GapFilling(T_air_ac_1h,T_air_FluxNet_1h);
LW_in_ac_1h_filled = GapFilling(LW_in_ac_1h,LW_in_FluxNet_1h);
SW_in_ac_1h_filled = GapFilling(SW_in_ac_1h,SW_in_FluxNet_1h);
P_air_ac_1h_filled = GapFilling(P_air_ac_1h,Pair_FluxNet_1h);
Wind_ac_1h_filled = GapFilling(Wind_ac_1h,Wind_FluxNet_1h);
for l=1:length(Wind_ac_1h_filled)
    Wind_ac_1h_filled(l) = max(Wind_ac_1h_filled(l),0);
end
RH_ac_1h_filled = GapFilling(RH_ac_1h,RH_FluxNet_1h);
for l=1:length(RH_ac_1h_filled)
    RH_ac_1h_filled(l) = max(min(RH_ac_1h_filled(l),100),0);
end
Prec_ac_1h_filled = GapFilling(Prec_ac_1h,Prec_FluxNet_1h);
for l=1:length(Prec_ac_1h_filled)
    Prec_ac_1h_filled(l) = max(Prec_ac_1h_filled(l),0);
end
 

%--------------------------------------------------------------------------
%----------------------------  miscellaneous  ----------------------------%
%--------------------------------------------------------------------------
% effective emissivity of the sky
ema_eff_1h = nan(size(T_air_ac_1h_filled));
for t=1:length(time_1h)
    ema_eff_1h(t) = LW_in_ac_1h_filled(t)/(5.67*10^(-8)*(T_air_ac_1h_filled(t)+273.15)^4);
end
% ratio of insolation
r_SW_FluxNet_1h = nan(size(SW_in_FluxNet_1h));
for l=1:length(SW_in_FluxNet_1h)
    if SW_in_pot_FluxNet_1h(l) >= 1
        r_SW_FluxNet_1h(l) = SW_in_FluxNet_1h(l)/SW_in_pot_FluxNet_1h(l);
    end
    r_SW_FluxNet_1h(l) = min(r_SW_FluxNet_1h(l),1);
end

%{
Weirdly, time seems to be shifted by 1h. Although not symmetrical,
insolation peaks at 13:00, starts at 6:00 and ends at 20:00. Therefore, we
shift time_1h backwards by 1h, and now it also corresponds to solar angles
calculated by the models.
%}
time_1h = time_1h-1/24;
save('ForcingData_ToyModel_Seehornwald.mat','time_1h',...
    'T_air_ac_1h_filled','RH_ac_1h_filled','P_air_ac_1h_filled',...
    'Prec_ac_1h_filled','Wind_ac_1h_filled','SW_in_ac_1h_filled',...
    'LW_in_ac_1h_filled','LW_in_bc_1h','ema_eff_1h','r_SW_FluxNet_1h',...
    'LW_out_bc_1h','T_surf_1h_filled','SW_out_bc_1h','SW_in_bc_1h',...
    'alb_gr_1h','z_snow_forest_1h','frac_snow_1h','SWC_1h')
toc