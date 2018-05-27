tic
clear
close all

%-------------------------------------------------------------------------%
%-----------------------------  import data  -----------------------------%
%-------------------------------------------------------------------------%

%---------------------  meteorological forcing data  ---------------------%
%{
"Gaps were filled using CRUNCEP V8."
More extensive data can be found in Borden_Met_2013.csv
%}
bla = importdata('Borden_2013_Forcing.prn');
hours = bla(:,1);
mins = bla(:,2);
DoY = bla(:,3);
Year = bla(:,4);
time_met = nan(size(DoY));
for l=1:length(DoY)
    [day,month] = DoY_to_DayNMonth(Year(l),DoY(l));
    time_met(l) = datenum(Year(l),month,day,hours(l),mins(l),0);
end
SW_in_ac = bla(:,5);    % [W m^{-2}]
LW_in_ac = bla(:,6);    % [W m^{-2}]
Precip = bla(:,7);      % [kg m^{-2} s^{-1}] == [mm s^{-1}]
T_air_ac = bla(:,8);    % [°C]
SpecHum = bla(:,9);     % [kg kg^{-1}]
Wind = bla(:,10);       % [m s^{-1}]
Pair = bla(:,11);       % [Pa]
clear bla

% convert specific humidity into relative humidity using CLM calculations
RelHum = nan(size(SpecHum));
for l=1:length(SpecHum)
    tdc = min(50,max(-50,T_air_ac(l)));
    if tdc > 0      % saturation vapor pressure over water [Pa]
        e = 100*(6.107799961 + tdc*(4.436518521*10^(-1) + tdc*(1.428945805*10^(-2)...
            + tdc*(2.650648471*10^(-4) + tdc*(3.031240396*10^(-6)...
            + tdc*(2.034080948*10^(-8) + tdc*6.136820929*10^(-11)))))));
    else                    % saturation vapor pressure over ice [Pa]
        e = 100*(6.109177956 + tdc*(5.034698970*10^(-1) + tdc*(1.886013408*10^(-2)...
            + tdc*(4.176223716*10^(-4) + tdc*(5.824720280*10^(-6)...
            + tdc*(4.838803174*10^(-8) + tdc*1.838826904*10^(-10)))))));
    end
    qsat = 0.622*e/(Pair(l)-0.378*e);
    RelHum(l) = SpecHum(l)/qsat * 100;
%{
A handful of values are slightly above 100%, but the last 11 values of the
year jump to around 120%, likely due to a sudden drop in air pressure.
Since the end of the year is not considered in our analysis, we just limit
relative humidity of these timesteps to 100%, otherwise they would be
emitted.
%}
    RelHum(l) = max(0,min(RelHum(l),100));
end

%{
Lead/lag between forcing file and csv files of 10 timesteps (5 hours),
exactly those that were unrealistic for relative humidity. Also, this lag
was seen when comparing to maple and aspen temperatures (down below).
However, timestamps do not feature this so that variables have to be used
from timestep 11 onwards but time used from timestep 1 onwards.
%}
% import data pre-gap filling
bla = importdata('BordenMet_2013.csv');
JulDay_pure = bla.data(:,2);
T_air_ac_33_pure = bla.data(:,5);
T_air_ac_41_pure = bla.data(:,7);
RelHum_33_pure = bla.data(:,6);
RelHum_41_pure = bla.data(:,8);
Wind_43_pure = bla.data(:,20);
SW_in_ac_34_pure = bla.data(:,14);
SW_in_ac_34_pure = NaN999(SW_in_ac_34_pure);
SW_in_ac_34_pure = NaN9999(SW_in_ac_34_pure);
for l=1:length(SW_in_ac_34_pure)
    if SW_in_ac_34_pure(l) < 0
        SW_in_ac_34_pure(l) = 0;
    end
end
SW_out_ac_34_pure = bla.data(:,15);
LW_in_ac_34_pure = bla.data(:,16);
LW_out_ac_34_pure = bla.data(:,17);
clear bla
T_air_ac_33_pure = NaN9999(T_air_ac_33_pure);
    T_air_ac_33_pure = T_air_ac_33_pure + 273.15;
RelHum_33_pure = NaN9999(RelHum_33_pure);
Wind_43_pure = NaN9999(Wind_43_pure);
SW_in_ac_34_pure = NaN9999(SW_in_ac_34_pure);
SW_out_ac_34_pure = NaN9999(SW_out_ac_34_pure);
LW_in_ac_34_pure = NaN9999(LW_in_ac_34_pure);
LW_out_ac_34_pure = NaN9999(LW_out_ac_34_pure);

%------------------------------  soil data  ------------------------------%
bla = importdata('BordenSoil_2013.csv');
JulDay = bla.data(:,2);
SoilTemp = bla.data(:,5:16);    % [°C]
SoilMoist = bla.data(:,27:end); % [m^3 m^{-3}]
%{
Soil temperature measurements:
1  - Site 1, 5cm
2  - Site 1, 5cm
3  - Site 2, 5cm
4  - Site 2, 5cm
5  - Site 2, 10cm
6  - Site 2, 20cm
7  - Site 2, 50cm
8  - Site 2, 100cm
9  - Site 1, 10cm
10 - Site 1, 20cm
11 - Site 1, 50cm
12 - Site 2, 100cm      ... Site 1?
All measurements feature gaps in unison, although more than meteorological
forcing. Some sensors might be drifting as 10cm and 20cm temperatures for
site 2 "seem to warm and don't fit the expected profile very well". The
same can be said about measurement #2 at site 1.
Therefore, we use site 1 excluding measurement #2 to determine fraction of
frozen soil.
Soil moisture measurements:
1  - Site A, 2cm
2  - Site A, 5cm
3  - Site A, 10cm       ... lots of gaps
4  - Site A, 20cm
5  - Site A, 50cm
6  - Site A, 100cm
7  - Site B, 50cm
8  - Site B, 2cm
9  - Site B, 10cm       ... almost completely empty
10 - Site B, 20cm
11 - Site B, 50cm       ... lots of gaps
12 - Site B, 100cm      ... completely empty
We use site A for obvious reasons and exclude layer 10cm when featuring
gaps.
%}
SoilTemp_filter = horzcat(SoilTemp(:,1),SoilTemp(:,9:12));
SoilMoist_filter = SoilMoist(:,1:6);
T_Maple = bla.data(:,18:20);    % 2m, 5m and 10m
T_Aspen = bla.data(:,21:23);    % 2m, 5m and 9m
clear bla
for i=1:length(SoilTemp_filter(1,:))
    SoilTemp_filter(:,i) = NaN9999(SoilTemp_filter(:,i));
end
    SoilTemp_filter = SoilTemp_filter + 273.15;
T_topsoil = SoilTemp_filter(:,1);
for i=1:length(SoilMoist_filter(1,:))
    SoilMoist_filter(:,i) = NaN9999(SoilMoist_filter(:,i));
end
for i=1:3
    T_Maple(:,i) = NaN9999(T_Maple(:,i));
    T_Aspen(:,i) = NaN9999(T_Aspen(:,i));
end
    T_Maple = T_Maple + 273.15;
    T_Aspen = T_Aspen + 273.15;
    
%-----------------------------  forest data  -----------------------------%
%{
There are two pyranometers above the canopy, one Eppley and one Kipp&Zonen.
However, the Eppley one had not been calibrated in years so that the
Kipp&Zonen one should be used and indeed, the measurements are the same as
the ones in BordenMet_2013.csv (from timestep 11 onwards).
There are several pyranometers beneath the canopy, which had been calibrated
before being deployed but not since. Also, a Kipp&Zonen CNR1 had been
deployed both for measuring incoming and outgoing sub-canopy radiation, but
since the technician was replaced they had not been checked for frost or
snow, although checked that they were level.
%}
bla = xlsread('Preliminary Transmissivity Borden_2013.xlsx');
% SW_in_ac_forest = bla(2:end,5);
    % Kipp&Zonen sub-canopy radiometers
SW_in_bc_forest = bla(2:end,10);
SW_out_bc_forest = bla(2:end,11);
LW_in_bc_forest = bla(2:end,12);
LW_out_bc_forest = bla(2:end,13);
    % array of different pyranometers
SW_in_bc_forest_array = bla(2:end,14:25);
SW_in_bc_forest_all = horzcat(SW_in_bc_forest,SW_in_bc_forest_array);
SW_in_bc_forest_array_avg = nan(size(SW_in_bc_forest));
SW_in_bc_forest_all_avg = nan(size(SW_in_bc_forest));
for l=1:length(SW_in_bc_forest)
    for i=1:length(SW_in_bc_forest_array(1,:))
        if SW_in_bc_forest_array(l,i) < 0
            SW_in_bc_forest_array(l,i) = nan;
        end
    end
    for i=1:length(SW_in_bc_forest_all(1,:))
        if SW_in_bc_forest_all(l,i) < 0
            SW_in_bc_forest_all(l,i) = nan;
        end
    end
    SW_in_bc_forest_array_avg(l) = nanmean(SW_in_bc_forest_array(l,:));
    SW_in_bc_forest_all_avg(l) = nanmean(SW_in_bc_forest_all(l,:));
end

% calculate surface albedos
alb_KnZ = nan(size(SW_in_bc_forest));
alb_array = nan(size(SW_in_bc_forest));
alb_all = nan(size(SW_in_bc_forest));
for l=1:length(SW_in_bc_forest)
    if SW_in_bc_forest(l) >= 1
        alb_KnZ(l) = SW_out_bc_forest(l)/SW_in_bc_forest(l);
    end
    if SW_in_bc_forest_array_avg(l) >= 1
        alb_array(l) = SW_out_bc_forest(l)/SW_in_bc_forest_array_avg(l);
    end
    if SW_in_bc_forest_all_avg(l) >= 1
        alb_all(l) = SW_out_bc_forest(l)/SW_in_bc_forest_array_avg(l);
    end
end
%{
The Kipp&Zonen pyranometer yielded consistently higher insolation than the
other pyranometers, which results in lower albedo values. As there is only
one pyranometer for reflected sub-canopy shortwave radiation we use the K&Z
pyranometer for incoming shortwave radiation as well.
%}
clear bla


%-------------------------------------------------------------------------%
%---------------------------  Post-Processing  ---------------------------%
%-------------------------------------------------------------------------%

%----------------  check sub-canopy LWR for snow & frost  ----------------%
%{
Suspicious values of constant LWR around 314.3 W m^{-2} set to nan.
%}
% LW_in_bc_forest(274:277) = nan;
LW_in_bc_forest(1325:1374) = nan;
LW_in_bc_forest(1977:1981) = nan;
LW_in_bc_forest(2149:2160) = nan;
% LW_in_bc_forest(2368:2371) = nan;
LW_in_bc_forest(2552:2568) = nan;
% LW_in_bc_forest(2615:2617) = nan;
LW_in_bc_forest(2733:2752) = nan;
LW_in_bc_forest(2778:2803) = nan;
LW_in_bc_forest(3177:3180) = nan;
LW_in_bc_forest(4864:4880) = nan;

%-----------------  unphysical values for depth of 10cm  -----------------%
SoilMoist_filter(5175,3) = nan;
SoilMoist_filter(6176,3) = nan;
SoilMoist_filter(7903,3) = nan;
SoilMoist_filter(8309,3) = nan;
SoilMoist_filter(10318,3) = nan;
SoilMoist_filter(13754:13755,3) = nan;
SoilMoist_filter(15521,3) = nan;

%-----------------------  average soil vertically  -----------------------%
% soil moisture
%{
Depths: 2cm, 5cm, 10cm, 20cm, 50cm, 100cm.
One level includes nans when others do not, so we have to calculate
non-nan depth for every time step.
%}
levels = [0.035 0.04 0.075 0.2 0.4 0.5];
temp1 = nan(size(SoilMoist_filter));
temp2 = nan(size(SoilMoist_filter));
SoilMoist_filter_avg = nan(size(SoilMoist_filter(:,1)));
for l=1:length(SoilMoist_filter_avg)
    for i=1:length(levels)
        temp1(l,i) = SoilMoist_filter(l,i)*levels(i);
        if ~isnan(temp1(l,i))
            temp2(l,i) = levels(i);
        end
    end
    SoilMoist_filter_avg(l) = nansum(temp1(l,:))/nansum(temp2(l,:));
end
clear temp1 temp2

% fraction of frozen soil - Depths: 5cm, 10cm, 20cm, 50cm, 100cm.
boundaries = [0.05 0.1 0.2 0.5 1];
levels = boundaries(2:end)-boundaries(1:end-1);
temp = nan(length(T_topsoil),length(levels));
T_freez = 273.15;
frac_soil_frozen = nan(size(T_topsoil));
for l=1:length(T_topsoil)
    for i=1:length(levels)
        if SoilTemp_filter(l,i) <= T_freez && SoilTemp_filter(l,i+1) <= T_freez
            temp(l,i) = levels(i);
        elseif SoilTemp_filter(l,i) > T_freez && SoilTemp_filter(l,i+1) > T_freez
            temp(l,i) = 0;
        else
            mintemp = min(SoilTemp_filter(l,i),SoilTemp_filter(l,i+1));
            maxtemp = max(SoilTemp_filter(l,i),SoilTemp_filter(l,i+1));
            temp(l,i) = (T_freez-mintemp)/(maxtemp-mintemp) * levels(i);
        end
    end
    frac_soil_frozen(l) = sum(temp(l,:))/sum(levels);
end
clear temp mintemp maxtemp

%----------------------  calculate hourly averages  ----------------------%
i=0;
for l=12:2:length(JulDay_pure)
    i=i+1;
    JulDay_1h(i) = JulDay_pure(l-10);
    time_1h(i) = time_met(l-9); % starts at 0:00 not 1:00
    Wind_1h(i) = mean(Wind_43_pure(l-1:l));
    RH_1h(i) = mean(RelHum_33_pure(l-1:l));
    T_air_ac_1h(i) = mean(T_air_ac_33_pure(l-1:l));
    LW_in_ac_1h(i) = mean(LW_in_ac_34_pure(l-1:l));
    SW_in_ac_1h(i) = mean(SW_in_ac_34_pure(l-1:l));
    SW_in_bc_1h(i) = mean(SW_in_bc_forest(l-1-10:l-10));
    SW_in_bc_all_1h(i) = mean(SW_in_bc_forest_all_avg(l-1-10:l-10));
    SW_out_bc_1h(i) = mean(SW_out_bc_forest(l-1-10:l-10));
    alb_surf_1h(i) = mean(alb_KnZ(l-1-10:l-10));
    alb_surf_all_1h(i) = mean(alb_all(l-1-10:l-10));
    LW_in_bc_1h(i) = mean(LW_in_bc_forest(l-1-10:l-10));
    LW_out_bc_1h(i) = mean(LW_out_bc_forest(l-1-10:l-10));
    frac_soil_frozen_1h(i) = mean(frac_soil_frozen(l-1:l));
    SoilMoist_filter_avg_1h(i) = mean(SoilMoist_filter_avg(l-1:l));
    Precip_1h(i) = mean(Precip(l-1-10:l-10));   % already rate !!!!
    for j=1:3
        T_Maple_1h(i,j) = mean(T_Maple(l-1:l,j));
        T_Aspen_1h(i,j) = mean(T_Aspen(l-1:l,j));
    end
end

%----------------------------  interpolation  ----------------------------%
%{
As was done by Essery et al. (2016) for data at Sodankylä, gaps of 4 hours
or shorter are filled by linear interpolation. At least for shortwave
radiation an approach based on diurnal cycles from previous days could be
applied.
%}
Wind_1h_filled = LinearInterp_4h(Wind_1h);
RH_1h_filled = LinearInterp_4h(RH_1h);
T_air_ac_1h_filled = LinearInterp_4h(T_air_ac_1h);
LW_in_ac_1h_filled = LinearInterp_4h(LW_in_ac_1h);
SW_in_ac_1h_filled = LinearInterp_4h(SW_in_ac_1h);
SW_in_bc_1h_filled = LinearInterp_4h(SW_in_bc_1h);
SW_in_bc_all_1h_filled = LinearInterp_4h(SW_in_bc_all_1h);
SW_out_bc_1h_filled = LinearInterp_4h(SW_out_bc_1h);
LW_out_bc_1h_filled = LinearInterp_4h(LW_out_bc_1h);
frac_soil_frozen_1h_filled = LinearInterp_4h(frac_soil_frozen_1h);
SoilMoist_filter_avg_1h_filled = LinearInterp_4h(SoilMoist_filter_avg_1h);
%{
Wind measurements features two longer gaps (8 and 6 hours) that other met
variables don't have. Both occur during periods of low wind speeds and we
interpolate linearly.
%}
for t=1394:1401
    Wind_1h_filled(t) = Wind_1h_filled(t-1)...
        + (Wind_1h_filled(1402)-Wind_1h_filled(1393))/9;
end
for t=2433:2438
    Wind_1h_filled(t) = Wind_1h_filled(t-1)...
        + (Wind_1h_filled(2439)-Wind_1h_filled(2432))/7;
end

% also interpolate longer gaps (7 - 8 hours) that only exist for soil
for t=1232:1238
    frac_soil_frozen_1h_filled(t) = frac_soil_frozen_1h_filled(t-1)...
        + (frac_soil_frozen_1h_filled(1239)-frac_soil_frozen_1h_filled(1231))/8;
    SoilMoist_filter_avg_1h_filled(t) = SoilMoist_filter_avg_1h_filled(t-1)...
        + (SoilMoist_filter_avg_1h_filled(1239)-SoilMoist_filter_avg_1h_filled(1231))/8;
end
for t=1269:1275
    frac_soil_frozen_1h_filled(t) = frac_soil_frozen_1h_filled(t-1)...
        + (frac_soil_frozen_1h_filled(1276)-frac_soil_frozen_1h_filled(1268))/8;
    SoilMoist_filter_avg_1h_filled(t) = SoilMoist_filter_avg_1h_filled(t-1)...
        + (SoilMoist_filter_avg_1h_filled(1277)-SoilMoist_filter_avg_1h_filled(1268))/9;
end
SoilMoist_filter_avg_1h_filled(1276) = SoilMoist_filter_avg_1h_filled(1275)...
    + (SoilMoist_filter_avg_1h_filled(1277)-SoilMoist_filter_avg_1h_filled(1268))/9;

%{
There are still two gaps of several hours (1 February 0:00 to 9:00, 25 April
0:00 to 8:00) that all meteorological variables feature, which we exclude
from simulations. Simulation stops before these gaps and continues
afterwards with variables required from previous timestep taken from before
the gap. Since evaluation will start 1:00 there is enough time for spin-up.
%}
gap_1 = 744:753;
gap_2 = 2736:2744;

save('ForcingData_ToyModel_Borden.mat','Precip_1h','RH_1h_filled',...
    'Wind_1h_filled','T_air_ac_1h_filled','LW_in_ac_1h_filled',...
    'SW_in_ac_1h_filled','SW_in_bc_1h_filled','SW_in_bc_all_1h_filled',...
    'alb_surf_1h','alb_surf_all_1h','SW_out_bc_1h_filled','LW_out_bc_1h_filled',...
    'LW_in_bc_1h','T_Maple_1h','T_Aspen_1h','frac_soil_frozen_1h_filled',...
    'SoilMoist_filter_avg_1h_filled','JulDay_1h','time_1h','gap_1','gap_2')
toc