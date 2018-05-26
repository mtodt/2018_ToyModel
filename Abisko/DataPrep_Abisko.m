tic
clear
close all

%-------------------------------------------------------------------------%
%-----------------------------  forest data  -----------------------------%
ab_data = xlsread('AllAbiskoData 190716.xls');
%{
Cut out period of available data. Start limited by sub-canopy air
temperature - from 10 March 13:35. End limited by radiation measurements in
open - to 5 April 8:35.
%}

% open = ground-level measurements
LW_in_open = ab_data(846:8274,91);
    % from Nick's script: "BF3 data already corrected for Abisko"
SW_in_open = ab_data(846:8274,89);
SW_in_diff = ab_data(846:8274,90);
f_SW_diff = ones(size(SW_in_diff));
for l=1:length(f_SW_diff)
    if SW_in_open(l) > 0
        f_SW_diff(l) = min(1,SW_in_diff(l)/SW_in_open(l));
    end
end

% continuous site
SW_in_bc_C = ab_data(846:8274,73:82);
LW_in_bc_C = ab_data(846:8274,83:86);
Rad_in_bc_C_flag = ab_data(846:8274,87);
for l=1:length(Rad_in_bc_C_flag)
    if Rad_in_bc_C_flag(l) == 3
        LW_in_bc_C(l,:) = nan;
        SW_in_bc_C(l,:) = nan;
    end
end
LW_in_bc_C_avg = nan(size(squeeze(LW_in_bc_C(:,1))));
SW_in_bc_C_avg = nan(size(squeeze(SW_in_bc_C(:,1))));
for l=1:length(Rad_in_bc_C_flag)
    LW_in_bc_C_avg(l) = mean(LW_in_bc_C(l,:));
    SW_in_bc_C_avg(l) = mean(SW_in_bc_C(l,:));
end
T_air_sub_C = ab_data(846:8274,88)+273.15;
T_veg_TC_C = ab_data(846:8274,8:71)+273.15;
%{
Based on descriptions in Nick's script, some instruments are dismissed for
several reasons:
- unreliable,
- imbalance N-S (or E-W, if available), and
- dead wood.
Lichen on birch might not b a sensible comparison as it should not be
considered in CLM4.5 for albedo, etc. (?), but instruments were imbalanced
on that tree anyway.
Some instruments were moved on 31 March/1 April, however, these do not
affect selection, as all of them are dismissed before and after.
%}
T_veg_TC_C_QC = horzcat(T_veg_TC_C(:,1:32),T_veg_TC_C(:,41:56));
T_veg_TC_C_flag = ab_data(846:8274,72);
T_veg_TC_C_avg = nan(size(T_veg_TC_C_flag));
for l=1:length(T_veg_TC_C_flag)
    if T_veg_TC_C_flag(l) == 3
        T_veg_TC_C_QC(l,:) = nan;
    end
    T_veg_TC_C_avg(l) = mean(T_veg_TC_C_QC(l,:));
end

% hourly averages
%{
Translating measurements from temporal resolution of 5 minutes into hourly
values. Coverages starts at 13:35 on 10 March and ends at 8:35 on 5 April,
with hourly values set to the end of the averaging period, i.e. 15:00 for
the first time step, which starts at 14:05. Evaluation period should then
be set to 1:00 on 11 March until 0:00 on 5 April.
%}
i=0;
for t=6+12:12:length(LW_in_open)-7
    i=i+1;
    SW_in_open_1h(i) = mean(SW_in_open(t-11:t));
    f_SW_diff_1h(i) = mean(f_SW_diff(t-11:t));
    SW_in_bc_C_avg_1h(i) = mean(SW_in_bc_C_avg(t-11:t));
    LW_in_open_1h(i) = mean(LW_in_open(t-11:t));
    LW_in_bc_C_avg_1h(i) = mean(LW_in_bc_C_avg(t-11:t));
    for n=1:4
        LW_in_bc_C_1h(i,n) = mean(LW_in_bc_C(t-11:t,n));
    end
    T_air_sub_C_1h(i) = mean(T_air_sub_C(t-11:t));
    T_veg_TC_C_avg_1h(i) = mean(T_veg_TC_C_avg(t-11:t));
    time_1h(i) = datenum(ab_data(845+t,1),ab_data(845+t,2),ab_data(845+t,3),...
        ab_data(845+t,5),ab_data(845+t,6),ab_data(845+t,7));
end
%{
Open area radiation measurements feature a gap of about 50 min. Since this
occurs in late morning a linear interpolation of hourly values is
applied. These timesteps should then be excluded from evaluation.
%}
for t=93:94
    LW_in_open_1h(t) = LW_in_open_1h(t-1) + (LW_in_open_1h(95)-LW_in_open_1h(92))/3;
    SW_in_open_1h(t) = SW_in_open_1h(t-1) + (SW_in_open_1h(95)-SW_in_open_1h(92))/3;
end


%-------------------------------------------------------------------------%
%-----------------------------  ground data  -----------------------------%
% snow surface temperature from sub-canopy air temperature capped at 0C
T_surf_C = nan(size(T_air_sub_C_1h));
for t=1:length(time_1h)
    if ~isnan(T_air_sub_C_1h(t))
        T_surf_C(t) = min(T_air_sub_C_1h(t),273.15);
    end
end

%{
from Annika Kristoffersson:
"You will also find precipitation and snow depth data from our manual
readings. The precipitation from the AWS is not reliable, specially during winter."
Manual measurements daily not hourly, so we use precip from the AWS for
timing throughout a day but use manual measurements to determine the
magnitude.
%}
bla = xlsread('SnowDepth.xlsx');
z_snow_manual = bla(10:36,2);
z_snow_manual = z_snow_manual/100;  % cm -> m
clear bla
% set daily snow depth for every hour of that day
z_snow = nan(size(T_surf_C));
z_snow(1:10) = z_snow_manual(1);
z_snow(end-7:end) = z_snow_manual(end);
i=1;    % 10 March already set
for t=10+24:24:length(z_snow)-8
    i=i+1;
    z_snow(t-23:t) = z_snow_manual(i);
end

%{
Soil moisture measurements are not available at Abisko, so a proxy is
created using ground water level. However, these measurements are only
taken once in 2 two weeks. So, hourly time steps in between are linearly
interpolated, as snow cover was persistent (and soil frozen?) so that large
changes are unlikely. For each manual measurement maximum and minimum
values are given as well, calculated from automatic measurements (multiple
times per day, but not available?). The extrema are used to provide a scale
for percentages. All-time minima and maxima are used to calculate a range
representing 0% and 100% soil moisture.
Two locations of ground water level measurements are next to the
researchstation at Abisko, and the average is used here.
%}
Max_1 = -2.2;   Min_1 = -6.55;
Max_2 = 0.19;   Min_2 = -5.08;
GWL = nan(46*24,1);
GWL(12) = ((-6.22-Min_1)/(Max_1-Min_1) + (-4.16-Min_2)/(Max_2-Min_2))/2;  % 1 March 12:00
GWL(348) = ((-6.5-Min_1)/(Max_1-Min_1) + (-4.37-Min_2)/(Max_2-Min_2))/2; % 15 March 12:00
GWL(756) = ((-6.5-Min_1)/(Max_1-Min_1) + (-4.75-Min_2)/(Max_2-Min_2))/2;   % 1 April 12:00
GWL(1092) = ((-6.45-Min_1)/(Max_1-Min_1) + (-4-Min_2)/(Max_2-Min_2))/2;% 15 April 12:00
for t=13:347
    GWL(t) = GWL(t-1) + (GWL(348)-GWL(12))/(348-12-1);
end
for t=349:755
    GWL(t) = GWL(t-1) + (GWL(756)-GWL(348))/(756-348-1);
end
for t=757:1091
    GWL(t) = GWL(t-1) + (GWL(1092)-GWL(756))/(1092-756-1);
end
GWL = GWL(9*24+15:(31+4)*24+8);


%-------------------------------------------------------------------------%
%-------------------------  meteorological data  -------------------------%
%{
AWS data: 10-min averages from 1 March 0:10 to 30 April 2011 0:00
Cut out 14:10 10 March to 8:00 5 April 2011.
2 different air temperatures (on at RH-sensor). Usually 0.5C warmer at
RH-sensor, but also up to 2C.
%}
bla = importdata('Abisko_AWS_Logger1.dat');
T_air_RH = bla.data(1381:5088,18)+273.15;;
RH = bla.data(1381:5088,17);
Wind = bla.data(1381:5088,5);
clear bla
bla = importdata('Abisko_AWS_Logger2.dat');
T_air = bla.data(1381:5088,5)+273.15;
T_soil_5 = bla.data(1381:5088,6)+273.15;
T_soil_20 = bla.data(1381:5088,7)+273.15;
T_soil_50 = bla.data(1381:5088,8)+273.15;
T_soil_100 = bla.data(1381:5088,9)+273.15;
Precip_AWS = bla.data(1381:5088,13);
clear bla
T_soil = nan(length(T_soil_5),4);
T_soil(:,1) = T_soil_5;
T_soil(:,2) = T_soil_20;
T_soil(:,3) = T_soil_50;
T_soil(:,4) = T_soil_100;
clear T_soil_5 T_soil_20 T_soil_50 T_soil_100

% calculate hourly averages
i=0;
for t=6:6:length(RH)
    i=i+1;
    T_air_RH_1h(i) = mean(T_air_RH(t-5:t));
    RH_1h(i) = mean(RH(t-5:t));
    Wind_1h(i) = mean(Wind(t-5:t));
    T_air_1h(i) = mean(T_air(t-5:t));
    for j=1:4
        T_soil_1h(i,j) = mean(T_soil(t-5:t,j));
    end
    Precip_AWS_1h(i) = sum(Precip_AWS(t-5:t));
end

bla = xlsread('SnowDepth.xlsx');
Precip_manual = bla(10:36,1);
clear bla
for p=1:length(Precip_manual)
    if isnan(Precip_manual(p))
        Precip_manual(p) = 0;
    end
end

%{
from Annika Kristoffersson:
"You will also find precipitation and snow depth data from our manual
readings. The precipitation from the AWS is not reliable, specially during winter."
Manual measurements daily not hourly, so we use precip from the AWS for
timing throughout a day but use manual measurements to determine the
magnitude. For 10 March and 5 April coverage is not full 24h, however,
precipitation on those days is 0. For precipitation events not covered by
AWS, precip is spread throughout the day.
%}
Precip_adjust = zeros(size(Precip_AWS_1h));
i=1;    % 10 March set to 0
for t=10+24:24:length(Precip_adjust)-8
    i=i+1;
    Precip_sum = sum(Precip_AWS_1h(t-23:t));
    if Precip_sum > 0
        for h=1:24
            Precip_adjust(t-(h-1)) = Precip_AWS_1h(t-(h-1))/Precip_sum * Precip_manual(i);
        end
    elseif Precip_sum == 0 && Precip_manual(i) > 0
        Precip_adjust(t-23:t) = Precip_manual(i)/24;
    end
end

save('ForcingData_ToyModel_Abisko.mat','SW_in_open_1h','f_SW_diff_1h',...
    'LW_in_open_1h','LW_in_bc_C_1h','SW_in_bc_C_avg_1h','T_veg_TC_C_avg_1h',...
    'T_surf_C','T_soil_1h','GWL','z_snow','Wind_1h','T_air_RH_1h','RH_1h',...
    'Precip_adjust','time_1h')
toc