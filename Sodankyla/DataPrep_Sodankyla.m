tic
clear
close all
%-------------------------------------------------------------------------%
%-----------------------------  forest data  -----------------------------%
sod_data = xlsread('AllSodankylaData 270616.xls');

% open = top of met tower (but not 4Xm one)
SW_in_open = sod_data(280:end,142);
LW_in_open = sod_data(280:end,143);
%{
from Nick's script:
"Following Reid et al. 2014 (Appendix C) need to correct the BF3 data for
Sodankyla"
%}
sod_open_BF3tot = sod_data(280:end,140)/1.081;    % correcting total irradiance
sod_open_BF3diff = sod_data(280:end,141)/1.107;   % correcting diffuse irradiance
f_SW_diff = ones(size(sod_open_BF3diff));
for l=1:length(f_SW_diff)
    if sod_open_BF3tot(l) > 0
        f_SW_diff(l) = min(1,sod_open_BF3diff(l)/sod_open_BF3tot(l));
    end
end

% continuous site
SW_in_bc_C = sod_data(280:end,144:153);
LW_in_bc_C = sod_data(280:end,154:157);
T_air_sub_C = sod_data(280:end,73)+273.15;
T_veg_IR_C = sod_data(280:end,8)+273.15;
T_veg_TC_C = sod_data(280:end,9:72)+273.15;
%{
Based on descriptions in Nick's script, some instruments are dismissed for
several reasons:
- unreliable, and
- imbalance N-S (or E-W, if available).
%}
T_veg_TC_C_QC = horzcat(T_veg_TC_C(:,1:28),T_veg_TC_C(:,30:35),...
    T_veg_TC_C(:,38:39),T_veg_TC_C(:,43:44),T_veg_TC_C(:,47:60));

% roving sites
SW_in_bc_R4 = sod_data(280:end,158:167);
LW_in_bc_R4 = sod_data(280:end,168:171);
T_air_sub_R4 = sod_data(280:end,139)+273.15;
T_veg_IR_R4 = sod_data(280:end,74)+273.15;
T_veg_TC_R4 = sod_data(280:end,75:138)+273.15;
%{
Based on descriptions in Nick's script, some instruments are dismissed for
several reasons:
- unreliable,
- imbalance N-S (or E-W, if available), and
- birches for site R4 as we assume a needleleaf evergreen PFT since
  fractions of different vegetation types are not available.
%}
T_veg_TC_R4_QC = horzcat(T_veg_TC_R4(:,1:4),T_veg_TC_R4(:,11:14),T_veg_TC_R4(:,35:46));

%{
Site averaging:
- 4 LW-radiometers per site - but maybe not used because of varying VAI
- 10 SW-radiometers per site - averaged since not required for model
- 64 thermocouplets per site, although varying number active
%}
SW_in_bc_C_avg = nan(size(SW_in_bc_C(:,1)));
SW_in_bc_R4_avg = nan(size(SW_in_bc_R4(:,1)));
LW_in_bc_C_avg = nan(size(LW_in_bc_C(:,1)));
LW_in_bc_R4_avg = nan(size(LW_in_bc_R4(:,1)));
T_veg_TC_C_avg = nan(size(T_veg_TC_C_QC(:,1)));
T_veg_TC_R4_avg = nan(size(T_veg_TC_R4(:,1)));
for l=1:length(LW_in_open)
    SW_in_bc_C_avg(l) = mean(SW_in_bc_C(l,:));
    SW_in_bc_R4_avg(l) = mean(SW_in_bc_R4(l,:));
    
    LW_in_bc_C_avg(l) = mean(LW_in_bc_C(l,:));
    LW_in_bc_R4_avg(l) = mean(LW_in_bc_R4(l,:));
    
    T_veg_TC_C_avg(l) = nanmean(T_veg_TC_C_QC(l,:));
    T_veg_TC_R4_avg(l) = nanmean(T_veg_TC_R4_QC(l,:));
end

% hourly averages
%{
Translating measurements from temporal resolution of 5 minutes into hourly
values. Coverages starts at 13:00 on 9 March and ends at 8:00 on 25 April,
with hourly values set to the end of the averaging period, i.e. 15:00 for
the first time step, which starts at 13:05. Evaluation period should then
be set to 1:00 on 10 March until 0:00 on 25 April.
%}
i=0;
for t=13:12:length(LW_in_open)
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
    T_veg_IR_C_1h(i) = mean(T_veg_IR_C(t-11:t));
    T_veg_TC_C_avg_1h(i) = mean(T_veg_TC_C_avg(t-11:t));
    SW_in_bc_R4_avg_1h(i) = mean(SW_in_bc_R4_avg(t-11:t));
    LW_in_bc_R4_avg_1h(i) = mean(LW_in_bc_R4_avg(t-11:t));
    for n=1:4
        LW_in_bc_R4_1h(i,n) = mean(LW_in_bc_R4(t-11:t,n));
    end
    T_air_sub_R4_1h(i) = mean(T_air_sub_R4(t-11:t));
    T_veg_IR_R4_1h(i) = mean(T_veg_IR_R4(t-11:t));
    T_veg_TC_R4_avg_1h(i) = mean(T_veg_TC_R4_avg(t-11:t));
    time_1h(i) = datenum(sod_data(279+t,1),sod_data(279+t,2),sod_data(279+t,3),...
        sod_data(279+t,5),sod_data(279+t,6),sod_data(279+t,7));
end

%{
Open area radiation measurements feature a gap of about 90 min. Since this
occurs in early afternoon a linear interpolation of hourly values is
applied. These timesteps should then be excluded from evaluation.
%}
for t=457:459
    LW_in_open_1h(t) = LW_in_open_1h(t-1) + (LW_in_open_1h(460)-LW_in_open_1h(456))/4;
    SW_in_open_1h(t) = SW_in_open_1h(t-1) + (SW_in_open_1h(460)-SW_in_open_1h(456))/4;
end


%-------------------------------------------------------------------------%
%-----------------------------  ground data  -----------------------------%

% snow surface temperature from sub-canopy air temperature capped at 0C
T_surf_C = nan(size(T_air_sub_C_1h)); T_surf_R4 = nan(size(T_air_sub_R4_1h));
for t=1:length(time_1h)
    if ~isnan(T_air_sub_C_1h(t))
        T_surf_C(t) = min(T_air_sub_C_1h(t),273.15);
    end
    if ~isnan(T_air_sub_R4_1h(t))
        T_surf_R4(t) = min(T_air_sub_R4_1h(t),273.15);
    end
end

%{
Soil moisture from IOA site (as fraction, not vol % as stated in
description) as well as soil temperature.
Three locations but 2 & 3 only 5cm and 10cm, location 1 down to 80cm.
Already hourly averages. In accordance with sub-canopy measurements,
observations are cut out for period 14:00 on 9 March until 8:00 on 25 April.
This allows for a spin-up period of 10 hours, which should be sufficient
for the Toy Model.
Strangely, soil moisture and temperature increase rapidly around 12-13
April, although supposedly still snow cover. Then again, not necessarily
same location.
%}
bla = importdata('IOA_in_situ_2011.csv');
%bla = xlsread('IOA_in_situ_2011.csv','AK:AO');
SM = bla.data(4575:5697,36:40);
T_soil = bla.data(4575:5697,45:49)+273.15;
clear bla

% average soil moisture vertically
SM_avg = nan(size(squeeze(SM(:,1))));
for l=1:length(SM_avg)
    SM_avg(l) = (SM(l,1)*7.5 + SM(l,2)*7.5 + SM(l,3)*15 + SM(l,4)*30 + SM(l,5)*40)/100;
end

% calculate fraction of frozen soil
boundaries = [0.05 0.1 0.2 0.4 0.8];
levels = boundaries(2:end)-boundaries(1:end-1);
temp = nan(length(T_soil(:,1)),length(levels));
T_freez = 273.15;
frac_soil_frozen = nan(size(T_soil(:,1)));
for t=1:length(frac_soil_frozen)
    for i=1:length(levels)
        if T_soil(t,i) <= T_freez && T_soil(t,i+1) <= T_freez
            temp(t,i) = levels(i);
        elseif T_soil(t,i) > T_freez && T_soil(t,i+1) > T_freez
            temp(t,i) = 0;
        else
            mintemp = min(T_soil(t,i),T_soil(t,i+1));
            maxtemp = max(T_soil(t,i),T_soil(t,i+1));
            temp(t,i) = (T_freez-mintemp)/(maxtemp-mintemp) * levels(i);
        end
    end
    frac_soil_frozen(t) = sum(temp(t,:))/sum(levels);
end
clear temp mintemp maxtemp


%-------------------------------------------------------------------------%
%-------------------------  meteorological data  -------------------------%
%{
Import air temperature, relative humidity and wind speed from
meteorological tower at 47m height, already cut out for desired period.
However, hourly data start with 0:00 and end with 23:00, although usually
aggregated from previous hour. (Shift by 1 hour?)
In accordance with sub-canopy measurements, observations are cut out for
period 14:00 on 9 March until 8:00 on 25 April. This allows for a spin-up
period of 10 hours, which should be enough for the Toy Model.
%}
bla = importdata('CoSDAS_winter_11_12\met_mast_data.txt');
Wind_mast = bla.data(3855:4977,7);              % [m s^{-1}]
T_air_mast = bla.data(3855:4977,11)+273.15;     % [K]
RH_mast =  bla.data(3855:4977,16);              % [%]
clear bla

% some gaps, but only individual hours - linear interpolation
for l=2:length(Wind_mast)-1
    if isnan(Wind_mast(l))
        Wind_mast(l) = Wind_mast(l-1) + (Wind_mast(l+1)-Wind_mast(l-1))/2;
    end
    if isnan(T_air_mast(l))
        T_air_mast(l) = T_air_mast(l-1) + (T_air_mast(l+1)-T_air_mast(l-1))/2;
    end
    if isnan(RH_mast(l))
        RH_mast(l) = RH_mast(l-1) + (RH_mast(l+1)-RH_mast(l-1))/2;
    end
end

bla = importdata('CoSDAS_winter_11_12\met_data.txt');
z_snow_open =  bla.data(3855:4977,12); % [cm]
    z_snow_open =  z_snow_open/100;
Prec_met_weigh =  bla.data(3855:4977,13); % [mm h^{-1}]
    for l=1:length(Prec_met_weigh)
        if Prec_met_weigh(l) < 0
            Prec_met_weigh(l) = 0;
        elseif isnan(Prec_met_weigh(l))
            Prec_met_weigh(l) = 0;
        end
    end
Prec_met_present =  bla.data(3855:4977,14); % [mm h^{-1}]
    for l=1:length(Prec_met_present)
        if Prec_met_present(l) < 0
            Prec_met_present(l) = 0;
        elseif isnan(Prec_met_present(l))
            Prec_met_present(l) = 0;
        end
    end
clear bla

% corr between dZsnow/dt and precip: 0.22 for _present and 0.31 for _weigh
save('ForcingData_ToyModel_Sodankyla.mat','SW_in_open_1h','f_SW_diff_1h',...
    'LW_in_open_1h','LW_in_bc_C_1h','LW_in_bc_R4_1h','SW_in_bc_C_avg_1h',...
    'SW_in_bc_R4_avg_1h','T_surf_C','T_surf_R4','T_veg_IR_C_1h','T_veg_TC_C_avg_1h',...
    'T_veg_IR_R4_1h','T_veg_TC_R4_avg_1h','frac_soil_frozen','SM_avg','z_snow_open',...
    'Wind_mast','T_air_mast','RH_mast','Prec_met_weigh','Prec_met_present','time_1h')
toc