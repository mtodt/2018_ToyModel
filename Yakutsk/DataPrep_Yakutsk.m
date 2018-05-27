tic
clear
close all
%-------------------------------------------------------------------------%
%-----------------------------  import data  -----------------------------%
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%------------------------------  open site  ------------------------------%
%{
Precipitation from Carbon Dioxide Information Analysis Center (CDIAC).
"Jakutsk" station at 62°01'N and 129°43'E, elevation 98m.
Format:
- station index
- year
- month
- day
- temperature group quality flag
- temperature minimum
- quality flag
- temperature mean
- quality flag
- temperature maximum
- quality flag
- daily total precipitation (liquid equivalent) [mm] "...snowfall being
  converted to a liquid total by melting the snow in the gauge."
- additional flag: 0 - "Measured precipitation amount equal to or higher
                        than 0.1mm."
                   1 - "Precipitation measured for a few days."
                   2 - "Precipitation measurements were made, but no
                        precipitation occurred (R = 0)."
                   3 - "Only precipitation traces were observed (< 0.1mm) R=0"
                   9 - "The value is rejected or no observations were made."
- quality flag: 0 - "The value is reliable."
                9 - "The value is rejected or no observations were made."
%}
bla = load('Precip.txt');
    % daily values -> select 1 January 1998 - 31 December 2000
Precip = bla(40179:41274,12);
Precip_flags = bla(40179:41274,13:14);
clear bla

% quality control
for l=1:length(Precip)
    if Precip_flags(l,1) == 9
        Precip(l) = nan;
    elseif Precip_flags(l,2) == 9
        Precip(l) = nan;
    end
end

% distribute daily values equally over 24 hours
Precip_1h = nan(length(Precip)*24,1);
p=0;
for l=1:length(Precip)
    p=p+24;
    Precip_1h(p-23:p) = Precip(l)/24;
end

% split into years as no continuous simulation possible
Precip_98_1h = Precip_1h(1:365*24);


%-------------------------------------------------------------------------%
%-----------------------------  forest data  -----------------------------%

%{
Import hourly data from Larch-1998.csv, however, the file unfortunately
contains row (and columns) of nans only, which forces us to use this more
tedious approach of reading it in.
%}
fileID = fopen('Larch-1998.csv');
formatSpec = '%f%f%f%s%f%s%f%s%f%s%f%s%f%f%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%s%s%f%s%f%s%f%f%f%f%f%f%f%f%f%f%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%f%s%s%f%s%f%f%f%f%f%f%f%f%f%s%s%[^\n\r]';
bla = textscan(fileID,formatSpec,'Delimiter',',','EmptyValue',NaN,'HeaderLines',1,'ReturnOnError', false);
fclose(fileID);
clear fileID formatSpec

% time
YEAR = bla{1,1};
Date = bla{1,2};
UTC = bla{1,3};
time_98_1h = nan(size(Date));
time_98_1h(1) = datenum(1998,1,1,1,0,0);
for t=2:length(time_98_1h)
    time_98_1h(t) = time_98_1h(t-1) + 1/24;
end
clear YEAR Date UTC

% wind speed [m s^{-1}]
Wind_98_1h = str2double(bla{1,4});
    % Measurements at 32m end in late October, however, at 27m do not.
Wind_98_1h = str2double(bla{1,6});

% relative humidity [%]
RH_98_1h = str2double(bla{1,18});
    % Measurements end on 30 June, however, do not in individual file.
RH_98_1h = csvread('Lv1_32m\Rh\LRh98hr-1.csv',2,5,[2,5,8761,5]);
RH_98_1h = NaN999(RH_98_1h);

% air temperature [°C]
T_air_98_1h = str2double(bla{1,26}); T_air_98_1h = T_air_98_1h + 273.15;

% surface temperature [°C]
T_surf_98_1h = str2double(bla{1,36}); T_surf_98_1h = T_surf_98_1h + 273.15;

% soil temperatures [°C] (at 0cm, 10cm, 20cm, 40cm, 60cm, 80cm, and 120cm ?)
T_soil_0_98_1h = str2double(bla{1,38});
T_soil_10_98_1h = str2double(bla{1,40});
T_soil_20_98_1h = str2double(bla{1,42});
T_soil_40_98_1h = str2double(bla{1,44});
T_soil_60_98_1h = str2double(bla{1,46});
T_soil_80_98_1h = str2double(bla{1,48});
T_soil_120_98_1h = str2double(bla{1,50});
T_soil_98_1h = horzcat(T_soil_0_98_1h,T_soil_10_98_1h,T_soil_20_98_1h,...
    T_soil_40_98_1h,T_soil_60_98_1h,T_soil_80_98_1h,T_soil_120_98_1h);
T_soil_98_1h = T_soil_98_1h + 273.15;
clear T_soil_0_98_1h T_soil_10_98_1h T_soil_20_98_1h T_soil_40_98_1h
clear T_soil_60_98_1h T_soil_80_98_1h T_soil_120_98_1h

% soil moisture [%] (at 10cm, 20cm, 40cm, 60cm, and 80cm)
SM_01m = bla{1,54};
SM_02m = bla{1,56};
SM_04m = bla{1,58};
SM_06m = bla{1,60};
SM_08m = bla{1,62};
SM_98_1h = horzcat(SM_01m,SM_02m,SM_04m,SM_06m,SM_08m);
clear SM_01m SM_02m SM_04m SM_06m SM_08m

% above-canopy radiation [W m^{-2}]
SW_in_ac_98_1h = str2double(bla{1,76});
SW_out_ac_98_1h = str2double(bla{1,80});
LW_in_ac_98_1h = str2double(bla{1,64});
LW_out_ac_98_1h = str2double(bla{1,68});
Rnet_ac_98_1h = str2double(bla{1,72});

% sub-canopy radiation [W m^{-2}]
SW_in_bc_98_1h = str2double(bla{1,78});     % ends in mid-November
SW_out_bc_98_1h = str2double(bla{1,82});
LW_in_bc_98_1h = str2double(bla{1,66});     % empty!!!
LW_out_bc_98_1h = str2double(bla{1,70});    % empty!!!
Rnet_bc_98_1h = str2double(bla{1,74});      % ends in early Augsut
clear bla


%-------------------------------------------------------------------------%
%----------------------------  interpolation  ----------------------------%
%-------------------------------------------------------------------------%
%{
As was done by Essery et al. (2016) for data at Sodankylä, gaps of 4 hours
or shorter are filled by linear interpolation. At least for shortwave
radiation an approach based on diurnal cycles from previous days could be
applied.
%}
% wind speed
Wind_98_1h_filled = LinearInterp_4h(Wind_98_1h);

% relative humidity
RH_98_1h_filled = LinearInterp_4h(RH_98_1h);

% air temperature
T_air_98_1h_filled = LinearInterp_4h(T_air_98_1h);

% atmospheric LWR
LW_in_ac_98_1h_filled = LinearInterp_4h(LW_in_ac_98_1h);

% SWR
SW_in_ac_98_1h_filled = LinearInterp_4h(SW_in_ac_98_1h);
SW_in_bc_98_1h_filled = LinearInterp_4h(SW_in_bc_98_1h);
SW_out_bc_98_1h_filled = LinearInterp_4h(SW_out_bc_98_1h);

% sub-canopy net radiation
Rnet_bc_98_1h_filled = LinearInterp_4h(Rnet_bc_98_1h);

% surface temperature
T_surf_98_1h_filled = LinearInterp_4h(T_surf_98_1h);

% soil temperature
T_soil_98_1h_filled = nan(size(T_soil_98_1h));
for i=1:length(T_soil_98_1h(1,:))
    T_soil_98_1h_filled(:,i) = LinearInterp_4h(T_soil_98_1h(:,i));
end

% soil moisture
SM_98_1h_filled = nan(size(SM_98_1h));
for i=1:length(SM_98_1h(1,:))
    SM_98_1h_filled(:,i) = LinearInterp_4h(SM_98_1h(:,i));
end

%{
There is still a gap of about 24 hours on 14/15 March 1998, which we
exclude from simulations. Simulation stops before gap and continues
afterwards with variables required from previous timestep taken from before
the gap. Since evaluation will start 1:00 there is enough time for spin-up.
%}
gap_98_1h = 1744:1764;


%-------------------------------------------------------------------------%
%--------------------------  prepare soil data  --------------------------%
%-------------------------------------------------------------------------%

%----------------  use top soil temperature individually  ----------------%
T_topsoil_98_1h_filled = T_soil_98_1h_filled(:,1);

%------------------  average soil variables vertically  ------------------%
% soil moisture
%{
One level includes nans for some time in 1998, so we have to calculate
non-nan depth for every time step.
%}
levels = [0.15 0.15 0.2 0.2 0.2];
temp1 = nan(size(SM_98_1h_filled));
temp2 = nan(size(SM_98_1h_filled));
SM_98_1h_filled_avg = nan(size(SM_98_1h_filled(:,1)));
for t=1:length(SM_98_1h_filled_avg)
    for i=1:length(levels)
        temp1(t,i) = SM_98_1h_filled(t,i)*levels(i);
        if ~isnan(temp1(t,i))
            temp2(t,i) = levels(i);
        end
    end
    SM_98_1h_filled_avg(t) = nansum(temp1(t,:))/nansum(temp2(t,:));
end
clear temp1 temp2

% soil temperature
%{
One level starts later than others in 1998, so we have to calculate
non-nan depth for every time step.
%}
boundaries = [0 0.1 0.2 0.4 0.6 0.8 1.2];
levels = boundaries(2:end)-boundaries(1:end-1);
temp = (T_soil_98_1h_filled(:,2:end)+T_soil_98_1h_filled(:,1:end-1))/2;
temp1 = nan(size(T_soil_98_1h_filled));
temp2 = nan(size(T_soil_98_1h_filled));
T_soil_98_1h_filled_avg = nan(size(T_soil_98_1h_filled(:,1)));
for t=1:length(T_soil_98_1h_filled_avg)
    for i=1:length(levels)
        temp1(t,i) = temp(t,i)*levels(i);
        if ~isnan(temp1(t,i))
            temp2(t,i) = levels(i);
        end
    end
    T_soil_98_1h_filled_avg(t) = nansum(temp1(t,:))/nansum(temp2(t,:));
end
clear temp1 temp2

% calculate fraction of frozen soil
T_proxy = T_soil_98_1h_filled;
for t=1:length(T_soil_98_1h_filled_avg)
    for i=2:length(boundaries)-1
        if isnan(T_soil_98_1h_filled(t,i))
            gradient = (T_soil_98_1h_filled(t,i+1)-T_soil_98_1h_filled(t,i-1))/(levels(i-1)+levels(i));
            T_proxy(t,i) = T_soil_98_1h_filled(t,i-1) + gradient*levels(i-1);
        end
    end
end
temp = nan(length(T_soil_98_1h_filled_avg),length(levels));
T_freez = 273.15;
frac_soil_frozen_98_1h = nan(size(T_soil_98_1h_filled_avg));
for t=1:length(T_soil_98_1h_filled_avg)
    for i=1:length(levels)
        if T_proxy(t,i) <= T_freez && T_proxy(t,i+1) <= T_freez
            temp(t,i) = levels(i);
        elseif T_proxy(t,i) > T_freez && T_proxy(t,i+1) > T_freez
            temp(t,i) = 0;
        else
            mintemp = min(T_proxy(t,i),T_proxy(t,i+1));
            maxtemp = max(T_proxy(t,i),T_proxy(t,i+1));
            temp(t,i) = (T_freez-mintemp)/(maxtemp-mintemp) * levels(i);
        end
    end
    frac_soil_frozen_98_1h(t) = sum(temp(t,:))/sum(levels);
end
clear T_proxy temp mintemp maxtemp


%-------------------------------------------------------------------------%
%------------------------------  save data  ------------------------------%
%-------------------------------------------------------------------------%
save('ForcingData_ToyModel_Yakutsk.mat','Precip_98_1h','RH_98_1h_filled',...
    'Wind_98_1h_filled','T_air_98_1h_filled','LW_in_ac_98_1h_filled',...
    'SW_in_ac_98_1h_filled','SW_in_bc_98_1h_filled','SW_out_bc_98_1h_filled',...
    'Rnet_bc_98_1h_filled','T_surf_98_1h_filled','T_topsoil_98_1h_filled',...
    'frac_soil_frozen_98_1h','SM_98_1h_filled_avg','time_98_1h','gap_98_1h')
toc