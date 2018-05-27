tic
clear
close all
%-------------------------------------------------------------------------%
%-----------------------------  import data  -----------------------------%
%-------------------------------------------------------------------------%

%------------------------  forest radiation data  ------------------------%
filename = 'C:\Users\w15044103\Documents\Modelling\Toy Model\Experiments\Cherskiy\netR.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%q%q%f%[^\n\r]';
fileID = fopen(filename,'r');
bla = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,...
    NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
DoY = bla{:,1};
year = bla{:,2};
hour = bla{:,3};
minute = bla{:,4};
RECORD = bla{:,5};
SW_in = bla{:,6};
SW_out = bla{:,7};
for l=1:length(DoY)
    if SW_in(l) < 0
        SW_in(l) = 0;
    end
    if SW_out(l) < 0
        SW_out(l) = 0;
    end
end
LW_in = bla{:,8};
LW_out = bla{:,9};
site = bla{:,10};
loc = bla{:,11};
sensorZ = bla{:,12};
clear filename delimiter startRow formatSpec fileID bla ans

i_ac = 0;
i_bc = 0;
for l=1:length(DoY)
    if strcmp(site(l),'hd') == 1    % only use "high density" stand
        if strcmp(loc(l),'canopy') == 1
            i_ac = i_ac + 1;
            [day,month] = DoY_to_DayNMonth(year(l),DoY(l));
            time_rad_ac(i_ac) = datenum(year(l),month,day,floor(hour(l)),minute(l),0);
            DoY_ac(i_ac) = DoY(l);
            JulDay_ac(i_ac) = DoY(l) + hour(l)/24;
            year_ac(i_ac) = year(l);
            SW_in_ac(i_ac) = SW_in(l);
            SW_out_ac(i_ac) = SW_out(l);
            LW_in_ac(i_ac) = LW_in(l);
            LW_out_ac(i_ac) = LW_out(l);
        elseif strcmp(loc(l),'understory') == 1
            i_bc = i_bc + 1;
            [day,month] = DoY_to_DayNMonth(year(l),DoY(l));
            time_rad_bc(i_bc) = datenum(year(l),month,day,floor(hour(l)),minute(l),0);
            DoY_bc(i_bc) = DoY(l);
            JulDay_bc(i_bc) = DoY(l) + hour(l)/24;
            year_bc(i_bc) = year(l);
            SW_in_bc(i_bc) = SW_in(l);
            SW_out_bc(i_bc) = SW_out(l);
            LW_in_bc(i_bc) = LW_in(l);
            LW_out_bc(i_bc) = LW_out(l);
        end
    end
end

% four more timesteps for sub-canopy measurements -> align
i=0;
for l=1:length(time_rad_ac)
    bla = find(time_rad_bc==time_rad_ac(l),1);
    if ~isnan(bla)
        i=i+1;
        idx(i,1) = l;       % sanity check
        idx(i,2) = bla;     % sanity check
        time_rad_clean(i) = time_rad_ac(l);
        DoY_rad_clean(i) = DoY_ac(l);
        JulDay_rad_clean(i) = JulDay_ac(l);
        year_rad_clean(i) = year_ac(l);
        SW_in_ac_clean(i) = SW_in_ac(l);
        SW_out_ac_clean(i) = SW_out_ac(l);
        LW_in_ac_clean(i) = LW_in_ac(l);
        LW_out_ac_clean(i) = LW_out_ac(l);
        SW_in_bc_clean(i) = SW_in_bc(bla);
        SW_out_bc_clean(i) = SW_out_bc(bla);
        LW_in_bc_clean(i) = LW_in_bc(bla);
        LW_out_bc_clean(i) = LW_out_bc(bla);
    end
end

% get value for soil albedo
alb_surf = SW_out_bc_clean./SW_in_bc_clean;
for l=1:length(alb_surf)
    if alb_surf(l) > 1
        alb_surf(l) = 1;
    end
end
%{
    % concatenate midday values from different periods
alb_surf_midday = horzcat(alb_surf(34:48:34+46*48),...
    alb_surf(2304:48:2304+63*48),alb_surf(5405:48:5405+136*48));
    % create histogram
spectrum_alb = 0:0.025:1;
spectrum_alb_xaxis = 0.0125:0.025:0.9875;
hist_alb = histogram(alb_surf_midday,spectrum_alb,'Normalization','probability');
hist_ground = hist_alb.Values;
fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
plot(spectrum_alb_xaxis,hist_ground,'k','LineWidth',2)
xlim([0 1])
% ylim([0 0.08])
xlabel('Albedo','FontSize',17,'FontWeight','bold')
ylabel('Probablity','FontSize',17,'FontWeight','bold')
title('Mid-day Ground Albedo at Cherskiy (limited to 1)','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0:0.1:1)
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 7];
fig.PaperSize = [8 7];
print(fig,'-dpdf','-r600','AlbedoGroundPDF_Cherskiy.pdf')
%}

%-------------------------  meteorological data  -------------------------%
filename = 'C:\Users\w15044103\Documents\Modelling\Toy Model\Experiments\Cherskiy\met.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
    'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
    
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,11]
% Converts strings in the input cell array to numbers. Replaced non-numeric
% strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
% Create a regular expression to detect and remove non-numeric prefixes and
% suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
% Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
% Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

rawNumericColumns = raw(:, [1,2,3,4,5,6,11]);
rawCellColumns = raw(:, [7,8,9,10]);
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
DoY = cell2mat(rawNumericColumns(:, 1));
year = cell2mat(rawNumericColumns(:, 2));
hour = cell2mat(rawNumericColumns(:, 3));
RH = cell2mat(rawNumericColumns(:, 4));
    RH = RH*100;    % convert to %
sensor_height = cell2mat(rawNumericColumns(:, 5));
T_air = cell2mat(rawNumericColumns(:, 7));
    T_air = T_air + 273.15; % convert to K
clear filename delimiter startRow formatSpec fileID dataArray ans raw 
clear col numericData rawData row regexstr result numbers 
clear invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R
minute_met = (hour - floor(hour))*60;
hour_met = hour;
DoY_met = DoY;
JulDay_met = DoY + hour/24;
year_met = year;
clear DoY hour year

% gap of one time step - linear interpolation
minute_met(15500) = 0;
hour_met(15500) = 16;
DoY_met(15500) = 138;
year_met(15500) = 2017;
T_air(15500) = sum(T_air(15499)+T_air(15501))/2;
RH(15500) = sum(RH(15499)+RH(15501))/2;

time_met = nan(size(DoY_met));
for l=1:length(DoY_met)
    [day,month] = DoY_to_DayNMonth(year_met(l),DoY_met(l));
    time_met(l) = datenum(year_met(l),month,day,floor(hour_met(l)),minute_met(l),0); 
end

%--------------------------  data from airport  --------------------------%
%{
Wind speed and precipitation aren't available at the forest site, so we use
data from the met station at the airport of Cherskiy. However, temporal
resolution is only 3 hours.
%}
filename = 'C:\Users\w15044103\Documents\Modelling\Toy Model\Experiments\Cherskiy\airport_met.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%f%f%f%f%f%q%[^\n\r]';
fileID = fopen(filename,'r');
bla = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,...
    NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
DoY = bla{:,1};
year = bla{:,2};
hour = bla{:,3};
Precip = bla{:,4};
Wind = bla{:,5};
clear filename delimiter startRow formatSpec fileID bla
minute_airport = (hour - floor(hour))*60;
hour_airport = hour;
DoY_airport = DoY;
JulDay_airport = DoY + hour/24;
year_airport = year;
clear DoY hour year
time_airport = nan(size(DoY_airport));
for l=1:length(DoY_airport)
    [day,month] = DoY_to_DayNMonth(year_airport(l),DoY_airport(l));
    time_airport(l) = datenum(year_airport(l),month,day,floor(hour_airport(l)),minute_airport(l),0); 
end

% time series for airport variables go backward in time -> flip upside down
time_airport = flipud(time_airport);
Precip = flipud(Precip);
Wind = flipud(Wind);

%----------------------------  soil moisture  ----------------------------%
filename = 'C:\Users\w15044103\Documents\Modelling\Toy Model\Experiments\Cherskiy\vwc_profile.csv';
delimiter = ',';
startRow = 2;
formatSpec = '%s%{MM/dd/yyyy}D%f%f%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
bla = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,...
    startRow-1, 'ReturnOnError', false);
fclose(fileID);
DoY = bla{:,3};
year = bla{:,4};
samplenumber = bla{:,5};
depthmoss = bla{:,6};
mossvwcmineral = bla{:,7};
organicvwcmineral = bla{:,8};
clear filename delimiter startRow formatSpec fileID bla
DoY_soil = DoY;
year_soil = year;
clear DoY year
time_soil = nan(size(DoY_soil));
for l=1:length(DoY_soil)
    [day,month] = DoY_to_DayNMonth(year_soil(l),DoY_soil(l));
    time_soil(l) = datenum(year_soil(l),month,day,12,0,0); 
end


%-------------------------------------------------------------------------%
%---------------------------  synchronization  ---------------------------%
%-------------------------------------------------------------------------%
%{
Radiation measurements consist of three periods, 29 June - 15 August 2016,
30 August - 3 November 2016 and 29 March - 13 August 2017. Since we're only
interested in snowmelt we choose the last part only (and cut it later to a 
shorter period).
%}
time_rad_2017 = time_rad_clean(5360:end);
JulDay_rad_2017 = JulDay_rad_clean(5360:end);
DoY_rad_2017 = DoY_rad_clean(5360:end);
year_rad_2017 = year_rad_clean(5360:end);
SW_in_ac_2017 = SW_in_ac_clean(5360:end);
SW_out_ac_2017 = SW_out_ac_clean(5360:end);
SW_in_bc_2017 = SW_in_bc_clean(5360:end);
SW_out_bc_2017 = SW_out_bc_clean(5360:end);
alb_surf_2017 = alb_surf(5360:end);
LW_in_ac_2017 = LW_in_ac_clean(5360:end);
LW_out_ac_2017 = LW_out_ac_clean(5360:end);
LW_in_bc_2017 = LW_in_bc_clean(5360:end);
LW_out_bc_2017 = LW_out_bc_clean(5360:end);

% select same period for air temperature and relative humidity
bla = find(time_met==time_rad_2017(1));
T_air_2017 = T_air(bla:bla-1+length(time_rad_2017));
RH_2017 = RH(bla:bla-1+length(time_rad_2017));

%----------------------  calculate hourly averages  ----------------------%
%{
Sub-canopy shortwave radiation measurements indicate rapid melt for 18 to
20 May 2017 with already low albedo before (between 0.3 and 0.6) but a
clear drop for those days. Snow cover might already have been patchy
before, however, albedo rarely stayed above 0.7 even in October & November
2016. This might mean that we have to use a lower "fresh snow albedo" for
the Toy Model instead of trying to adjust via snow cover fraction.
The end of simulation will be set to 21 May 2017 0:00 - the last day could
be excluded from evaluation - and start is determined by available data.
Radiation data start at 13:30 on 29 March 2017, so that 11 hours of spin-up
are available.
%}
i=0;
for t=2:2:2518
    i=i+1;
    time_1h(i) = time_rad_2017(t);
    JulDay_1h(i) = JulDay_rad_2017(t);
    DoY_1h(i) = DoY_rad_2017(t);
    SW_in_ac_1h(i) = mean(SW_in_ac_2017(t-1:t));
    SW_out_ac_1h(i) = mean(SW_out_ac_2017(t-1:t));
    SW_in_bc_1h(i) = mean(SW_in_bc_2017(t-1:t));
    SW_out_bc_1h(i) = mean(SW_out_bc_2017(t-1:t));
    alb_surf_1h(i) = mean(alb_surf_2017(t-1:t));
    LW_in_ac_1h(i) = mean(LW_in_ac_2017(t-1:t));
    LW_out_ac_1h(i) = mean(LW_out_ac_2017(t-1:t));
    LW_in_bc_1h(i) = mean(LW_in_bc_2017(t-1:t));
    LW_out_bc_1h(i) = mean(LW_out_bc_2017(t-1:t));
    T_air_1h(i) = mean(T_air_2017(t-1:t));
    RH_1h(i)= mean(RH_2017(t-1:t));
end

%{
Radiation and met forcing were measured every 30 mins, however, data from
the airport are only available every 3 hours. We assume values for
precipitation are sums and wind speeds are averaged over the previous 3
hours, and we set them constant for the two previous hours before each
value.
%}
Wind_1h = nan(size(time_1h));
Precip_1h = nan(size(time_1h));
for l=1:length(time_airport)
    bla = find(time_1h==time_airport(l));
    Wind_1h(bla) = Wind(l);        
    Precip_1h(bla) = Precip(l);
end
for t=4:3:length(Wind_1h)-1
    Wind_1h(t-2:t-1) = Wind_1h(t);
end
bla = find(time_airport==time_1h(end-1));
Wind_1h(end) = Wind(bla+1);
for t=4:3:length(Precip_1h)-1
    Precip_1h(t-2:t) = Precip_1h(t)/3;
end
Precip_1h(1) = Precip_1h(1)/3;
Precip_1h(end) = Precip(bla+1)/3;


%-------------------------------------------------------------------------%
%-----------------------  quality-controlling LWR  -----------------------%
%-------------------------------------------------------------------------%
%{
Check for constant values lasting several hours.
Only one case for above-canopy sensor, values ~295 W m^{-2} which fits air
temperature. Since we can't set met forcing to nan, as that would stop the
simulation, we already prescribe evaluation at this point.
%}
QualMeas = ones(size(time_1h));
QualMeas(1172:1179) = nan;


save('ForcingData_ToyModel_Cherskiy.mat','time_1h','JulDay_1h','DoY_1h',...
    'SW_in_ac_1h','SW_out_ac_1h','SW_in_bc_1h','SW_out_bc_1h','alb_surf_1h',...
    'LW_in_ac_1h','LW_out_ac_1h','LW_in_bc_1h','LW_out_bc_1h',...
    'T_air_1h','RH_1h','Wind_1h','Precip_1h','QualMeas')

toc