%{
Collect ingredients necessary for multi-linear regression.
- emissivity of the sky & solar angle -> description of meteorological
  conditions (insolation doesn't translate from simulation to simulation
  because of differences in latitude and season, and it's also depending on
  cloudiness)
%}
% select evaluation period
EP1 = 12:443;       % 10 March 1:00 to 28 March 0:00
EP2 = 468:923;      % 29 March 1:00 to 17 April 0:00

% select variables
LAI_Sodankyla = LAI;
PAI_Sodankyla = LAI + SAI;
Hourly_Time_Sodankyla = horzcat(time_1h(EP1(1):EP1(end)),time_1h(EP2(1):EP2(end)));
Cosine_Solar_Angle_Sodankyla = horzcat(Cosine_Solar_Angle(EP1(1):EP1(end)),...
    Cosine_Solar_Angle(EP2(1):EP2(end)));
SWR_Incoming_Sodankyla = horzcat(SW_In_Above_Vegetation(EP1(1):EP1(end)),...
    SW_In_Above_Vegetation(EP2(1):EP2(end)));
LWR_Atmosphere_Sodankyla = horzcat(LW_In_Above_Vegetation(EP1(1):EP1(end)),...
    LW_In_Above_Vegetation(EP2(1):EP2(end)));

bla = 5.67*10^(-8)*Air_Temperature_Above_Vegetation.^4;
blub = LW_In_Above_Vegetation'./bla;
Sky_Emissivity_Sodankyla = vertcat(blub(EP1(1):EP1(end)),blub(EP2(1):EP2(end)));
clear bla blub

LWR_Subcanopy_Observations_Sodankyla = vertcat(LW_in_bc_C_1h(EP1(1):EP1(end),:),...
    LW_in_bc_C_1h(EP2(1):EP2(end),:));
LWR_Subcanopy_CLM_Sodankyla = vertcat(LW_in_bc_CLM(EP1(1):EP1(end),:),...
    LW_in_bc_CLM(EP2(1):EP2(end),:));


save('MLRdata_Sodankyla.mat','LAI_Sodankyla','PAI_Sodankyla','Hourly_Time_Sodankyla',...
    'Cosine_Solar_Angle_Sodankyla','SWR_Incoming_Sodankyla','LWR_Atmosphere_Sodankyla',...
    'Sky_Emissivity_Sodankyla','LWR_Subcanopy_Observations_Sodankyla','LWR_Subcanopy_CLM_Sodankyla')