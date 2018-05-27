%{
Collect ingredients necessary for multi-linear regression.
Emissivity of the sky & solar angle -> description of meteorological
conditions (insolation doesn't translate from simulation to simulation
because of differences in latitude and season, and it's also depending on
cloudiness).
%}
% select evaluation period
EP1 = 12:1163;
EP2 = 1188:length(time_1h);

% select variables
LAI_Cherskiy = LAI;
PAI_Cherskiy = LAI + SAI;
Hourly_Time_Cherskiy = horzcat(time_1h(EP1(1):EP1(end)),time_1h(EP2(1):EP2(end)));
Cosine_Solar_Angle_Cherskiy = horzcat(Cosine_Solar_Angle(EP1(1):EP1(end)),...
    Cosine_Solar_Angle(EP2(1):EP2(end)));
SWR_Incoming_Cherskiy = horzcat(SW_In_Above_Vegetation(EP1(1):EP1(end)),...
    SW_In_Above_Vegetation(EP2(1):EP2(end)));
LWR_Atmosphere_Cherskiy = horzcat(LW_In_Above_Vegetation(EP1(1):EP1(end)),...
    LW_In_Above_Vegetation(EP2(1):EP2(end)));

bla = 5.67*10^(-8)*Air_Temperature_Above_Vegetation.^4;
blub = LW_In_Above_Vegetation./bla;
Sky_Emissivity_Cherskiy = horzcat(blub(EP1(1):EP1(end)),blub(EP2(1):EP2(end)));
clear bla blub

LWR_Subcanopy_Observations_Cherskiy = horzcat(LW_in_bc_1h(EP1(1):EP1(end)),...
    LW_in_bc_1h(EP2(1):EP2(end)));
LWR_Subcanopy_CLM_Cherskiy = horzcat(LW_in_bc_CLM(EP1(1):EP1(end)),...
    LW_in_bc_CLM(EP2(1):EP2(end)));

save('MLRdata_Cherskiy.mat','LAI_Cherskiy','PAI_Cherskiy','Hourly_Time_Cherskiy',...
    'Cosine_Solar_Angle_Cherskiy','SWR_Incoming_Cherskiy','LWR_Atmosphere_Cherskiy',...
    'Sky_Emissivity_Cherskiy','LWR_Subcanopy_Observations_Cherskiy','LWR_Subcanopy_CLM_Cherskiy')