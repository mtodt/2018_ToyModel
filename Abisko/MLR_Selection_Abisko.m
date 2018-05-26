%{
Collect ingredients necessary for multi-linear regression.
Emissivity of the sky & solar angle -> description of meteorological
conditions (insolation doesn't translate from simulation to simulation
because of differences in latitude and season, and it's also depending on
cloudiness).
%}
% select evaluation period
EP1 = 11:34;
EP2 = 107:130;
EP3 = 179:226;
EP4 = 395:442;
EP5 = 515:586;

% select variables
LAI_Abisko = LAI;
PAI_Abisko = LAI + SAI;
Hourly_Time_Abisko = horzcat(time_1h(EP1(1):EP1(end)),time_1h(EP2(1):EP2(end)),...
    time_1h(EP3(1):EP3(end)),time_1h(EP4(1):EP4(end)),time_1h(EP5(1):EP5(end)));
Cosine_Solar_Angle_Abisko = horzcat(Cosine_Solar_Angle(EP1(1):EP1(end)),Cosine_Solar_Angle(EP2(1):EP2(end)),...
    Cosine_Solar_Angle(EP3(1):EP3(end)),Cosine_Solar_Angle(EP4(1):EP4(end)),Cosine_Solar_Angle(EP5(1):EP5(end)));
SWR_Incoming_Abisko = horzcat(SW_In_Above_Vegetation(EP1(1):EP1(end)),...
    SW_In_Above_Vegetation(EP2(1):EP2(end)),SW_In_Above_Vegetation(EP3(1):EP3(end)),...
    SW_In_Above_Vegetation(EP4(1):EP4(end)),SW_In_Above_Vegetation(EP5(1):EP5(end)));
LWR_Atmosphere_Abisko = horzcat(LW_In_Above_Vegetation(EP1(1):EP1(end)),...
    LW_In_Above_Vegetation(EP2(1):EP2(end)),LW_In_Above_Vegetation(EP3(1):EP3(end)),...
    LW_In_Above_Vegetation(EP4(1):EP4(end)),LW_In_Above_Vegetation(EP5(1):EP5(end)));

bla = 5.67*10^(-8)*Air_Temperature_Above_Vegetation.^4;
blub = LW_In_Above_Vegetation./bla;
Sky_Emissivity_Abisko = horzcat(blub(EP1(1):EP1(end)),blub(EP2(1):EP2(end)),...
    blub(EP3(1):EP3(end)),blub(EP4(1):EP4(end)),blub(EP5(1):EP5(end)));
clear bla blub

LWR_Subcanopy_Observations_Abisko = vertcat(LW_in_bc_C_1h(EP1(1):EP1(end),:),LW_in_bc_C_1h(EP2(1):EP2(end),:),...
    LW_in_bc_C_1h(EP3(1):EP3(end),:),LW_in_bc_C_1h(EP4(1):EP4(end),:),LW_in_bc_C_1h(EP5(1):EP5(end),:));
LWR_Subcanopy_CLM_Abisko = vertcat(LW_in_bc_CLM(EP1(1):EP1(end),:),LW_in_bc_CLM(EP2(1):EP2(end),:),...
    LW_in_bc_CLM(EP3(1):EP3(end),:),LW_in_bc_CLM(EP4(1):EP4(end),:),LW_in_bc_CLM(EP5(1):EP5(end),:));

save('MLRdata_Abisko.mat','LAI_Abisko','PAI_Abisko','Hourly_Time_Abisko',...
    'Cosine_Solar_Angle_Abisko','SWR_Incoming_Abisko','LWR_Atmosphere_Abisko',...
    'Sky_Emissivity_Abisko','LWR_Subcanopy_Observations_Abisko','LWR_Subcanopy_CLM_Abisko')