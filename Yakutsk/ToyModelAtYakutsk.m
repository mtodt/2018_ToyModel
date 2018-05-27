tic
clear
close all
%{
The Toy Model is designed in a way to facilitate general usage, requiring
input of parameters for vegetation and ground as well as meteorological
forcing. This file gives an overview of how to drive the Toy Model.
%}
% run DataPrep_Yakutsk.m first
load ForcingData_ToyModel_Yakutsk.mat

%-------------------------------------------------------------------------%
%--------------------  Location & General Parameters  --------------------%
Longitude = 129+37/60+08/60/60;
Latitude = 62+15/60+18/60/60;
Elevation = 220;
Height_of_Hum_Obs_98 = 31.4;
Height_of_Temp_Obs = 31.4;
Height_of_Wind_Obs = 27;
%{
Radiation measurements were conducted at 32m and 1.8m, the latter likely
because of understory vegetation. Labelled 1.2m in 2000, though.
%}
dt = time_98_1h(2)-time_98_1h(1);
Timestep_CLM = dt*24*3600;  % for CLM timestep required in [s]
Timestep_SP = dt*24;        % for Alptal Precipitation in mm/hour and timestep 1 hour -> dt = 1


%-------------------------------------------------------------------------%
%------------------------  Vegetation Parameters  ------------------------%
%{
List of PFTs used in CLM4.5 (abbreviations for boreal forests):
1 - not vegetated
2 - needleleaf evergreen temperate tree (NETTs)
3 - needleleaf evergreen boreal tree (NEBTs)
4 - needleleaf deciduous boreal tree (NDBTs)
5 - broadleaf evergreen tropical tree
6 - broadleaf evergreen temperate tree
7 - broadleaf deciduous tropical tree
8 - broadleaf deciduous temperate tree (BDTTs)
9 - broadleaf deciduous boreal tree (BDBTs)
10 - broadleaf evergreen shrub
11 - broadleaf deciduous temperate shrub
12 - broadleaf deciduous boreal shrub (BDBSs)
13 - c3 arctic grass (C3AGs)
14 - c3 non-arctic grass
15 - c4 grass
16 - c3 crop
17 - c3 irrigated
The following PFTs are optional, requiring activating the crop module.
18 - corn
19 - spring temperate cereal
20 - winter temperate cereal
21 - soybean
Note that Fortran starts counting with 0, so that the PFT indices used
within the CLM model code are 0 to 20!
%}
PFT = 4;    % larch -> needleleaf deciduous boreal tree
Vegetation_Height = 18;
%{
Ohta et al. (2001):
"The value of the plant area index (PAI) varied between 3.71 in the foliated
season and 1.71 in the leafless season. The PAI was obtained from analysis
of fish-eye photographs and was confirmed by litter fall observations.
...
After snowmelt, while there was still no foliage..."
%}
LAImin = 0;
LAImax = 2;
SAI = 1.71;
LAI = 0;

%-----------------------  only for SNOWPACK-2LHM  ------------------------%
Stand_Density = 840;    % [trees ha^{-1}]
Tree_Diameter = 0.2558; % [m], linear regression from trees used for sap flow
% basal area calculated [m^2 m^{-2}]
Basal_Area_of_Forest_Stand = pi*0.25*Tree_Diameter^2 * Stand_Density/100/100;

% calibrated for Alptal (?) and kept for other sites
Fraction_LAI_Top = 0.5;         % fraction of LAI attributed to uppermost layer
K_LAI = 0.75;               % radiation transmissivity parameter [0.4-0.8]
Fraction_Trunk_Height = 1;    % fraction of tree height occupied by trunks
    % '-> 0.2 for Alptal, 1 since leafless
    

%-------------------------------------------------------------------------%
%--------------------------  Ground Parameters  --------------------------%
% vegetation understory that determines albedo!
%{
Ohta et al. (2001):
"There is an evergreen broad-leaved Vaccinium understory vegetation 0.05m
high. Its leaf density is high, but the PAI (or LAI) was not measured.
...
Moreover, the albedo above the canopy was not high (0.22–0.27), even when
the larch stand had no foliage and the forest floor was snow-covered.
After snowmelt, while there was still no foliage, the albedo was
0.11–0.13."
However, that albedo might be above-canopy, as calculation from sub-canopy
SW radiation yields slightly higher albedo values of about 0.19 in summer.
%}
Soil_Albedo_Direct = 0.19;
Soil_Albedo_Diffuse = 0.19;
Sand_Fraction = 0.408048;
Clay_Fraction = 0.261906;
Fraction_Organic_Matter = 0.0934;
Fresh_Snow_Albedo = 0.8;


%-------------------------------------------------------------------------%
%--------------------------------  Flags  --------------------------------%
% no near-infrared measurements
NIR_available = 0;
%{
Soil moisture or water content measurements aren't always available so that
a proxy needs to be used. This proxy is then not necessarily a partial
volume but usually just a (unitless) fraction, which has to be treated
differently - saturation water content is scaled by the fraction to get
water content. Flag values are either 'proxy', if a unitless frction, or
'observation' if a partial volume.
%}
SoilWaterFlag = 'observation';
%{
CLM4.5 features a multiple-layer canopy but only for photosynthesis-related
calculations. This influences the energy balance only via the stomatal
resistance. By setting the canopy layers to 1, the "big leaf"
parameterisation is applied, which should be equal to the previous CLM4
version.
%}
CLM45_canopy_layers = 1;
%{
SNOWPACK's biomass parameterisation is transferable to CLM4.5, which is
enabled by setting HM_in_CLM to 'yes'.
%}
HM_in_CLM = 'no';
%{
Permitting insolation to reach the trunk level in SNOWPACK is optional,
decided by 'yes' or 'no'.
%}
SW_to_Trunk = 'yes';


%-------------------------------------------------------------------------%
%--------------------  Meteorological & Other Forcing  -------------------%
%-------------------------------------------------------------------------%
% calculate diffuse fraction from emissivity of the sky (same height!)
em_sky_98 = nan(size(T_air_98_1h_filled));
for t=1:length(time_98_1h)
    em_sky_98(t) = LW_in_ac_98_1h_filled(t)/(5.67*10^(-8)*T_air_98_1h_filled(t)^4);
end
frac_diff_98 = nan(size(em_sky_98));
for t=1:length(time_98_1h)
    frac_diff_98(t) = (em_sky_98(t)-min(em_sky_98))...
        /(max(em_sky_98)-min(em_sky_98));
end

% meteorology
SW_In_Above_Vegetation = SW_in_ac_98_1h_filled;
Diffuse_Fraction = frac_diff_98;
LW_In_Above_Vegetation = LW_in_ac_98_1h_filled;
Air_Temperature_Above_Vegetation = T_air_98_1h_filled;
Wind_Speed_Above_Vegetation = Wind_98_1h_filled;
Relative_Humidity_Above_Vegetation = RH_98_1h_filled;
% Conversion from mm/h to mm/s as required for CLM4.5
Precipitation = Precip_98_1h/3600;

% ground
Soil_Moisture = SM_98_1h_filled_avg/100;    % conversion from %

% determined by eye from top soil temperature
Snow_Depth_Below_Vegetation = nan(size(T_topsoil_98_1h_filled));
    Snow_Depth_Below_Vegetation(1:3195) = 0.2;
    Snow_Depth_Below_Vegetation(3196:end) = 0;
Surface_Temperature = T_surf_98_1h_filled;
Fraction_Frozen_Soil = frac_soil_frozen_98_1h;


%-------------------------------------------------------------------------%
%---------------------------  Initializations  ---------------------------%
Snowfall = zeros(size(Precipitation));
Rainfall = zeros(size(Precipitation));
Snow_Age = zeros(size(Snow_Depth_Below_Vegetation));
alb_surf_98 = nan(size(time_98_1h));
JulDay = nan(size(time_98_1h));

% vegetation hydrology
CanInt_98 = nan(size(time_98_1h));        % canopy interception [mm]
CanInt_98(1049) = 0;
CanWetFrac_98 = nan(size(time_98_1h));    % fraction of wet canopy
CanWetFrac_98(1049) = 0;

% vegetation temperature
T_veg_CLM_98 = nan(size(time_98_1h));        % CLM4.5 vegetation temperature [K]
T_veg_CLM_98(1049) = Air_Temperature_Above_Vegetation(1049) ...
    + 0.0098*(Height_of_Temp_Obs-Vegetation_Height);
T_leaf_SP_98 = nan(size(time_98_1h));        % SP-2LHM leaf layer temperature [K]
T_leaf_SP_98(1049) = Air_Temperature_Above_Vegetation(1049) ...
    + 0.0098*(Height_of_Temp_Obs-Vegetation_Height);
T_trunk_SP_98 = nan(size(time_98_1h));       % SP-2LHM trunk layer temperature [K]
T_trunk_SP_98(1049) = Surface_Temperature(1049);

% output variables
Cosine_Solar_Angle_98 = nan(size(time_98_1h));
SW_net_veg_CLM_98 = nan(size(time_98_1h));
SW_in_bc_dir_CLM_98 = nan(length(time_98_1h),2);
SW_in_bc_dif_CLM_98 = nan(length(time_98_1h),2);
LW_in_bc_CLM_98 = nan(size(time_98_1h));
LW_out_ac_CLM_98 = nan(size(time_98_1h));
T_air_wc_CLM_98 = nan(size(time_98_1h));
T_ref2m_CLM_98 = nan(size(time_98_1h));
EB_CLM_98 = nan(length(time_98_1h),17,41);
EB_SP_leaf_98 = nan(length(time_98_1h),16,7);
EB_SP_trunk_98 = nan(length(time_98_1h),16,7);
SW_net_can_SP_98 = nan(size(time_98_1h));
SW_net_trunk_SP_98 = nan(size(time_98_1h));
LW_in_bc_SP_98 = nan(size(time_98_1h));
LW_out_ac_SP_98 = nan(size(time_98_1h));
%{
Terms of Energy Balance matrix.
CLM4.5:
1)  net direct SW radiation
2)  net diffuse SW radiation
    1) + 2) = net SW radiation
3)  net LW radiation from atmosphere (gain) -> emv*forc_lwrad + (1-emv)*(1-emg)*emv*forc_lwrad
4)  net LW radiation from vegetation (loss) -> -2*emv*sb*TV_old^4 + emv*(1-emg)*emv*sb*TV_old^4
5)  net LW radiation from second vegetation layer (gain) -> 0 for CLM
6)  net LW radiation from ground (gain)     -> emv*emg*sb*Tgrnd^4
7)  net interaction with second vegetation layer -> 0 for CLM
    3) + 4) + [5) - 7)] + 6) = net LW radiation
    1) + 2) + 3) + 4) + [5) - 7)] + 6) + 7) = net radiation
    '-> complicated because of SNOWPACK
8)  net conductive heat flux -> 0 for CLM
9)  net sensible heat flux
10) net latent heat flux
11) dSW / dTveg -> per definition 0
12) dLW / dTveg -> derived from 4)
13) dTT / dTveg -> interaction with second vegetation layer 0 for CLM
14) dHM / dTveg -> conductive heat flux 0 for CLM
15) dSH / dTveg
16) dLH / dTveg
SNOWPACK-2LHM:
1)  net direct SW radiation
2)  net diffuse SW radiation
    1) + 2) = net SW radiation
3)  net LW radiation from atmosphere (gain)
4)  net LW radiation from vegetation (loss) -> emitted by respective layer
5)  net LW radiation from second vegetation layer (gain) -> absorbed from second layer
6)  net LW radiation from ground (gain)
7)  net interaction with second vegetation layer -> due to theory and
    derivation 7) already included in 5), therefore this:
    3) + 4) + [5) - 7)] + 6) = net LW radiation
    1) + 2) + 3) + 4) + [5) - 7)] + 6) + 7) = net radiation
8)  net conductive heat flux
9)  net sensible heat flux
10) net latent heat flux
11) dSW / dTveg -> per definition 0
12) dLW / dTveg -> derived from 4)
13) dTT / dTveg -> 0 for trunk layer because EB derived from this term for leaf layer
14) dHM / dTveg
15) dSH / dTveg
16) dLH / dTveg
%}


%-------------------------------------------------------------------------%
%----------------------------  Toy Model run  ----------------------------%
for t=1050:gap_98_1h(1)-1
% Julian Day
    date = datestr(time_98_1h(t),'yyyy-mm-dd-HH-MM-SS');
    Jultemp = datenum(str2double(date(1:4)),str2double(date(6:7)),str2double(date(9:10)),...
        str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));
    JulDay(t) = Jultemp-datenum(str2double(date(1:4))-1,12,31,0,0,0);
    
% direct-diffuse partitioning
    SW_In_Diffuse = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
    SW_In_Direct = SW_In_Above_Vegetation(t)-SW_In_Diffuse;
    
% rain-snow partitioning
%{
Essery et al. (2016) used threshold temperature of 2C and scaling factors
to determine snowfall. This is applied here as well. However, their air
temperature and precipitation measurements were conducted at the same
height, which is why we adjust threshold temperature using CLM4.5's lapse
rate.
Another solution would be to apply the algortihm used for Hyytiala during
SnowMIP2.
%}
    Threshold_Temperature = Air_Temperature_Above_Vegetation(t) + 0.0098*Height_of_Temp_Obs;
    if  Threshold_Temperature < 273.15+2
        Snowfall(t) = Precipitation(t);
    else
        Rainfall(t) = Precipitation(t);
    end
    
    if Snow_Depth_Below_Vegetation(t) > 0
        Snow_Fraction = 1;
    end
    
    if Snow_Depth_Below_Vegetation(t) > 0 && Snowfall(t) == 0
        Snow_Age(t) = Snow_Age(t-1) + dt;
    end
    
% surface albedo (constant snow cover)
    alb_surf_98(t) = (1-Snow_Fraction)*Soil_Albedo_Direct + Snow_Fraction*...
    ((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3);

% call Driver
    [Cosine_Solar_Angle_98(t),CanInt_98(t),CanWetFrac_98(t),SW_net_veg_CLM_98(t),LW_in_bc_CLM_98(t),...
        LW_out_ac_CLM_98(t),T_veg_CLM_98(t),T_air_wc_CLM_98(t),T_ref2m_CLM_98(t),...
        EB_CLM_98(t,:,:),EB_SP_leaf_98(t,:,:),EB_SP_trunk_98(t,:,:),...
        SW_net_can_SP_98(t),SW_net_trunk_SP_98(t),LW_in_bc_SP_98(t),LW_out_ac_SP_98(t),...
        T_leaf_SP_98(t),T_trunk_SP_98(t)] = ...
    SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay(t),...
        Timestep_CLM,Timestep_SP,PFT,LAI,SAI,Vegetation_Height,...
        Basal_Area_of_Forest_Stand,Tree_Diameter,Fraction_LAI_Top,K_LAI,...
        Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs_98,...
        Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
        SW_In_Direct,SW_In_Diffuse,NIR_available,LW_In_Above_Vegetation(t),...
        Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
        Relative_Humidity_Above_Vegetation(t),Rainfall(t),Snowfall(t),...
        Snow_Depth_Below_Vegetation(t),Snow_Fraction,Fraction_Frozen_Soil(t),alb_surf_98(t),...
        Surface_Temperature(t),T_veg_CLM_98(t-1),T_leaf_SP_98(t-1),T_trunk_SP_98(t-1),...
        Sand_Fraction,Clay_Fraction,Fraction_Organic_Matter,Soil_Moisture(t),...
        SoilWaterFlag,CanWetFrac_98(t-1),CanInt_98(t-1));
end

% continue after gap
CanInt_98(gap_98_1h(end)) = CanInt_98(gap_98_1h(1)-1);
CanWetFrac_98(gap_98_1h(end)) = CanWetFrac_98(gap_98_1h(1)-1);
T_veg_CLM_98(gap_98_1h(end)) = T_veg_CLM_98(gap_98_1h(1)-1);
T_leaf_SP_98(gap_98_1h(end)) = T_leaf_SP_98(gap_98_1h(1)-1);
T_trunk_SP_98(gap_98_1h(end)) = T_trunk_SP_98(gap_98_1h(1)-1);

for t=gap_98_1h(end)+1:length(time_98_1h)
% Julian Day
    date = datestr(time_98_1h(t),'yyyy-mm-dd-HH-MM-SS');
    JulDay(t) = datenum(0000,str2double(date(6:7)),str2double(date(9:10)),...
        str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));
    
% direct-diffuse partitioning
    SW_In_Diffuse = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
    SW_In_Direct = SW_In_Above_Vegetation(t)-SW_In_Diffuse;
    
% rain-snow partitioning
%{
Essery et al. (2016) used threshold temperature of 2C and scaling factors
to determine snowfall. This is applied here as well. However, their air
temperature and precipitation measurements were conducted at the same
height, which is why we adjust threshold temperature using CLM4.5's lapse
rate.
Another solution would be to apply the algortihm used for Hyytiala during
SnowMIP2.
%}
    Threshold_Temperature = Air_Temperature_Above_Vegetation(t) + 0.0098*Height_of_Temp_Obs;
    if  Threshold_Temperature < 273.15+2
        Snowfall(t) = Precipitation(t);
    else
        Rainfall(t) = Precipitation(t);
    end
    
    if Snow_Depth_Below_Vegetation(t) > 0
        Snow_Fraction = 1;
    end
    
    if Snow_Depth_Below_Vegetation(t) > 0 && Snowfall(t) == 0
        Snow_Age(t) = Snow_Age(t-1) + dt;
    end
    
% surface albedo (constant snow cover)
    alb_surf_98(t) = (1-Snow_Fraction)*Soil_Albedo_Direct + Snow_Fraction*...
    ((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3);

% call Driver
    [Cosine_Solar_Angle_98(t),CanInt_98(t),CanWetFrac_98(t),SW_net_veg_CLM_98(t),LW_in_bc_CLM_98(t),...
        LW_out_ac_CLM_98(t),T_veg_CLM_98(t),T_air_wc_CLM_98(t),T_ref2m_CLM_98(t),...
        EB_CLM_98(t,:,:),EB_SP_leaf_98(t,:,:),EB_SP_trunk_98(t,:,:),...
        SW_net_can_SP_98(t),SW_net_trunk_SP_98(t),LW_in_bc_SP_98(t),LW_out_ac_SP_98(t),...
        T_leaf_SP_98(t),T_trunk_SP_98(t)] = ...
    SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay(t),...
        Timestep_CLM,Timestep_SP,PFT,LAI,SAI,Vegetation_Height,...
        Basal_Area_of_Forest_Stand,Tree_Diameter,Fraction_LAI_Top,K_LAI,...
        Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs_98,...
        Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
        SW_In_Direct,SW_In_Diffuse,NIR_available,LW_In_Above_Vegetation(t),...
        Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
        Relative_Humidity_Above_Vegetation(t),Rainfall(t),Snowfall(t),...
        Snow_Depth_Below_Vegetation(t),Snow_Fraction,Fraction_Frozen_Soil(t),alb_surf_98(t),...
        Surface_Temperature(t),T_veg_CLM_98(t-1),T_leaf_SP_98(t-1),T_trunk_SP_98(t-1),...
        Sand_Fraction,Clay_Fraction,Fraction_Organic_Matter,Soil_Moisture(t),...
        SoilWaterFlag,CanWetFrac_98(t-1),CanInt_98(t-1));
end

toc