tic
clear
close all
%{
The Toy Model is designed in a way to facilitate general usage, requiring
input of parameters for vegetation and ground as well as meteorological
forcing. This file gives an overview of how to drive the Toy Model.
%}
% run DataPrep_Sodankyla.m first
load ForcingData_ToyModel_Sodankyla.mat

%-------------------------------------------------------------------------%
%--------------------  Location & General Parameters  --------------------%
Longitude = 26+38/60;
Latitude = 67+22/60;
Elevation = 179;
Height_of_Hum_Obs = 48;
Height_of_Temp_Obs = 48;
Height_of_Wind_Obs = 47;
%{
Wind, temperature and relative humidity available from met tower wit
highest level 47m (or 48m). Since locations of sub-canopy measuements
change, this level is chosen to use observations least affected by
vegetation.
Radiation also at top of met tower, but different one.
%}
dt = time_1h(2)-time_1h(1);
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
PFT = 3;

%{
from Nick's script: PAI from hemi photos using van Gardinen et al 1999 - 
see Reid et al. 2014 data were cut and pasted from Tim's spreadsheets 
LwDays 080714.xlsx (LW positions) and SVPAIgraphs.xlsx (SW positions).
%}
% PAI at each of the 4 LW positions at each site, veg from Reid et al. (2014)
VAI = [1.08648 1.22392 1.11532 1.14996];      % tall pine, high density

%{
Duration of measurements (in completed days) 39 days, however, Reid et al.
(2014) give different numbers after quality control is 23 days.
%}
Fraction_LAI_Of_VAI = 0.78; % value from CLM4.5
Vegetation_Height = 18;
LAI = Fraction_LAI_Of_VAI * VAI_C;
SAI = (1 - Fraction_LAI_Of_VAI) * VAI_C;

%-----------------------  only for SNOWPACK-2LHM  ------------------------%
Basal_Area_of_Forest_Stand = 0.002452;
Tree_Diameter = 0.1155;

% calibrated for Alptal (?) and kept for other sites
Fraction_LAI_Top = 0.5;         % fraction of LAI attributed to uppermost layer
K_LAI = 0.75;                   % radiation transmissivity parameter [0.4-0.8]
Fraction_Trunk_Height = 0.5;    % fraction of tree height occupied by trunks
    % '-> 0.2 for Alptal, but larger trunk fraction from photos


%-------------------------------------------------------------------------%
%--------------------------  Ground Parameters  --------------------------%
% constant snow cover therefore soil albedo redundant
Soil_Albedo_Direct = nan;
Soil_Albedo_Diffuse =  nan;
Snow_Fraction = 1;
Sand_Fraction = 0.607785;
Clay_Fraction = 0.121581;
Fraction_Organic_Matter = 0.2494;
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
% meteorology
SW_In_Above_Vegetation = SW_in_open_1h;
Diffuse_Fraction = f_SW_diff_1h;
LW_In_Above_Vegetation = LW_in_open_1h;
Air_Temperature_Above_Vegetation = T_air_mast;
Wind_Speed_Above_Vegetation = Wind_mast;
Relative_Humidity_Above_Vegetation = RH_mast;
%{
Essery et al. (2016) used optical precipitation measurements, which should
be the Vaisala FD12P instrument corresponding to Prec_met_present.
Therefore, we use this variable here.
Conversion from mm/h to mm/s as required for CLM4.5
%}
Precipitation = Prec_met_present/3600;

% ground
Soil_Water = SM_avg;
%{
Fig. 6 by Essery et al. (2016) shows snow depth in the forest being about
75% of that in the open. SInce the actual snow depth is not important for
the Toy Model, this crude adjustment is used.
%}
Snow_Depth_Below_Vegetation = z_snow_open*0.75;
Fraction_Frozen_Soil = frac_soil_frozen;
Surface_Temperature = T_surf_C;


%-------------------------------------------------------------------------%
%---------------------------  Initializations  ---------------------------%
Snowfall = zeros(size(Precipitation));
Rainfall = zeros(size(Precipitation));
Snow_Age = zeros(size(Snow_Depth_Below_Vegetation));
Surface_Albedo = nan(size(time_1h));
JulDay = nan(size(time_1h));
        
% vegetation hydrology
CanInt = nan(length(time_1h),4);        % canopy interception [mm]
CanInt(1,:) = 0;
CanWetFrac = nan(length(time_1h),4);    % fraction of wet canopy
CanWetFrac(1,:) = 0;

% vegetation temperature
T_veg_CLM = nan(length(time_1h),4);        % CLM4.5 vegetation temperature [K]
T_veg_CLM(1,:) = Air_Temperature_Above_Vegetation(1) ...
    + 0.0098*(Height_of_Temp_Obs-Vegetation_Height);
T_leaf_SP = nan(length(time_1h),4);        % SP-2LHM leaf layer temperature [K]
T_leaf_SP(1,:) = Air_Temperature_Above_Vegetation(1) ...
    + 0.0098*(Height_of_Temp_Obs-Vegetation_Height);
T_trunk_SP = nan(length(time_1h),4);       % SP-2LHM trunk layer temperature [K]
T_trunk_SP(1,:) = Surface_Temperature(1);

% output variables
Cosine_Solar_Angle = nan(size(time_1h));
SW_net_veg_CLM = nan(length(time_1h),4);
LW_in_bc_CLM = nan(length(time_1h),4);
LW_out_ac_CLM = nan(length(time_1h),4);
T_air_wc_CLM = nan(length(time_1h),4);
T_ref2m_CLM = nan(length(time_1h),4);
EB_CLM = nan(length(time_1h),4,17,41);
EB_SP_leaf = nan(length(time_1h),4,16,7);
EB_SP_trunk = nan(length(time_1h),4,16,7);
SW_net_can_SP = nan(length(time_1h),4);
SW_net_trunk_SP = nan(length(time_1h),4);
LW_in_bc_SP = nan(length(time_1h),4);
LW_out_ac_SP = nan(length(time_1h),4);
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
for i=1:4       % 4 different LW radiometers with different VAIs
    t_start = 1;
    for t=t_start+1:length(time_1h)
% Julian Day
        date = datestr(time_1h(t),'yyyy-mm-dd-HH-MM-SS');
        JulDay(t) = datenum(0000,str2double(date(6:7)),str2double(date(9:10)),...
            str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));
        
% direct-diffuse partitioning
        SW_In_Diffuse(1) = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
        SW_In_Diffuse(2) = 0;
        SW_In_Direct(1) = SW_In_Above_Vegetation(t)-SW_In_Diffuse(1);
        SW_In_Direct(2) = 0;
        
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
        Threshold_Temperature = Air_Temperature_Above_Vegetation(t) ...
            + 0.0098*Height_of_Temp_Obs;
        if  Threshold_Temperature < 273.15+2
            Snowfall(t) = Precipitation(t);
        else
            Rainfall(t) = Precipitation(t);
        end
        
        if Snow_Depth_Below_Vegetation(t) > 0 && Snowfall(t) == 0
            Snow_Age(t) = Snow_Age(t-1) + dt;
        end
    
% surface albedo (constant snow cover)
        Surface_Albedo(t) = (Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3;
        
% call Driver
        [Cosine_Solar_Angle(t),CanInt(t,i),CanWetFrac(t,i),SW_net_veg_CLM(t,i),LW_in_bc_CLM(t,i),...
            LW_out_ac_CLM(t,i),T_veg_CLM(t,i),T_air_wc_CLM(t,i),T_ref2m_CLM(t,i),...
            EB_CLM(t,i,:,:),EB_SP_leaf(t,i,:,:),EB_SP_trunk(t,i,:,:),...
            SW_net_can_SP(t,i),SW_net_trunk_SP(t,i),LW_in_bc_SP(t,i),LW_out_ac_SP(t,i),...
            T_leaf_SP(t,i),T_trunk_SP(t,i)] = ...
        SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay(t),...
            Timestep_CLM,Timestep_SP,PFT,LAI(i),SAI(i),Vegetation_Height,...
            Basal_Area_of_Forest_Stand,Tree_Diameter,Fraction_LAI_Top,K_LAI,...
            Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs,...
            Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
            SW_In_Direct,SW_In_Diffuse,NIR_available,LW_In_Above_Vegetation(t),...
            Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
            Relative_Humidity_Above_Vegetation(t),Rainfall(t),Snowfall(t),...
            Snow_Depth_Below_Vegetation(t),Snow_Fraction,Fraction_Frozen_Soil(t),Surface_Albedo(t),...
            Surface_Temperature(t),T_veg_CLM(t-1,i),T_leaf_SP(t-1,i),T_trunk_SP(t-1,i),...
            Sand_Fraction,Clay_Fraction,Fraction_Organic_Matter,Soil_Water(t),...
            SoilWaterFlag,CanWetFrac(t-1,i),CanInt(t-1,i));
    end
end
toc