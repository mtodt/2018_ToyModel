tic
clear
close all
%{
The Toy Model is designed in a way to facilitate general usage, requiring
input of parameters for vegetation and ground as well as meteorological
forcing. This file gives an overview of how to drive the Toy Model.
%}
% run DataPrep_Alptal.m first
load ForcingData_ToyModel_Alptal.mat

%-------------------------------------------------------------------------%
%--------------------  Location & General Parameters  --------------------%
%-------------------------------------------------------------------------%
Longitude = 8+48/60;
Latitude = 47+3/60;
Elevation = 1220;   % open area on same altitude as measurements above forest
Height_of_Hum_Obs = 35;     % actually open site but same altitude as tower
Height_of_Temp_Obs = 35;    % actually open site but same altitude as tower
Height_of_Wind_Obs = 35;
dt = time_1h_all(2,2)-time_1h_all(1,2);
Timestep_CLM = dt*24*3600;  % for CLM timestep required in [s]
Timestep_SP = dt*24;        % for Alptal Precipitation in mm/hour and timestep 1 hour -> dt = 1


%-------------------------------------------------------------------------%
%------------------------  Vegetation Parameters  ------------------------%
%-------------------------------------------------------------------------%
PFT = 3;    % CLM4.5: 12.7% temperate, 9.8% boreal
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
Vegetation_Height = 25;
PAI = 4.1;  % average "LAI" along the rail
Fraction_LAI_Of_VAI = 0.79;   % average for Jan-Mar & both PFTs
LAI = Fraction_LAI_Of_VAI*PAI;
SAI = (1-Fraction_LAI_Of_VAI)*PAI;

%-----------------------  only for SNOWPACK-2LHM  ------------------------%
Basal_Area_of_Forest_Stand = 0.004; % [m2/m2], value for closed canopies like Alptal
Tree_Diameter = 1;                  % average canopy (tree) diameter [m]
% calibrated for Alptal and kept here
Fraction_LAI_Top = 0.5;         % fraction of LAI attributed to uppermost layer
K_LAI = 0.75;                   % radiation transmissivity parameter [0.4-0.8]
Fraction_Trunk_Height = 0.2;    % fraction of tree height occupied by trunks


%-------------------------------------------------------------------------%
%--------------------------  Ground Parameters  --------------------------%
%-------------------------------------------------------------------------%
Soil_Albedo_Direct = 0.11;
Soil_Albedo_Diffuse = 0.11;
%{
Since we're using a range of sites, for most of which there probably are no
soil composition information, it's easier to use data from CLM4.5 surface
dataset. Seehornwald corresponds to grid cell (33,584) in 0.31x0.23 surface
file. Averaged over whole soil column using layer thickness.
%}
Sand_Fraction = 0.484415;
Clay_Fraction = 0.239169;
Fraction_Organic_Matter = 0.0657;
% snow albedo for surface albedo parameterization
Fresh_Snow_Albedo = 0.8;


%-------------------------------------------------------------------------%
%--------------------------------  Flags  --------------------------------%
%-------------------------------------------------------------------------%
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
SoilWaterFlag = 'proxy';
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
% meteorology
SW_In_Above_Vegetation = SW_in_ac_1h_all;
Diffuse_Fraction = frac_diff_1h_all;
LW_In_Above_Vegetation = LW_in_ac_1h_all;
Air_Temperature_Above_Vegetation = T_air_ac_1h_all;
Wind_Speed_Above_Vegetation = Wind_1h_all;
Relative_Humidity_Above_Vegetation = RH_open_1h_all;
Precipitation = Precip_1h_all/3600; % conversion into mm s^{-1}

% ground
Soil_Water = SWProxy_1h_all;
Surface_Temperature = T_surf_1h_all;
Snow_Depth_Below_Vegetation = z_snow_forest_1h_all;
Snow_Fraction = frac_snow_1h_all;
%{
There is only one soil temperature for Alptal, at 20cm. This temperature
never drops below 274K during 2004-2007, so we assume soil is unfrozen
constantly.
%}
Fraction_Frozen_Soil = 0;


%-------------------------------------------------------------------------%
%---------------------------  Initializations  ---------------------------%
%-------------------------------------------------------------------------%
JulDay = nan(size(time_1h_all));
Snow_Age = zeros(size(Snow_Depth_Below_Vegetation));
Surface_Albedo = nan(size(Snow_Age));

% vegetation hydrology
CanInt = nan(size(time_1h_all));        % canopy interception [mm]
CanInt(start_2004,1) = 0;
CanInt(1,2:4) = 0;
CanWetFrac = nan(size(time_1h_all));    % fraction of wet canopy
CanWetFrac(start_2004,1) = 0;
CanWetFrac(1,2:4) = 0;

% vegetation temperature
T_veg_CLM = nan(size(time_1h_all));        % CLM4.5 vegetation temperature [K]
T_veg_CLM(start_2004,1) = Air_Temperature_Above_Vegetation(start_2004,1);
T_veg_CLM(1,2:4) = Air_Temperature_Above_Vegetation(1,2:4);
T_leaf_SP = nan(size(time_1h_all));        % SP-2LHM leaf layer temperature [K]
T_leaf_SP(start_2004,1) = Air_Temperature_Above_Vegetation(start_2004,1);
T_leaf_SP(1,2:4) = Air_Temperature_Above_Vegetation(1,2:4);
T_trunk_SP = nan(size(time_1h_all));       % SP-2LHM trunk layer temperature [K]
T_trunk_SP(start_2004,1) = Surface_Temperature(start_2004,1);
T_trunk_SP(1,2:4) = Surface_Temperature(1,2:4);

% output variables
Cosine_Solar_Angle = nan(size(time_1h_all));
SW_net_veg_CLM = nan(size(time_1h_all));
LW_in_bc_CLM = nan(size(time_1h_all));
LW_out_ac_CLM = nan(size(time_1h_all));
T_air_wc_CLM = nan(size(time_1h_all));
T_ref2m_CLM = nan(size(time_1h_all));
EB_CLM = nan(length(time_1h_all(:,1)),length(time_1h_all(1,:)),17,41);
EB_SP_leaf = nan(length(time_1h_all(:,1)),length(time_1h_all(1,:)),16,7);
EB_SP_trunk = nan(length(time_1h_all(:,1)),length(time_1h_all(1,:)),16,7);
SW_net_can_SP = nan(size(time_1h_all));
SW_net_trunk_SP = nan(size(time_1h_all));
LW_in_bc_SP = nan(size(time_1h_all));
LW_out_ac_SP = nan(size(time_1h_all));
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
%-------------------------------------------------------------------------%
for t=start_2004+1:length(time_1h_all(:,1))
    date = datestr(time_1h_all(t,1),'yyyy-mm-dd-HH-MM-SS');
    Jultemp = datenum(str2double(date(1:4)),str2double(date(6:7)),str2double(date(9:10)),...
        str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));
    JulDay(t,1) = Jultemp-datenum(str2double(date(1:4))-1,12,31,0,0,0);
    
%{
Determine fraction of direct & diffuse insolation. In the absence of
specific measurements we use the effective emissivity of the sky, which is
a proxy for cloudiness. We scale emissivity to a range of 0 to 1 using its
own minimum and maximum (over all 4 years).
Note that the impact of this fraction on the sub-canopy LW radiation is
marginal.
%}
    forc_sol_dir(1) = SW_In_Above_Vegetation(t,1)*(1-Diffuse_Fraction(t,1));
    forc_sol_dir(2) = 0;    % no nir radiation
    forc_sol_dif(1) = SW_In_Above_Vegetation(t,1)*Diffuse_Fraction(t,1);
    forc_sol_dif(2) = 0;    % no nir radiation
        
%{
Differentiating between snow and rain using the algorithm given by Rutter
et al. (2009).
%}
    frac_prec_snow = RainSnowPartitioningAlptal(Air_Temperature_Above_Vegetation(t,1));
    forc_snow = Precipitation(t,1)*frac_prec_snow;
    forc_rain = Precipitation(t,1)*(1-frac_prec_snow);
        
% snow age for surface albedo parameterization
    if Snow_Depth_Below_Vegetation(t,1) > 0 && forc_snow == 0
        Snow_Age(t,1) = Snow_Age(t-1,1) + dt;   % snow age in days for ageing
    elseif Snow_Depth_Below_Vegetation(t,1) < Snow_Depth_Below_Vegetation(t-1,1)
        Snow_Age(t,1) = Snow_Age(t-1,1) + dt;   % snow age in days for ageing
    end
    
% surface albedo parameterization via snow ageing
    Surface_Albedo(t,1) = (1-Snow_Fraction(t,1)) * Soil_Albedo_Direct...
        + Snow_Fraction(t,1)*((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t,1)/7) + 0.3);
    
[Cosine_Solar_Angle(t,1),CanInt(t,1),CanWetFrac(t,1),SW_net_veg_CLM(t,1),LW_in_bc_CLM(t,1),...
    LW_out_ac_CLM(t,1),T_veg_CLM(t,1),T_air_wc_CLM(t,1),T_ref2m_CLM(t,1),...
    EB_CLM(t,1,:,:),EB_SP_leaf(t,1,:,:),EB_SP_trunk(t,1,:,:),...
    SW_net_can_SP(t,1),SW_net_trunk_SP(t,1),LW_in_bc_SP(t,1),...
    LW_out_ac_SP(t,1),T_leaf_SP(t,1),T_trunk_SP(t,1)] = ...
SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay(t,1),...
    Timestep_CLM,Timestep_SP,PFT,LAI,SAI,Vegetation_Height,...
    Basal_Area_of_Forest_Stand,Tree_Diameter,Fraction_LAI_Top,K_LAI,...
    Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs,...
    Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
    forc_sol_dir,forc_sol_dif,NIR_available,LW_In_Above_Vegetation(t,1),...
    Air_Temperature_Above_Vegetation(t,1),Wind_Speed_Above_Vegetation(t,1),...
    Relative_Humidity_Above_Vegetation(t,1),forc_rain,forc_snow,...
    Snow_Depth_Below_Vegetation(t,1),Snow_Fraction(t,1),Fraction_Frozen_Soil,...
    Surface_Albedo(t,1),Surface_Temperature(t,1),T_veg_CLM(t-1,1),...
    T_leaf_SP(t-1,1),T_trunk_SP(t-1,1),Sand_Fraction,Clay_Fraction,...
    Fraction_Organic_Matter,Soil_Water(t,1),SoilWaterFlag,CanWetFrac(t-1,1),CanInt(t-1,1));
    
end
for i=2:4
    for t=2:length(time_1h_all(:,i))
        if isnan(time_1h_all(t,i))
            break
        end
        date = datestr(time_1h_all(t,i),'yyyy-mm-dd-HH-MM-SS');
        Jultemp = datenum(str2double(date(1:4)),str2double(date(6:7)),str2double(date(9:10)),...
            str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));
        JulDay(t,i) = Jultemp-datenum(str2double(date(1:4))-1,12,31,0,0,0);
        
%{
Determine fraction of direct & diffuse insolation. In the absence of
specific measurements we use the effective emissivity of the sky, which is
a proxy for cloudiness. We scale emissivity to a range of 0 to 1 using its
own minimum and maximum (over all 4 years).
Note that the impact of this fraction on the sub-canopy LW radiation is
marginal.
%}
        forc_sol_dir(1) = SW_In_Above_Vegetation(t,i)*(1-Diffuse_Fraction(t,i));
        forc_sol_dir(2) = 0;    % no nir radiation
        forc_sol_dif(1) = SW_In_Above_Vegetation(t,i)*Diffuse_Fraction(t,i);
        forc_sol_dif(2) = 0;    % no nir radiation
        
%{
Differentiating between snow and rain using the algorithm given by Rutter
et al. (2009).
%}
        frac_prec_snow = RainSnowPartitioningAlptal(Air_Temperature_Above_Vegetation(t,i));
        forc_snow = Precipitation(t,i)*frac_prec_snow;
        forc_rain = Precipitation(t,i)*(1-frac_prec_snow);
        
% snow age for surface albedo parameterization
        if Snow_Depth_Below_Vegetation(t,i) > 0 && forc_snow == 0
            Snow_Age(t,i) = Snow_Age(t-1,i) + dt;   % snow age in days for ageing
        elseif Snow_Depth_Below_Vegetation(t,i) < Snow_Depth_Below_Vegetation(t-1,i)
            Snow_Age(t,i) = Snow_Age(t-1,i) + dt;   % snow age in days for ageing
        end
        
% surface albedo parameterization via snow ageing
    Surface_Albedo(t,i) = (1-Snow_Fraction(t,i)) * Soil_Albedo_Direct...
        + Snow_Fraction(t,i)*((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t,i)/7) + 0.3);
        
    [Cosine_Solar_Angle(t,i),CanInt(t,i),CanWetFrac(t,i),SW_net_veg_CLM(t,i),LW_in_bc_CLM(t,i),...
        LW_out_ac_CLM(t,i),T_veg_CLM(t,i),T_air_wc_CLM(t,i),T_ref2m_CLM(t,i),...
        EB_CLM(t,i,:,:),EB_SP_leaf(t,i,:,:),EB_SP_trunk(t,i,:,:),...
        SW_net_can_SP(t,i),SW_net_trunk_SP(t,i),LW_in_bc_SP(t,i),...
        LW_out_ac_SP(t,i),T_leaf_SP(t,i),T_trunk_SP(t,i)] = ...
    SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay(t,i),...
        Timestep_CLM,Timestep_SP,PFT,LAI,SAI,Vegetation_Height,...
        Basal_Area_of_Forest_Stand,Tree_Diameter,Fraction_LAI_Top,K_LAI,...
        Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs,...
        Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
        forc_sol_dir,forc_sol_dif,NIR_available,LW_In_Above_Vegetation(t,i),...
        Air_Temperature_Above_Vegetation(t,i),Wind_Speed_Above_Vegetation(t,i),...
        Relative_Humidity_Above_Vegetation(t,i),forc_rain,forc_snow,...
        Snow_Depth_Below_Vegetation(t,i),Snow_Fraction(t,i),Fraction_Frozen_Soil,...
        Surface_Albedo(t,i),Surface_Temperature(t,i),T_veg_CLM(t-1,i),...
        T_leaf_SP(t-1,i),T_trunk_SP(t-1,i),Sand_Fraction,Clay_Fraction,...
        Fraction_Organic_Matter,Soil_Water(t,i),SoilWaterFlag,CanWetFrac(t-1,i),CanInt(t-1,i));
    end
end
toc