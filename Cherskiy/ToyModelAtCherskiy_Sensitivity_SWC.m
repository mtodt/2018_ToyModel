tic
clear
close all
%{
The Toy Model is designed in a way to facilitate general usage, requiring
input of parameters for vegetation and ground as well as meteorological
forcing. This file gives an overview of how to drive the Toy Model.
%}
% run DataPrep_Cherskiy.m first
load ForcingData_ToyModel_Cherskiy.mat

%-------------------------------------------------------------------------%
%--------------------  Location & General Parameters  --------------------%
%-------------------------------------------------------------------------%
Longitude = 161.45482;
Latitude = 68.75574;
Elevation = 39;
Height_of_Hum_Obs = 7;
Height_of_Temp_Obs = 7;
Height_of_Wind_Obs = 11+5;  % sensor height in open plus vegetation height
dt = time_1h(2)-time_1h(1);
Timestep_CLM = dt*24*3600;  % for CLM timestep required in [s]
Timestep_SP = dt*24;        % for Alptal Precipitation in mm/hour and timestep 1 hour -> dt = 1
JulDay = JulDay_1h;


%-------------------------------------------------------------------------%
%------------------------  Vegetation Parameters  ------------------------%
%-------------------------------------------------------------------------%
PFT = 4;                % larch -> needleleaf deciduous
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
Vegetation_Height = 5;
%{
There is only one LAI measurement, in late July - early August. To
calculate SAI and changes in LAI throughout the snowmelt season, we use
monthly LAI and SAI values from CLM's surface dataset.
Unfortunately, this yields unrealistic model results, which clearly
indicate that SAI is too low (general bias in model results). However,
we can calculate SAI geometrically assuming trees are cones. This
underestimates surface area due to not considering branches, but not all of 
the trees might be 5m tall. These simplifications should then cancel each 
other out to some degree.
For those calculations we need vegetation height, breast height, diameter
at breast height and number of trees per m^2.
%}
LAI = 0;    % observation on site
slant = sqrt(Vegetation_Height^2 * (1 + 0.017^2/(4*(Vegetation_Height-1.3)^2))); 
LSA = slant*pi*0.017*Vegetation_Height/(2*(Vegetation_Height-1.3));
SAI = 3.7*LSA;                      % number of trees per m^2, see below

%-----------------------  only for SNOWPACK-2LHM  ------------------------%
Tree_Diameter = 0.017;                % diameter at breast height [m]
Tree_Basal_Area = 12.9/100/100;       % basal area, conversion from cm
Trees_per_m2 = 3.7;
Basal_Area_of_Forest_Stand = Tree_Basal_Area*Trees_per_m2;

% calibrated for Alptal and kept here
Fraction_LAI_Top = 0.5;         % fraction of LAI attributed to uppermost layer
K_LAI = 0.75;                   % radiation transmissivity parameter [0.4-0.8]
Fraction_Trunk_Height = 0.2;    % fraction of tree height occupied by trunks


%-------------------------------------------------------------------------%
%--------------------------  Ground Parameters  --------------------------%
%-------------------------------------------------------------------------%
% estimated soil albedo from surface albedo measurements
Soil_Albedo_Direct = 0.0875;    % value from sub-canopy SWR measurements
Soil_Albedo_Diffuse = Soil_Albedo_Direct;
Clay_Fraction = 0.186888;
Sand_Fraction = 0.458740;
Fraction_Organic_Matter = 0.1230;
% snow albedo for surface albedo parameterization
Fresh_Snow_Albedo = 0.45;        % value from sub-canopy SWR measurements


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

%-----------------------------  meteorology  -----------------------------%
% calculate diffuse fraction from emissivity of the sky (same height!)
boltz = 5.67*10^(-8);
em_sky = nan(size(T_air_1h));
for t=1:length(T_air_1h)
    em_sky(t) = LW_in_ac_1h(t)/(boltz*T_air_1h(t)^4);
end
frac_diff = (em_sky-min(em_sky))/(max(em_sky)-min(em_sky));

SW_In_Above_Vegetation = SW_in_ac_1h;
Diffuse_Fraction = frac_diff;
LW_In_Above_Vegetation = LW_in_ac_1h;
Air_Temperature_Above_Vegetation = T_air_1h;
Wind_Speed_Above_Vegetation = Wind_1h;
Relative_Humidity_Above_Vegetation = RH_1h;
% Conversion from mm/h to mm/s as required for CLM4.5
Precipitation = Precip_1h/3600;

%-------------------------------  surface  -------------------------------%
%{
Soil is frozen most of the time, especially since we're not looking at
snow-free periods. So we assume soil is frozen completely.
Of course, snowmelt would have an impact on that.
%}
Fraction_Frozen_Soil = 1;

%{
Soil water content measurements are only available in late June and early
July 2017, so far-removed from snow cover and with soil thawed to some
degree. Multiple measurements are available per day and we average them for
a spatial representation. The lowest of these values, after several days
without precipitation, are at or less than 0.3 while individual lowest
values are between 0.15 and 0.2 (apart from one being 0.051).
Since we assume soil is consistently and completely frozen, soil moisture
has less of an impact, and we assume a low (considering post-snowmelt 
period) value.
%}
Soil_Water = nan(11,1);
for j=1:11
    Soil_Water(j) = min(1,max(0,(0.8+(j-1)*0.04)*0.15));
end
% Soil_Water = 0.15;

% snow fraction determined by eye from surface albedo measurements
Snow_Fraction = nan(size(time_1h));
Snow_Fraction(1:1187) = 1;
for s=1188:length(time_1h)
    Snow_Fraction(s) = Snow_Fraction(s-1) - 1/(length(time_1h)-1187);
end

% determine snow depth from surface albedo (since only presence necessary)
%{
There are no snow depth measurements but snow depth itself doesn't play a
major role. For SNOWPACK, 0.03m is a critical value for roughness length
(and 0 for choice of surface temperature). For CLM4.5, 0.05m is a critical
value for litter layer coverage by snow. That's already it, then, because
roughness length in CLM4.5 is determined by snow cover fraction, where 0 is
the critical value.
%}
Snow_Depth_Below_Vegetation = nan(size(Snow_Fraction));
for s=1:length(Snow_Depth_Below_Vegetation)
    if Snow_Fraction(s) > 0
        Snow_Depth_Below_Vegetation(s) = 0.2;
    end
end

% calculate surface temperature from outgoing LWR
Surface_Temperature = nan(size(LW_out_bc_1h));
for l=1:length(Surface_Temperature)
    em_gr = Snow_Fraction(l)*0.97 + (1-Snow_Fraction(l))*0.96;
    Surface_Temperature(l) = nthroot(LW_out_bc_1h(l)/(em_gr*boltz),4);
end


%-------------------------------------------------------------------------%
%---------------------------  Initializations  ---------------------------%
%-------------------------------------------------------------------------%
Rainfall = zeros(size(time_1h));
Snowfall = zeros(size(time_1h));
Snow_Age = zeros(size(Snow_Depth_Below_Vegetation));
Surface_Albedo = nan(size(time_1h));

% vegetation hydrology
CanInt = nan(size(time_1h));        % canopy interception [mm]
CanInt(1) = 0;
CanWetFrac = nan(size(time_1h));    % fraction of wet canopy
CanWetFrac(1) = 0;

% vegetation temperature
T_veg_CLM = nan(size(time_1h));        % CLM4.5 vegetation temperature [K]
T_veg_CLM(1) = Air_Temperature_Above_Vegetation(1);
T_leaf_SP = nan(size(time_1h));        % SP-2LHM leaf layer temperature [K]
T_leaf_SP(1) = Air_Temperature_Above_Vegetation(1);
T_trunk_SP = nan(size(time_1h));       % SP-2LHM trunk layer temperature [K]
T_trunk_SP(1) = Surface_Temperature(1);

% output variables
Cosine_Solar_Angle = nan(size(time_1h));
SW_net_veg_CLM = nan(size(time_1h));
LW_in_bc_CLM = nan(length(time_1h),11);
LW_out_ac_CLM = nan(size(time_1h));
T_air_wc_CLM = nan(size(time_1h));
T_ref2m_CLM = nan(size(time_1h));
EB_CLM = nan(length(time_1h),17,41);
EB_SP_leaf = nan(length(time_1h),16,7);
EB_SP_trunk = nan(length(time_1h),16,7);
SW_net_can_SP = nan(size(time_1h));
SW_net_trunk_SP = nan(size(time_1h));
LW_in_bc_SP = nan(length(time_1h),11);
LW_out_ac_SP = nan(size(time_1h));
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
for j=1:11
for t=2:length(time_1h)
% direct-diffuse partitioning
    SW_In_Diffuse(1) = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
    SW_In_Diffuse(2) = 0;    % no nir radiation
    SW_In_Direct(1) = SW_In_Above_Vegetation(t)*(1-Diffuse_Fraction(t));
    SW_In_Direct(2) = 0;    % no nir radiation

% rain-snow partitioning
%{
Essery et al. (2016) used threshold temperature of 2C and scaling factors
to determine snowfall. This is applied here as well.Another solution would 
be to apply the algortihm used for Hyytiala during SnowMIP2.
%}
    Threshold_Temperature = Air_Temperature_Above_Vegetation(t);
    if  Threshold_Temperature < 273.15+2
        Snowfall(t) = Precipitation(t);
    else
        Rainfall(t) = Precipitation(t);
    end

    if Snow_Depth_Below_Vegetation(t) > 0 && Snowfall(t) == 0
        Snow_Age(t) = Snow_Age(t-1) + dt;
    end
    
% surface albedo
    Surface_Albedo(t) = (1-Snow_Fraction(t)) * Soil_Albedo_Direct...
        + Snow_Fraction(t)*((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3);
        
% call Driver
    [Cosine_Solar_Angle(t),CanInt(t),CanWetFrac(t),SW_net_veg_CLM(t),LW_in_bc_CLM(t,j),...
        LW_out_ac_CLM(t),T_veg_CLM(t),T_air_wc_CLM(t),T_ref2m_CLM(t),...
        EB_CLM(t,:,:),EB_SP_leaf(t,:,:),EB_SP_trunk(t,:,:),...
        SW_net_can_SP(t),SW_net_trunk_SP(t),LW_in_bc_SP(t,j),LW_out_ac_SP(t),...
        T_leaf_SP(t),T_trunk_SP(t)] = ...
    SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay(t),...
        Timestep_CLM,Timestep_SP,PFT,LAI,SAI,Vegetation_Height,...
        Basal_Area_of_Forest_Stand,Tree_Diameter,Fraction_LAI_Top,K_LAI,...
        Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs,...
        Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
        SW_In_Direct,SW_In_Diffuse,NIR_available,LW_In_Above_Vegetation(t),...
        Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
        Relative_Humidity_Above_Vegetation(t),Rainfall(t),Snowfall(t),...
        Snow_Depth_Below_Vegetation(t),Snow_Fraction(t),Fraction_Frozen_Soil,...
        Surface_Albedo(t),Surface_Temperature(t),T_veg_CLM(t-1),...
        T_leaf_SP(t-1),T_trunk_SP(t-1),Sand_Fraction,Clay_Fraction,...
        Fraction_Organic_Matter,Soil_Water(j),SoilWaterFlag,CanWetFrac(t-1),CanInt(t-1));    
end
end

save('LWsub_Cherskiy_SWC.mat','LW_in_bc_1h','LW_in_bc_CLM','LW_in_bc_SP')

toc