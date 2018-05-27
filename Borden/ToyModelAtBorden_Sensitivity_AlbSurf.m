tic
clear
close all
%{
The Toy Model is designed in a way to facilitate general usage, requiring
input of parameters for vegetation and ground as well as meteorological
forcing. This file gives an overview of how to drive the Toy Model.
%}
% run DataPrep_Borden.m first
load ForcingData_ToyModel_Borden.mat


%-------------------------------------------------------------------------%
%--------------------  Location & General Parameters  --------------------%
%-------------------------------------------------------------------------%
% Borden Forest Research Station ON Canada
Longitude = -(79+56/60);  % should not matter whether negative or >180°
Latitude = 44+19/60;
Elevation = 222;    % altitude of military base next to the forest
%{
Relative humidity and air temperature are chosen from a different height so
that these are in accordance with radiation measurements.
%}
Height_of_Hum_Obs = 33;
Height_of_Temp_Obs = 33;
Height_of_Wind_Obs = 42.7;
dt = JulDay_1h(2)-JulDay_1h(1);
Timestep_CLM = dt*24*3600;  % for CLM timestep required in [s]
Timestep_SP = dt*24;        % for Alptal Precipitation in mm/hour and timestep 1 hour -> dt = 1
boltz = 5.67*10^(-8);


%-------------------------------------------------------------------------%
%------------------------  Vegetation Parameters  ------------------------%
%-------------------------------------------------------------------------%
PFT = [2 8]; % mixed forest in temperate part of ecotone
PFT_Fractions = [0.1867 0.8133];
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
Vegetation_Height = 22; % average (Teklemariam et al., 2009; Croft et al., (2015))

%-----------------  estimate LAI & SAI from observation  -----------------%
%{
"At DOY 141, LAI values were at 98% of their seasonal maximum..." (Croft et
al., 2015) - referring to 2013
mid-growing season: LAI = 4.6 -> stayed at that value until ~September,
which is in accordance with LAI fractions for that grid cell.
Private communication with Paul Bartlett: max PAI about 5.6.
%}
SAI = [0.48 1.1];
LAI = [1.93 0.05]; % 0 doesn't work for latent heat flux

%-----------------------  only for SNOWPACK-2LHM  ------------------------%
Stand_Density = 2996;       % [trees ha^{-1}] including understory >1.5m
living_tree_stem_density = 0.3; % [m^{-2}]
Tree_Diameter = [0.123 0.0682];    % based on Neumann et al. (1989)
Basal_Area_of_Forest_Stand(1) = pi*(Tree_Diameter(1)/2)^2 * living_tree_stem_density;
Basal_Area_of_Forest_Stand(2) = pi*(Tree_Diameter(2)/2)^2 * living_tree_stem_density;

% calibrated for Alptal (?) and kept for other sites
Fraction_LAI_Top = 0.5;         % fraction of LAI attributed to uppermost layer
K_LAI = 0.75;               % radiation transmissivity parameter [0.4-0.8]
Fraction_Trunk_Height = 0.2;    % fraction of tree height occupied by trunks


%-------------------------------------------------------------------------%
%--------------------------  Ground Parameters  --------------------------%
%-------------------------------------------------------------------------%
Soil_Albedo_Direct = 0.195;     % derived from histogram for surface albedo
Soil_Albedo_Diffuse = Soil_Albedo_Direct;
%{
blub1 = alb_surf_1h;
blub2 = alb_surf_all_1h;
for l=1:length(alb_surf_1h)
    if alb_surf_1h(l) > 1
        blub1(l) = nan;
    end
    if alb_surf_all_1h(l) > 1
        blub2(l) = nan;
    end
end
hist(blub1,100)
hist(blub2,100)
%}
Sand_Fraction = 0.710686;
Clay_Fraction = 0.029556;
Fraction_Organic_Matter = 0.057566;
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
em_sky = nan(size(T_air_ac_1h_filled));
for t=1:length(T_air_ac_1h_filled)
    em_sky(t) = LW_in_ac_1h_filled(t)/(boltz*T_air_ac_1h_filled(t)^4);
end
frac_diff = (em_sky-min(em_sky))/(max(em_sky)-min(em_sky));

SW_In_Above_Vegetation = SW_in_ac_1h_filled;
Diffuse_Fraction = frac_diff;
LW_In_Above_Vegetation = LW_in_ac_1h_filled;
Air_Temperature_Above_Vegetation = T_air_ac_1h_filled;
Wind_Speed_Above_Vegetation = Wind_1h_filled;
Relative_Humidity_Above_Vegetation = RH_1h_filled;
Precipitation = Precip_1h;  % already in mm s^{-1}

%-------------------------------  surface  -------------------------------%
Soil_Water = SoilMoist_filter_avg_1h_filled;
Fraction_Frozen_Soil = frac_soil_frozen_1h_filled;

% snow fraction determined by eye from surface albedo measurements
Snow_Fraction = nan(size(JulDay_1h));
Snow_Fraction(1:1900) = 1;  % maybe dip from 1500 - 1900 with peak at 1700
Snow_Fraction(2251:end) = 0;
for s=1901:2250
    Snow_Fraction(s) = Snow_Fraction(s-1) - 1/350;
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
Surface_Temperature = nan(size(LW_out_bc_1h_filled));
for l=1:length(LW_out_bc_1h_filled)
    em_gr = Snow_Fraction(l)*0.97 + (1-Snow_Fraction(l))*0.96;
    Surface_Temperature(l) = nthroot(LW_out_bc_1h_filled(l)/(em_gr*boltz),4);
end


%-------------------------------------------------------------------------%
%---------------------------  Initializations  ---------------------------%
%-------------------------------------------------------------------------%
Snow_Age = zeros(size(Snow_Depth_Below_Vegetation));
alb_surf_mod = nan(size(JulDay_1h));
Surface_Albedo = nan(size(JulDay_1h));

% vegetation hydrology
CanInt = nan(length(JulDay_1h),2);        % canopy interception [mm]
CanInt(10,:) = 0;
CanWetFrac = nan(length(JulDay_1h),2);    % fraction of wet canopy
CanWetFrac(10,:) = 0;

% vegetation temperature
T_veg_CLM = nan(length(JulDay_1h),2);        % CLM4.5 vegetation temperature [K]
T_veg_CLM(10,:) = Air_Temperature_Above_Vegetation(10) ...
    + 0.0098*(Height_of_Temp_Obs-Vegetation_Height);
T_leaf_SP = nan(length(JulDay_1h),2);        % SP-2LHM leaf layer temperature [K]
T_leaf_SP(10,:) = Air_Temperature_Above_Vegetation(10) ...
    + 0.0098*(Height_of_Temp_Obs-Vegetation_Height);
T_trunk_SP = nan(length(JulDay_1h),2);       % SP-2LHM trunk layer temperature [K]
T_trunk_SP(10,:) = Surface_Temperature(10);

% output variables
Cosine_Solar_Angle = nan(size(JulDay_1h));
SW_net_veg_CLM = nan(length(JulDay_1h),2);
LW_in_bc_CLM = nan(length(JulDay_1h),2,11);
LW_in_bc_CLM_comb = nan(length(JulDay_1h),11);
LW_out_ac_CLM = nan(length(JulDay_1h),2);
T_air_wc_CLM = nan(length(JulDay_1h),2);
T_ref2m_CLM = nan(length(JulDay_1h),2);
EB_CLM = nan(length(JulDay_1h),17,41,2);
EB_SP_leaf = nan(length(JulDay_1h),16,7,2);
EB_SP_trunk = nan(length(JulDay_1h),16,7,2);
SW_net_can_SP = nan(length(JulDay_1h),2);
SW_net_trunk_SP = nan(length(JulDay_1h),2);
LW_in_bc_SP = nan(length(JulDay_1h),2,11);
LW_out_ac_SP = nan(length(JulDay_1h),2);
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
for t=11:gap_1(1)-1
%{
Determine fraction of direct & diffuse insolation. In the absence of
specific measurements we use the effective emissivity of the sky, which is
a proxy for cloudiness. We scale emissivity to a range of 0 to 1 using its
own minimum and maximum (over all 4 years).
Note that the impact of this fraction on the sub-canopy LW radiation is
marginal.
%}
    forc_sol_dir(1) = SW_In_Above_Vegetation(t)*(1-Diffuse_Fraction(t));
    forc_sol_dir(2) = 0;    % no nir radiation
    forc_sol_dif(1) = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
    forc_sol_dif(2) = 0;    % no nir radiation
    
%{
Essery et al. (2016) used threshold temperature of 2C and scaling factors
to determine snowfall. This is applied here as well. Air temperature
measurements were conducted at a tower but height is uncertain for
precipitation measurements - but likely on ground level, which is why we
adjust threshold temperature using CLM4.5's lapse rate.
Another solution would be to apply an algorithm used for SnowMIP2.
%}
    Threshold_Temperature = Air_Temperature_Above_Vegetation(t) + 0.0098*Height_of_Temp_Obs;
    if  Threshold_Temperature < 273.15+2
        forc_snow = Precipitation(t);
        forc_rain = 0;
    else
        forc_snow = 0;
        forc_rain = Precipitation(t);
    end
        
% snow age for surface albedo parameterization
    if Snow_Depth_Below_Vegetation(t) > 0 && forc_snow == 0
        Snow_Age(t) = Snow_Age(t-1) + dt;   % snow age in days for ageing
    else
        Snow_Age(t) = Snow_Age(t-1) + dt;   % snow age in days for ageing
    end
% surface albedo parameterization via snow ageing
    Surface_Albedo(t) = (1-Snow_Fraction(t)) * Soil_Albedo_Direct...
        + Snow_Fraction(t)*((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3);
    Surface_Albedo(t) = min(1,max(0,(0.8+(j-1)*0.04)*Surface_Albedo(t)));
    
    for i=1:2
        [Cosine_Solar_Angle(t),CanInt(t,i),CanWetFrac(t,i),...
            SW_net_veg_CLM(t,i),LW_in_bc_CLM(t,i,j),LW_out_ac_CLM(t,i),...
            T_veg_CLM(t,i),T_air_wc_CLM(t,i),T_ref2m_CLM(t,i),...
            EB_CLM(t,:,:,i),EB_SP_leaf(t,:,:,i),EB_SP_trunk(t,:,:,i),...
            SW_net_can_SP(t,i),SW_net_trunk_SP(t,i),LW_in_bc_SP(t,i,j),...
            LW_out_ac_SP(t,i),T_leaf_SP(t,i),T_trunk_SP(t,i)] = ...
        SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay_1h(t),...
            Timestep_CLM,Timestep_SP,PFT(i),LAI(i),SAI(i),Vegetation_Height,...
            Basal_Area_of_Forest_Stand(i),Tree_Diameter(i),Fraction_LAI_Top,K_LAI,...
            Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs,...
            Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
            forc_sol_dir,forc_sol_dif,NIR_available,LW_In_Above_Vegetation(t),...
            Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
            Relative_Humidity_Above_Vegetation(t),forc_rain,forc_snow,...
            Snow_Depth_Below_Vegetation(t),Snow_Fraction(t),Fraction_Frozen_Soil(t),...
            Surface_Albedo(t),Surface_Temperature(t),T_veg_CLM(t-1,i),...
            T_leaf_SP(t-1,i),T_trunk_SP(t-1,i),Sand_Fraction,Clay_Fraction,...
            Fraction_Organic_Matter,Soil_Water(t),SoilWaterFlag,CanWetFrac(t-1,i),CanInt(t-1,i));
    end
    LW_in_bc_CLM_comb(t,j) = LW_in_bc_CLM(t,1,j) * PFT_Fractions(1) ...
        + LW_in_bc_CLM(t,2,j) * PFT_Fractions(2);
end
CanInt(gap_1(end)+1,:) = CanInt(gap_1(1)-1,:);
CanWetFrac(gap_1(end)+1,:) = CanWetFrac(gap_1(1)-1,:);
T_veg_CLM(gap_1(end)+1,:) = T_veg_CLM(gap_1(1)-1,:);
T_leaf_SP(gap_1(end)+1,:) = T_leaf_SP(gap_1(1)-1,:);
T_trunk_SP(gap_1(end)+1,:) = T_trunk_SP(gap_1(1)-1,:);

for t=gap_1(end)+2:2500
%{
Determine fraction of direct & diffuse insolation. In the absence of
specific measurements we use the effective emissivity of the sky, which is
a proxy for cloudiness. We scale emissivity to a range of 0 to 1 using its
own minimum and maximum (over all 4 years).
Note that the impact of this fraction on the sub-canopy LW radiation is
marginal.
%}
    forc_sol_dir(1) = SW_In_Above_Vegetation(t)*(1-Diffuse_Fraction(t));
    forc_sol_dir(2) = 0;    % no nir radiation
    forc_sol_dif(1) = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
    forc_sol_dif(2) = 0;    % no nir radiation
    
%{
Essery et al. (2016) used threshold temperature of 2C and scaling factors
to determine snowfall. This is applied here as well. Air temperature
measurements were conducted at a tower but height is uncertain for
precipitation measurements - but likely on ground level, which is why we
adjust threshold temperature using CLM4.5's lapse rate.
Another solution would be to apply an algorithm used for SNowMIP2.
SnowMIP2.
%}
    Threshold_Temperature = Air_Temperature_Above_Vegetation(t) + 0.0098*Height_of_Temp_Obs;
    if  Threshold_Temperature < 273.15+2
        forc_snow = Precipitation(t);
        forc_rain = 0;
    else
        forc_snow = 0;
        forc_rain = Precipitation(t);
    end
        
% snow age for surface albedo parameterization
    if Snow_Depth_Below_Vegetation(t) > 0 && forc_snow == 0
        Snow_Age(t) = Snow_Age(t-1) + dt;   % snow age in days for ageing
    else
        Snow_Age(t) = Snow_Age(t-1) + dt;   % snow age in days for ageing
    end
% surface albedo parameterization via snow ageing
    Surface_Albedo(t) = (1-Snow_Fraction(t)) * Soil_Albedo_Direct...
        + Snow_Fraction(t)*((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3);
    Surface_Albedo(t) = min(1,max(0,(0.8+(j-1)*0.04)*Surface_Albedo(t)));
    
    for i=1:2
        [Cosine_Solar_Angle(t),CanInt(t,i),CanWetFrac(t,i),...
            SW_net_veg_CLM(t,i),LW_in_bc_CLM(t,i,j),LW_out_ac_CLM(t,i),...
            T_veg_CLM(t,i),T_air_wc_CLM(t,i),T_ref2m_CLM(t,i),...
            EB_CLM(t,:,:,i),EB_SP_leaf(t,:,:,i),EB_SP_trunk(t,:,:,i),...
            SW_net_can_SP(t,i),SW_net_trunk_SP(t,i),LW_in_bc_SP(t,i,j),...
            LW_out_ac_SP(t,i),T_leaf_SP(t,i),T_trunk_SP(t,i)] = ...
        SNOWPACKwithinCLM45driver(Latitude,Longitude,Elevation,JulDay_1h(t),...
            Timestep_CLM,Timestep_SP,PFT(i),LAI(i),SAI(i),Vegetation_Height,...
            Basal_Area_of_Forest_Stand(i),Tree_Diameter(i),Fraction_LAI_Top,K_LAI,...
            Fraction_Trunk_Height,SW_to_Trunk,Height_of_Hum_Obs,...
            Height_of_Temp_Obs,Height_of_Wind_Obs,CLM45_canopy_layers,HM_in_CLM,...
            forc_sol_dir,forc_sol_dif,NIR_available,LW_In_Above_Vegetation(t),...
            Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
            Relative_Humidity_Above_Vegetation(t),forc_rain,forc_snow,...
            Snow_Depth_Below_Vegetation(t),Snow_Fraction(t),Fraction_Frozen_Soil(t),...
            Surface_Albedo(t),Surface_Temperature(t),T_veg_CLM(t-1,i),...
            T_leaf_SP(t-1,i),T_trunk_SP(t-1,i),Sand_Fraction,Clay_Fraction,...
            Fraction_Organic_Matter,Soil_Water(t),SoilWaterFlag,CanWetFrac(t-1,i),CanInt(t-1,i));
    end
    LW_in_bc_CLM_comb(t,j) = LW_in_bc_CLM(t,1,j) * PFT_Fractions(1) ...
        + LW_in_bc_CLM(t,2,j) * PFT_Fractions(2);
end
end

save('LWsub_Borden_AlbSurf.mat','LW_in_bc_1h','LW_in_bc_CLM_comb','LW_in_bc_SP')

toc