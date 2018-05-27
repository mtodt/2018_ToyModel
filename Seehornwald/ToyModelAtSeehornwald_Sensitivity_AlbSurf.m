tic
clear
close all
%{
The Toy Model is designed in a way to facilitate general usage, requiring
input of parameters for vegetation and ground as well as meteorological
forcing. This file gives an overview of how to drive the Toy Model.
%}
% run DataPrep_Seehornwald.m first
load ForcingData_ToyModel_Seehornwald.mat

%-------------------------------------------------------------------------%
%--------------------  Location & General Parameters  --------------------%
%-------------------------------------------------------------------------%
Longitude = 9+51/60+21/60/60;
Latitude = 46+48/60+55/60/60;
Elevation = 1640;
Height_of_Hum_Obs = 35;
Height_of_Temp_Obs = 35;
Height_of_Wind_Obs = 35;
dt = time_1h(2)-time_1h(1);
Timestep_CLM = dt*24*3600;  % for CLM timestep required in [s]
Timestep_SP = dt*24;        % for Alptal Precipitation in mm/hour and timestep 1 hour -> dt = 1
Timestep_CLM = 3600;        % timestep not exact in time_1h
Timestep_SP = 1;            % timestep not exact in time_1h
JulDay = nan(size(time_1h));
for t=1:length(time_1h)
    date = datestr(time_1h(t),'yyyy-mm-dd-HH-MM-SS');
    Jultemp = datenum(str2double(date(1:4)),str2double(date(6:7)),str2double(date(9:10)),...
        str2double(date(12:13)),str2double(date(15:16)),str2double(date(18:19)));
    JulDay(t) = Jultemp-datenum(str2double(date(1:4))-1,12,31,0,0,0);
end


%-------------------------------------------------------------------------%
%------------------------  Vegetation Parameters  ------------------------%
%-------------------------------------------------------------------------%
%{
It would seem sensible to use temperate rather than boreal needleleaf
evergreen PFT, however, CLM4.5's surface dataset for 1/4° lists the following
fractions (PFTs with 0% not mentioned):
16%     not vegetated
3%      needleleaf evergreen temperate tree (NETTs)
21.5%   needleleaf evergreen boreal tree (NEBTs)
2%      broadleaf deciduous temperate tree (BDTTs)
4.5%    broadleaf deciduous boreal tree (BDBTs)
38%     broadleaf deciduous boreal shrub (BDBSs)
10.5%   c3 arctic grass (C3AGs)
0.5%    c3 non-arctic grass
3%      c3 crop
Boreal PFTs are generally preferred in comparison to their temperate
equivalents.
%}
PFT = 3;                % spruce -> needleleaf evergreen
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
Vegetation_Height = 25; % from Zweifel et al. (2016), Webster: max 27m
LAI = 3.9;
Fraction_LAI_Of_VAI = 0.7705; % value from CLM4.5
PAI = LAI/Fraction_LAI_Of_VAI;
SAI = PAI*(1-Fraction_LAI_Of_VAI);

%-----------------------  only for SNOWPACK-2LHM  ------------------------%
    % from Zweifel et al. (2016)
Tree_Diameter = 0.4;    % diameter at breast height [m]
    % from www.wsl.ch/en/forest/forest-development-and-monitoring/
    %      long-term-forest-ecosystem-research-lwf/sites.html
Large_Trees = 498;          % dbh > 0.12m in 2006
Large_Trees_Diameter = 0.47;% "quadratic average diameter of the 100 thickest trees per ha"
Plot_Area = 0.6*100*100;	% 0.6 ha plot size
Tree_Basal_Area = pi*(Large_Trees_Diameter/2)^2;
Basal_Area_of_Forest_Stand = Tree_Basal_Area*Large_Trees/Plot_Area;

% calibrated for Alptal and kept here
Fraction_LAI_Top = 0.5;         % fraction of LAI attributed to uppermost layer
K_LAI = 0.75;                   % radiation transmissivity parameter [0.4-0.8]
Fraction_Trunk_Height = 0.2;    % fraction of tree height occupied by trunks


%-------------------------------------------------------------------------%
%--------------------------  Ground Parameters  --------------------------%
%-------------------------------------------------------------------------%
%{
% Determine albedo of ground from observations. Use midday values only for
% claridity.
fig=figure(1);
set(gcf,'Position',get(0,'ScreenSize'))
plot(JulDay,alb_gr_1h,'.k')
datetick('x','mmm')
xlabel('Time','FontSize',17,'FontWeight','bold')
ylabel('Albedo','FontSize',17,'FontWeight','bold')
title('Sub-Canopy Ground Albedo at Seehornwald (limited to 1)','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 7];
fig.PaperSize = [8 7];
print(fig,'-dpdf','-r600','AlbedoGround_Seehornwald.pdf')
spectrum_alb = 0:0.025:1;
spectrum_alb_xaxis = 0.0125:0.025:0.9875;
hist_alb = histogram(alb_gr_1h(12:24:end-12),spectrum_alb,'Normalization','probability');
hist_ground = hist_alb.Values;
fig=figure(2);
set(gcf,'Position',get(0,'ScreenSize'))
plot(spectrum_alb_xaxis,hist_ground,'k','LineWidth',2)
xlim([0 1])
ylim([0 0.08])
xlabel('Albedo','FontSize',17,'FontWeight','bold')
ylabel('Probablity','FontSize',17,'FontWeight','bold')
title('Mid-day Ground Albedo at Seehornwald (limited to 1)','FontSize',17,'FontWeight','bold')
set(gca,'FontSize',17,'FontWeight','bold','LineWidth',2)
set(gca,'XTick',0:0.1:1)
box on
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 7];
fig.PaperSize = [8 7];
print(fig,'-dpdf','-r600','AlbedoGroundPDF_Seehornwald.pdf')
%}
Soil_Albedo_Direct = 0.1875;
Soil_Albedo_Diffuse = Soil_Albedo_Direct;
%{
Since we're using a range of sites, for most of which there probably are no
soil composition information, it's easier to use data from CLM4.5 surface
dataset. Seehornwald corresponds to grid cell (33,584) in 0.31x0.23 surface
file. Averaged over whole soil column using layer thickness.
%}
Clay_Fraction = 0.308479;
Sand_Fraction = 0.522282;
Fraction_Organic_Matter = 0.0760;
%{
Soil moisture measurements are "soil water content" in %, so probably
partial volume of soil. Which should then be used as variable wx in CLM4.5,
which is "partial volume of ice and water of surface layer".
%}
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
SW_In_Above_Vegetation = SW_in_ac_1h_filled;
Diffuse_Fraction = 1-r_SW_FluxNet_1h;   % using available insolation
LW_In_Above_Vegetation = LW_in_ac_1h_filled;
Air_Temperature_Above_Vegetation = T_air_ac_1h_filled+273.15;
Wind_Speed_Above_Vegetation = Wind_ac_1h_filled;
Relative_Humidity_Above_Vegetation = RH_ac_1h_filled;
% Conversion from mm/h to mm/s as required for CLM4.5
Precipitation = Prec_ac_1h_filled/3600;

%-------------------------------  surface  -------------------------------%
Soil_Water = SWC_1h/100;    % convert % into fraction
Fraction_Frozen_Soil = 0;   % guess based on some soil(?) temperature
Snow_Depth_Below_Vegetation = z_snow_forest_1h;
% snow fraction determined by eye from surface albedo measurements
Snow_Fraction = frac_snow_1h;
% calculate surface temperature from outgoing LWR
Surface_Temperature = T_surf_1h_filled;


%-------------------------------------------------------------------------%
%---------------------------  Initializations  ---------------------------%
%-------------------------------------------------------------------------%
Snowfall = zeros(size(Precipitation));
Rainfall = zeros(size(Precipitation));
Snow_Age = zeros(size(Snow_Depth_Below_Vegetation));
Surface_Albedo = nan(size(Snow_Age));

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
EvalSeehornwald = nan(size(time_1h));
for j=1:11
for t=2:length(time_1h)
    if ~isnan(SW_In_Above_Vegetation(t)) && ~isnan(Surface_Temperature(t))
%{
Determine fraction of direct & diffuse insolation. In the absence of
specific measurements there are two options: using the effective emissivity
of the sky, as done for Alptal, or using further information as is the case
for Seehornwald since the FluxNet tower also provides potential SW
radiation (how calculated?). The ratio of observed to potential SW
radiation is then used as the fraction of direct radiation. Note that the
impact of this fraction on the sub-canopy LW radiation is marginal (as shown
for Alptal).
%}
        forc_sol_dir(1) = SW_In_Above_Vegetation(t)*(1-Diffuse_Fraction(t));
        forc_sol_dir(2) = 0;    % no nir radiation
        forc_sol_dif(1) = SW_In_Above_Vegetation(t)*Diffuse_Fraction(t);
        forc_sol_dif(2) = 0;    % no nir radiation
        
%{
Differentiating between snow and rain. Since Seehornwald should be similar
to Alptal, the same calculation is used.
%}
        frac_prec_snow = RainSnowPartitioningAlptal(Air_Temperature_Above_Vegetation(t));
        forc_snow = Precipitation(t)*frac_prec_snow;
        forc_rain = Precipitation(t)*(1-frac_prec_snow);
        
% snow age for surface albedo parameterization
        if Snow_Depth_Below_Vegetation(t) > 0 && forc_snow == 0
            Snow_Age(t) = Snow_Age(t-1) + 1/24;   % snow age in days for ageing
        elseif Snow_Depth_Below_Vegetation(t) < Snow_Depth_Below_Vegetation(t-1)
            Snow_Age(t) = Snow_Age(t-1) + 1/24;   % snow age in days for ageing
        end
        
% surface albedo parameterization via snow ageing
        Surface_Albedo(t) = (1-Snow_Fraction(t)) * Soil_Albedo_Direct...
            + Snow_Fraction(t)*((Fresh_Snow_Albedo - 0.3)*exp(-Snow_Age(t)/7) + 0.3);
    Surface_Albedo(t) = min(1,max(0,(0.8+(j-1)*0.04)*Surface_Albedo(t)));
        
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
        forc_sol_dir,forc_sol_dif,NIR_available,LW_In_Above_Vegetation(t),...
        Air_Temperature_Above_Vegetation(t),Wind_Speed_Above_Vegetation(t),...
        Relative_Humidity_Above_Vegetation(t),forc_rain,forc_snow,...
        Snow_Depth_Below_Vegetation(t),Snow_Fraction(t),Fraction_Frozen_Soil,...
        Surface_Albedo(t),Surface_Temperature(t),T_veg_CLM(t-1),...
        T_leaf_SP(t-1),T_trunk_SP(t-1),Sand_Fraction,Clay_Fraction,...
        Fraction_Organic_Matter,Soil_Water(t),SoilWaterFlag,CanWetFrac(t-1),CanInt(t-1));
    
        EvalSeehornwald(t) = 1;
    else
        CanInt(t) = CanInt(t-1);
        CanWetFrac(t) = CanWetFrac(t-1);
        T_veg_CLM(t) = T_veg_CLM(t-1);
        T_leaf_SP(t) = T_leaf_SP(t-1);
        T_trunk_SP(t) = T_trunk_SP(t-1);
    
    end
end
end

save('LWsub_Seehornwald_AlbSurf.mat','LW_in_bc_1h','LW_in_bc_CLM','LW_in_bc_SP')
toc