function [qg,soilbeta,btran,btran2] = ...
    Biogeophysics1(PFT,SoilWater,SoilWaterFlag,om_frac,sand,clay,frac_sno,frac_soil_frozen,t_grnd,forc_pbot,forc_q)
%{
Variables
- watsat (volumetric soil water at saturation),
- sucsat (minimum soil suction [mm]),
- bsw (Clapp and Hornberger "b"),
- smpmin (restriction for min of soil potential [mm]) and
- hksat (hydraulic conductivity at saturation [mm H2O s^{-1}])
are calculated from soil composition data.
Calculations partially taken from iniTimeConst.f90 and CanopyFluxesMod.f90.
    
Variables om_watsat, om_sucsat, om_b, and om_hksat are depending on soil
layer. However, the values used here are vertical averages considering
changing layer thickness since the Toy Model doesn't use multiple soil
layers.
%}
grav = 9.80616;             % acceleration of gravity [m s^{-2}]

%------------------------  from iniTimeConst.f90  ------------------------%
om_watsat = 0.8424;         % initially 0.88 -> [0.83,0.93]
om_sucsat = 9.8978;         % initially 10.2 -> [10.1,10.3] - misread min/max
om_b = 10.8453;             % initially 5.35 -> [2.7,12]
om_hksat = 0.0349;          % initially 0.14 -> [0.0001,0.28]
watsat = (1-om_frac)*(0.489-0.00126*sand*100) + om_watsat*om_frac;
sucsat = (1-om_frac)*(10*10^(1.88-0.0131*sand*100)) + om_sucsat*om_frac;
bsw = (1-om_frac)*(2.91+0.159*clay*100) + om_frac*om_b;
xksat = 0.0070556*10^(-0.884+0.0153*sand*100);
if om_frac>0.5
    perc_frac = (1-0.5)^(-0.139) * (om_frac-0.5)^0.139;
else
    perc_frac = 0;
end
if om_frac<1
    uncon_hksat = ((1-om_frac)+(1-perc_frac)*om_frac)...
        /((1-om_frac)/xksat + ((1-perc_frac)*om_frac)/om_hksat);
else
    uncon_hksat = 0;
end
hksat = ((1-om_frac)+(1-perc_frac)*om_frac)*uncon_hksat + (perc_frac*om_frac)*om_hksat;
watfc = watsat*(0.1/(hksat*86400))^(1/(2*bsw+3));

%-----------------------  from Biogeophysics1.f90  -----------------------%
%{
Soil moisture or water content measurements might not always be available
and we have to use a proxy, e.g. based on ground water level. As these
proxies are just dimensionless fractions and not partial volumes saturation
water content is scaled by these fractions and so we add a flag to
differentiate between proxies and direct measurements.
%}
if strcmp(SoilWaterFlag,'proxy') == 1
    wx = SoilWater*watsat; % soil water content or moisture just unitless fraction, not partial volume
elseif strcmp(SoilWaterFlag,'observation') == 1
    wx = SoilWater;
end
fac = max(0.01,min(1,wx/watsat)); % wx relative to watsat <-'
smpmin = -1*10^8;
psit = max(smpmin,-sucsat*fac^(-bsw));
roverg = (6.02214*10^26*1.38065*10^(-23)/18.016)/grav * 1000; % [mm K^{-1}]
hr = exp(psit/roverg/t_grnd);
qred = (1-frac_sno)*hr + frac_sno;  % frac_h2osfc disregarded
if wx < watfc
    fac_fc = max(0.01,min(1,wx/watfc));
    soilbeta = (1-frac_sno)*0.25*(1-cos(pi*fac_fc))^2 + frac_sno;  % frac_h2osfc disregarded
else
    soilbeta = 1;
end
soilalpha = qred;
[eg,degdT,qsatg,qsatgdT] = QSat(t_grnd,forc_pbot);
qg = qred*qsatg;
dqgdT = qred*qsatgdT;
if qsatg > forc_q && forc_q > qred*qsatg   % only works for frac_sno either 1 or 0...
    qg = forc_q;                           % ...otherwise split into snow and soil cases
    dqgdT = 0;
end

%------------------------  from CanopyFluxes.f90  ------------------------%
% transpiration wetness factor (0 to 1)
btran  = 0;
btran2  = 0;
%{
"Effective porosity of soil, partial volume of ice and liquid (needed for
btran) and root resistance factors"
...are calculated here to get btran, which is necessary for the calculation
of transpiration, etc. later on. A value of btran=0 would limit the
calculations and thus be disadvantageous for the assessment of CLM4.5.
Therefore, btran is calculated out of proxy values describing the soil
water content.
%}
% soil water potential at full stomatal closure [mm]
smpsc = [0 -255000 -255000 -255000 -255000 -255000 -224000 -224000 ...
    -224000 -428000 -428000 -428000 -275000 -275000 -275000 -275000 ...
    -275000 -275000 -275000 -275000 -275000];
% soil water potential at full stomatal opening [mm]
smpso = [0 -66000 -66000 -66000 -66000 -66000 -35000 -35000 -35000 ...
    -83000 -83000 -83000 -74000 -74000 -74000 -74000 -74000 -74000 ...
    -74000 -74000 -74000];
tfrz = 273.15;
% conversion of water content proxy to fractions of solid and liquid water
if strcmp(SoilWaterFlag,'proxy') == 1
    h2osoi_vol = SoilWater*watsat;
elseif strcmp(SoilWaterFlag,'observation') == 1
    h2osoi_vol = SoilWater;
end
%{
We don't consider soil layers here, as measurements for different locations
would feature different measurement depths and prevent generalization.
However, we use a prescribed fraction of frozen soil, which is calculated
out of soil temperatures if available or assumed otherwise, i.e. 0 or 1 as
alpine sites pretty much unfrozen.
Note that CLM's unit for h2osoi_vol is [kg m^{-2}] and it gets converted
into a unitless value via layer depth and denisty of ice and water. As
measurements are given in [m^3 m^{-3}] we use this value directly.
%}
vol_ice_proxy = h2osoi_vol * frac_soil_frozen;
vol_liq_proxy = h2osoi_vol * (1-frac_soil_frozen);
% Root resistance factors
vol_ice = min(watsat,vol_ice_proxy);
eff_porosity = watsat-vol_ice;
vol_liq = min(eff_porosity,vol_liq_proxy);
rootfr = 0.9949; % since no differentiation of soil layers -> all of soil considered, but due to design total rootfraction not 1
%{
roota_par = 7;      % CLM rooting distribution parameter [1/m], value for NBTs, 6 for DBTs
rootb_par = 2;      % CLM rooting distribution parameter [1/m], value for all boreal trees
%}
if vol_liq <= 0 || frac_soil_frozen == 1
    rootr = 0;
else
    s_node = max(vol_liq/eff_porosity,0.01);
    smp_node = max(smpsc(PFT),-sucsat*s_node^(-bsw));
    rresis = min((eff_porosity/watsat)*(smp_node - smpsc(PFT))/(smpso(PFT) - smpsc(PFT)),1);
    rootr = rootfr*rresis;
    btran = btran + rootr;
    smp_node_lf = max(smpsc(PFT),-sucsat*(h2osoi_vol/watsat)^(-bsw)) ;
    btran2 = btran2 + rootfr*min((smp_node_lf - smpsc(PFT))/(smpso(PFT) - smpsc(PFT)),1);
end
%{
% Normalize root resistances to get layer contribution to ET
if btran > 0
    rootr = rootr/btran;
else
    rootr = 0;
end
%}
end