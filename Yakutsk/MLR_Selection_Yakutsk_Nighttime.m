%{
Collect ingredients necessary for multi-linear regression.
Emissivity of the sky & solar angle -> description of meteorological
conditions (insolation doesn't translate from simulation to simulation
because of differences in latitude and season, and it's also depending on
cloudiness).
%}

%----------------------  calculate validation data  ----------------------%
% limit unrealistic sub-canopy SWR
SW_out_bc_98_1h_filled_corrected = nan(size(SW_out_bc_98_1h_filled));
for t=1:length(time_98_1h)
    SW_out_bc_98_1h_filled_corrected(t) ...
        = min(SW_out_bc_98_1h_filled(t),SW_in_bc_98_1h_filled(t));
end
SW_net_bc_98_1h_filled_corrected = nan(size(Rnet_bc_98_1h_filled));
for l=1:length(SW_net_bc_98_1h_filled_corrected)
    SW_net_bc_98_1h_filled_corrected(l) ...
        = SW_in_bc_98_1h_filled(l) - SW_out_bc_98_1h_filled_corrected(l);
end

% calculate net sub-canopy LWR
LW_net_bc_98 = Rnet_bc_98_1h_filled - SW_net_bc_98_1h_filled_corrected;

% subtract outgoing LWR calculated from surface temperature
boltz = 5.67*10^(-8);
em_soil = 0.96; em_snow = 0.97;
em_gr_98 = nan(size(time_98_1h));
em_gr_98(1:3195) = em_snow; em_gr_98(3196:end) = em_soil;
LW_in_bc_98 = nan(size(time_98_1h));
for t=1:length(time_98_1h)
    LW_in_bc_98(t) = LW_net_bc_98(t) + em_gr_98(t)*boltz*T_surf_98_1h_filled(t)^4;
end

%----------------------  select evaluation periods  ----------------------%
%{
Since we prefer full days, so that diurnal cycles can be calculated and
diurnal biases prevented generally, and meltout can be determined
from top soil temperature for the morning of 14 May 1998, we exclude that
last day and set the end of evaluation to 14 May 1998 0:00, respectively.
It is more difficult for 2000, as there are clear diurnal cycles of top
soil temperature for 3-5 May 2000, however, their ranges are nowhere close
to that of 6 May and afterwards. Since both insolation and air temperature
are similar for all days from 3 May to 6 May we still include that period
and set the end of evaluation to 6 May 2000 0:00.
Evaluation period in 2000 begins on 19 April 1:00.
For 1998, measurements start in the afternoon of 13 February, so evaluation
starts on 14 February 1:00. However, there's a gap from 14 March 15:00 to
15 March 12:00, so that we limit evaluation to 14 February 1:00 to 14 March
0:00 and 16 March 1:00 to 14 May 1998 0:00.

Caveat:
Sub-canopy shortwave radiation is not reliable and only the earlier period
in 1998 seems somewhat sensible, so we only use this period for the
multi-linear regression.
%}
EP_98_1 = 1057:1728;
EP_98_2 = 1777:3192;

% select variables
LAI_Yakutsk = LAI;  % 0 throughout anyway
PAI_Yakutsk = LAI_Yakutsk + SAI;
Hourly_Time_Yakutsk = vertcat(time_98_1h(EP_98_1(1):EP_98_1(end)),...
    time_98_1h(EP_98_2(1):EP_98_2(end)));
Cosine_Solar_Angle_Yakutsk = vertcat(Cosine_Solar_Angle_98(EP_98_1(1):EP_98_1(end)),...
    Cosine_Solar_Angle_98(EP_98_2(1):EP_98_2(end)));
SWR_Incoming_Yakutsk = vertcat(SW_in_ac_98_1h_filled(EP_98_1(1):EP_98_1(end)),...
    SW_in_ac_98_1h_filled(EP_98_2(1):EP_98_2(end)));
LWR_Atmosphere_Yakutsk = vertcat(LW_in_ac_98_1h_filled(EP_98_1(1):EP_98_1(end)),...
    LW_in_ac_98_1h_filled(EP_98_2(1):EP_98_2(end)));

bla = 5.67*10^(-8)*T_air_98_1h_filled.^4;
blub = LW_in_ac_98_1h_filled./bla;
Sky_Emissivity_Yakutsk = vertcat(blub(EP_98_1(1):EP_98_1(end)),...
    blub(EP_98_2(1):EP_98_2(end)));
clear bla blub

LWR_Subcanopy_Observations_Yakutsk = vertcat(LW_in_bc_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_98(EP_98_2(1):EP_98_2(end)));
LWR_Subcanopy_CLM_Yakutsk = vertcat(LW_in_bc_CLM_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_CLM_98(EP_98_2(1):EP_98_2(end)));

i=0;
for l=1:length(SWR_Incoming_Yakutsk)
    if SWR_Incoming_Yakutsk(l) == 0
        i=i+1;
        Hourly_Time_Yakutsk_Night(i) = Hourly_Time_Yakutsk(l);
        Cosine_Solar_Angle_Yakutsk_Night(i) =  Cosine_Solar_Angle_Yakutsk(l);
        SWR_Incoming_Yakutsk_Night(i) =  SWR_Incoming_Yakutsk(l);
        LWR_Atmosphere_Yakutsk_Night(i) =  LWR_Atmosphere_Yakutsk(l);
        Sky_Emissivity_Yakutsk_Night(i) =  Sky_Emissivity_Yakutsk(l);
        LWR_Subcanopy_Observations_Yakutsk_Night(i) =  LWR_Subcanopy_Observations_Yakutsk(l);
        LWR_Subcanopy_CLM_Yakutsk_Night(i) =  LWR_Subcanopy_CLM_Yakutsk(l);
    end
end

save('MLRdata_Yakutsk.mat','LAI_Yakutsk','PAI_Yakutsk','Hourly_Time_Yakutsk_Night',...
    'Cosine_Solar_Angle_Yakutsk_Night','SWR_Incoming_Yakutsk_Night',...
    'LWR_Atmosphere_Yakutsk_Night','Sky_Emissivity_Yakutsk_Night',...
    'LWR_Subcanopy_Observations_Yakutsk_Night','LWR_Subcanopy_CLM_Yakutsk_Night')