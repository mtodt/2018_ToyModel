%{
1) Run DataPrep_Yakutsk.m.
2) Run ToyModelAtYakutsk.m with HM_in_CLM set to 'yes'.
3) Change LW_in_bc_CLM to LW_in_bc_CLM_HM and T_veg_CLM to T_veg_CLM_HM,
   and save both in CLM_Yakutsk_HM.mat.
   LW_in_bc_CLM_HM_98 = LW_in_bc_CLM_98;
   T_veg_CLM_HM_98 = T_veg_CLM_98;
   save('CLM_Yakutsk_HM.mat','LW_in_bc_CLM_HM_98','T_veg_CLM_HM_98')
4) Run ToyModelAtYakutsk.m again with HM_in_CLM set to 'no'.
%}
load CLM_Yakutsk_HM.mat

%-------------------------------------------------------------------------%
%----------------------  calculate validation data  ----------------------%
% limit unrealistic sub-canopy SWR
SW_out_bc_98_1h_filled_corrected = nan(size(SW_out_bc_98_1h_filled));
for t=1:length(Rnet_bc_98_1h_filled)
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
em_gr_98 = nan(size(Rnet_bc_98_1h_filled));
em_gr_98(1:3195) = em_snow; em_gr_98(3196:end) = em_soil;
LW_in_bc_98 = nan(size(Rnet_bc_98_1h_filled));
for t=1:length(Rnet_bc_98_1h_filled)
    LW_in_bc_98(t) = LW_net_bc_98(t) + em_gr_98(t)*boltz*T_surf_98_1h_filled(t)^4;
end

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
EP_98_1 = 1057:1728;
EP_98_2 = 1777:3192;

LW_sub_val_eval = vertcat(LW_in_bc_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_98(EP_98_2(1):EP_98_2(end)));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_CLM_98(EP_98_2(1):EP_98_2(end)));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_CLM_HM_98(EP_98_2(1):EP_98_2(end)));
LW_sub_SP_eval = vertcat(LW_in_bc_SP_98(EP_98_1(1):EP_98_1(end)),...
    LW_in_bc_SP_98(EP_98_2(1):EP_98_2(end)));

SWR_Incoming_Yakutsk = vertcat(SW_in_ac_98_1h_filled(EP_98_1(1):EP_98_1(end),:),...
    SW_in_ac_98_1h_filled(EP_98_2(1):EP_98_2(end),:));

i=0;
for l=1:length(SWR_Incoming_Yakutsk)
    if SWR_Incoming_Yakutsk(l) == 0
        i=i+1;
        LW_sub_val_eval_Night(i,:) = LW_sub_val_eval(l,:);
        LW_sub_CLM_eval_Night(i,:) = LW_sub_CLM_eval(l,:);
        LW_sub_CLMHM_eval_Night(i,:) = LW_sub_CLMHM_eval(l,:);
        LW_sub_SP_eval_Night(i,:) = LW_sub_SP_eval(l,:);
    end
end

%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%

RMSE_Yakutsk_CLM = RMSE(length(LW_sub_val_eval_Night),LW_sub_CLM_eval_Night,...
    LW_sub_val_eval_Night);
RMSE_Yakutsk_CLMHM = RMSE(length(LW_sub_val_eval_Night),LW_sub_CLMHM_eval_Night,...
    LW_sub_val_eval_Night);
MBD_Yakutsk_CLM = MBD(length(LW_sub_val_eval_Night),LW_sub_CLM_eval_Night,...
    LW_sub_val_eval_Night);
MBD_Yakutsk_CLMHM = MBD(length(LW_sub_val_eval_Night),LW_sub_CLMHM_eval_Night,...
    LW_sub_val_eval_Night);