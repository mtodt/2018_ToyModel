%{
1) Run DataPrep_Cherskiy.m.
2) Run ToyModelAtCherskiy.m with HM_in_CLM set to 'yes'.
3) Change LW_in_bc_CLM to LW_in_bc_CLM_HM and T_veg_CLM to T_veg_CLM_HM,
   and save both in CLM_Cherskiy_HM.mat.
   LW_in_bc_CLM_HM = LW_in_bc_CLM;
   T_veg_CLM_HM = T_veg_CLM;
   save('CLM_Cherskiy_HM.mat','LW_in_bc_CLM_HM','T_veg_CLM_HM')
4) Run ToyModelAtCherskiy.m again with HM_in_CLM set to 'no'.
%}
load CLM_Cherskiy_HM.mat

%-------------------------------------------------------------------------%
%----------------  select evaluation periods based on QC  ----------------%
%{
Measurements from upward looking radiometers have been checked for snow 
covering the sensors. This has apparently been the case for longwave 
radiation above the canopy during several hours (9:00 to 16:00) on 17 May 
2017, however, we used those measurements to continue driving the Toy 
Model. All of 17 May 2017 will be excluded from evaluation to guarantee 
false measurements don't affect the analysis.
Also, we only consider full days to prevent any biases so we additionally 
exclude the first couple of hours (29 March 14:00 to 30 March 2017 0:00).
%}
EP1 = 12:1163;
EP2 = 1188:length(time_1h);

LW_sub_val_eval = horzcat(LW_in_bc_1h(EP1(1):EP1(end)),...
    LW_in_bc_1h(EP2(1):EP2(end)));
LW_sub_CLM_eval = horzcat(LW_in_bc_CLM(EP1(1):EP1(end)),...
    LW_in_bc_CLM(EP2(1):EP2(end)));
LW_sub_CLMHM_eval = horzcat(LW_in_bc_CLM_HM(EP1(1):EP1(end)),...
    LW_in_bc_CLM_HM(EP2(1):EP2(end)));
LW_sub_SP_eval = horzcat(LW_in_bc_SP(EP1(1):EP1(end)),...
    LW_in_bc_SP(EP2(1):EP2(end)));

% calculate RMSE and MBD
RMSE_Cherskiy_CLM = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,...
    LW_sub_val_eval);
RMSE_Cherskiy_CLMHM = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,...
    LW_sub_val_eval);
MBD_Cherskiy_CLM = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,...
    LW_sub_val_eval);
MBD_Cherskiy_CLMHM = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,...
    LW_sub_val_eval);
