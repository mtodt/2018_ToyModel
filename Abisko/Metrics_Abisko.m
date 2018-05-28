%{
1) Run DataPrep_Abisko.m.
2) Run TodModelAtAbisko.m with HM_in_CLM set to 'yes'.
3) Change LW_in_bc_CLM to LW_in_bc_CLM_HM and T_veg_CLM to T_veg_CLM_HM,
   and save both in CLM_Abisko_HM.mat.
   LW_in_bc_CLM_HM = LW_in_bc_CLM;
   T_veg_CLM_HM = T_veg_CLM;
   save('CLM_Abisko_HM.mat','LW_in_bc_CLM_HM','T_veg_CLM_HM')
4) Run TodModelAtAbisko.m again with HM_in_CLM set to 'no'.
%}
load CLM_Abisko_HM.mat

%-------------------------------------------------------------------------%
%----------------  select evaluation periods based on QC  ----------------%
%{
Data have been quality-controlled using flags in the original excel
spreadsheet. To prevent any biases we only consider full days, so we
exclude any residuals after quality control.
%}
    
%---------------  9 days of sub-canopy longwave radiation  ---------------%
EP1 = 11:34;     % evaluation period 1 - 1 day
EP2 = 107:130;   % evaluation period 2 - 1 day
EP3 = 179:226;   % evaluation period 3 - 2 days
EP4 = 395:442;   % evaluation period 4 - 2 days
EP5 = 515:586;   % evaluation period 5 - 3 days

LW_sub_val_eval = vertcat(LW_in_bc_C_1h(EP1(1):EP1(end),:),...
    LW_in_bc_C_1h(EP2(1):EP2(end),:),...
    LW_in_bc_C_1h(EP3(1):EP3(end),:),...
    LW_in_bc_C_1h(EP4(1):EP4(end),:),...
    LW_in_bc_C_1h(EP5(1):EP5(end),:));
LW_sub_val_eval_all = vertcat(LW_sub_val_eval(:,1),LW_sub_val_eval(:,2),...
    LW_sub_val_eval(:,3),LW_sub_val_eval(:,4));
LW_sub_val_eval_avg = nan(length(LW_sub_val_eval(:,1)),1);
for l=1:length(LW_sub_val_eval_avg)
    LW_sub_val_eval_avg(l) = mean(LW_sub_val_eval(l,:));
end

LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP1(1):EP1(end),:),...
    LW_in_bc_CLM(EP2(1):EP2(end),:),...
    LW_in_bc_CLM(EP3(1):EP3(end),:),...
    LW_in_bc_CLM(EP4(1):EP4(end),:),...
    LW_in_bc_CLM(EP5(1):EP5(end),:));
LW_sub_CLM_eval_all = vertcat(LW_sub_CLM_eval(:,1),LW_sub_CLM_eval(:,2),...
    LW_sub_CLM_eval(:,3),LW_sub_CLM_eval(:,4));
LW_sub_CLM_eval_avg = nan(length(LW_sub_CLM_eval(:,1)),1);
for l=1:length(LW_sub_CLM_eval_avg)
    LW_sub_CLM_eval_avg(l) = mean(LW_sub_CLM_eval(l,:));
end

LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP1(1):EP1(end),:),...
    LW_in_bc_CLM_HM(EP2(1):EP2(end),:),...
    LW_in_bc_CLM_HM(EP3(1):EP3(end),:),...
    LW_in_bc_CLM_HM(EP4(1):EP4(end),:),...
    LW_in_bc_CLM_HM(EP5(1):EP5(end),:));
LW_sub_CLMHM_eval_all = vertcat(LW_sub_CLMHM_eval(:,1),LW_sub_CLMHM_eval(:,2),...
    LW_sub_CLMHM_eval(:,3),LW_sub_CLMHM_eval(:,4));
LW_sub_CLMHM_eval_avg = nan(length(LW_sub_CLMHM_eval(:,1)),1);
for l=1:length(LW_sub_CLMHM_eval_avg)
    LW_sub_CLMHM_eval_avg(l) = mean(LW_sub_CLMHM_eval(l,:));
end

LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP1(1):EP1(end),:),...
    LW_in_bc_SP(EP2(1):EP2(end),:),...
    LW_in_bc_SP(EP3(1):EP3(end),:),...
    LW_in_bc_SP(EP4(1):EP4(end),:),...
    LW_in_bc_SP(EP5(1):EP5(end),:));
LW_sub_SP_eval_all = vertcat(LW_sub_SP_eval(:,1),LW_sub_SP_eval(:,2),...
    LW_sub_SP_eval(:,3),LW_sub_SP_eval(:,4));
LW_sub_SP_eval_avg = nan(length(LW_sub_SP_eval(:,1)),1);
for l=1:length(LW_sub_SP_eval_avg)
    LW_sub_SP_eval_avg(l) = mean(LW_sub_SP_eval(l,:));
end

%-------------------------------------------------------------------------%
%------------------------  calculate RMSE and MBD  ------------------------

RMSE_Abisko_CLM = nan(1,6); MBD_Abisko_CLM = nan(1,6);
RMSE_Abisko_CLMHM = nan(1,6); MBD_Abisko_CLMHM = nan(1,6);
for i=1:4
    RMSE_Abisko_CLM(i) = RMSE(length(LW_sub_val_eval(:,i)),LW_sub_CLM_eval(:,i),...
        LW_sub_val_eval(:,i));
    MBD_Abisko_CLM(i) = MBD(length(LW_sub_val_eval(:,i)),LW_sub_CLM_eval(:,i),...
        LW_sub_val_eval(:,i));
    RMSE_Abisko_CLMHM(i) = RMSE(length(LW_sub_val_eval(:,i)),LW_sub_CLMHM_eval(:,i),...
        LW_sub_val_eval(:,i));
    MBD_Abisko_CLMHM(i) = MBD(length(LW_sub_val_eval(:,i)),LW_sub_CLMHM_eval(:,i),...
        LW_sub_val_eval(:,i));
end
RMSE_Abisko_CLM(5) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLM_eval_all,...
    LW_sub_val_eval_all);
MBD_Abisko_CLM(5) = MBD(length(LW_sub_val_eval_all),LW_sub_CLM_eval_all,...
    LW_sub_val_eval_all);
RMSE_Abisko_CLMHM(5) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLMHM_eval_all,...
    LW_sub_val_eval_all);
MBD_Abisko_CLMHM(5) = MBD(length(LW_sub_val_eval_all),LW_sub_CLMHM_eval_all,...
    LW_sub_val_eval_all);
RMSE_Abisko_CLM(6) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg,...
    LW_sub_val_eval_avg);
MBD_Abisko_CLM(6) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg,...
    LW_sub_val_eval_avg);
RMSE_Abisko_CLMHM(6) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLMHM_eval_avg,...
    LW_sub_val_eval_avg);
MBD_Abisko_CLMHM(6) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLMHM_eval_avg,...
    LW_sub_val_eval_avg);
