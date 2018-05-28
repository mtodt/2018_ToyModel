%{
1) Run DataPrep_Sodankyla.m.
2) Run TodModelAtSodankyla.m with HM_in_CLM set to 'yes'.
3) Change LW_in_bc_CLM to LW_in_bc_CLM_HM and T_veg_CLM to T_veg_CLM_HM,
   and save both in CLM_Sodankyla_HM.mat.
        LW_in_bc_CLM_HM = LW_in_bc_CLM;
        T_veg_CLM_HM = T_veg_CLM;
        save('CLM_Sodankyla_HM.mat','LW_in_bc_CLM_HM','T_veg_CLM_HM')
8) Run TodModelAtSodankyla.m again with HM_in_CLM set to 'no'.
%}
load CLM_Sodankyla_C_HM.mat

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
%{
Reid et al. (2014) give different values of available duration of
measurements for pre- and post-quality control. However, according to Nick
Rutter (co-author of previously mentioned paper), radiometers at Sodankyla
were checked daily so that none of the measurements should be questionable.
Location C:
Measurements for evaluation available from 9 March 2012 15:00 to 17 April
2012 14:00. However, time steps 457-459 need to be excluded due to
interpolation of LWR and SWR, which corresponds to 14:00 - 16:00 on 28
March 2012. To dismiss these and keep equal amounts of data for every hour,
evaluation period is set to 10 March 1:00 to 17 April 0:00 and all of 28
March is skipped.
Location R4:
Instruments had been moved to site R4 on 31 March 2012 14:00 and were
removed on 16 April 2012 5:00, so that evaluation period is set to 1 April
1:00 to 16 April 2012 0:00.
%}
% C: select 10 March 1:00 to 28 March 0:00
EP_C_1 = 12:443;
% C: add 29 March 1:00 to 17 April 0:00
EP_C_2 = 468:923;

LW_In_Above_Vegetation = LW_In_Above_Vegetation';
LW_atm_eval = vertcat(LW_In_Above_Vegetation(EP_C_1(1):EP_C_1(end)),...
    LW_In_Above_Vegetation(EP_C_2(1):EP_C_2(end)));
LW_atm_eval_x4 = vertcat(LW_atm_eval,LW_atm_eval,LW_atm_eval,LW_atm_eval);

LW_sub_val_eval = vertcat(LW_in_bc_1h_C(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_1h_C(EP_C_2(1):EP_C_2(end),:));
LW_sub_val_eval_all = vertcat(LW_sub_val_eval(:,1),LW_sub_val_eval(:,2),...
    LW_sub_val_eval(:,3),LW_sub_val_eval(:,4));
LW_sub_val_eval_avg = nan(length(LW_sub_val_eval(:,1)),1);
for l=1:length(LW_sub_val_eval_avg)
    LW_sub_val_eval_avg(l) = mean(LW_sub_val_eval(l,:));
end

LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_CLM(EP_C_2(1):EP_C_2(end),:));
LW_sub_CLM_eval_all = vertcat(LW_sub_CLM_eval(:,1),LW_sub_CLM_eval(:,2),...
    LW_sub_CLM_eval(:,3),LW_sub_CLM_eval(:,4));
LW_sub_CLM_eval_avg = nan(length(LW_sub_CLM_eval(:,1)),1);
for l=1:length(LW_sub_CLM_eval_avg)
    LW_sub_CLM_eval_avg(l) = mean(LW_sub_CLM_eval(l,:));
end

LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_CLM_HM(EP_C_2(1):EP_C_2(end),:));
LW_sub_CLMHM_eval_all = vertcat(LW_sub_CLMHM_eval(:,1),LW_sub_CLMHM_eval(:,2),...
    LW_sub_CLMHM_eval(:,3),LW_sub_CLMHM_eval(:,4));
LW_sub_CLMHM_eval_avg = nan(length(LW_sub_CLMHM_eval(:,1)),1);
for l=1:length(LW_sub_CLMHM_eval_avg)
    LW_sub_CLMHM_eval_avg(l) = mean(LW_sub_CLMHM_eval(l,:));
end

LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_SP(EP_C_2(1):EP_C_2(end),:));
LW_sub_SP_eval_all = vertcat(LW_sub_SP_eval(:,1),LW_sub_SP_eval(:,2),...
    LW_sub_SP_eval(:,3),LW_sub_SP_eval(:,4));
LW_sub_SP_eval_avg = nan(length(LW_sub_SP_eval(:,1)),1);
for l=1:length(LW_sub_SP_eval_avg)
    LW_sub_SP_eval_avg(l) = mean(LW_sub_SP_eval(l,:));
end

%-------------------------------------------------------------------------%
%------------------------  calculate RMSE and MBD  ------------------------
RMSE_Sodankyla_SP = nan(1,6); MBD_Sodankyla_SP = nan(1,6);
RMSE_Sodankyla_CLM = nan(1,6); MBD_Sodankyla_CLM = nan(1,6);
RMSE_Sodankyla_CLMHM = nan(1,6); MBD_Sodankyla_CLMHM = nan(1,6);
for i=1:4
    RMSE_Sodankyla_SP(i) = RMSE(length(LW_sub_val_eval(:,i)),LW_sub_SP_eval(:,i),...
        LW_sub_val_eval(:,i));
    MBD_Sodankyla_SP(i) = MBD(length(LW_sub_val_eval(:,i)),LW_sub_SP_eval(:,i),...
        LW_sub_val_eval(:,i));
    RMSE_Sodankyla_CLM(i) = RMSE(length(LW_sub_val_eval(:,i)),LW_sub_CLM_eval(:,i),...
        LW_sub_val_eval(:,i));
    MBD_Sodankyla_CLM(i) = MBD(length(LW_sub_val_eval(:,i)),LW_sub_CLM_eval(:,i),...
        LW_sub_val_eval(:,i));
    RMSE_Sodankyla_CLMHM(i) = RMSE(length(LW_sub_val_eval(:,i)),LW_sub_CLMHM_eval(:,i),...
        LW_sub_val_eval(:,i));
    MBD_Sodankyla_CLMHM(i) = MBD(length(LW_sub_val_eval(:,i)),LW_sub_CLMHM_eval(:,i),...
        LW_sub_val_eval(:,i));
end
RMSE_Sodankyla_SP(5) = RMSE(length(LW_sub_val_eval_all),LW_sub_SP_eval_all,...
    LW_sub_val_eval_all);
MBD_Sodankyla_SP(5) = MBD(length(LW_sub_val_eval_all),LW_sub_SP_eval_all,...
    LW_sub_val_eval_all);
RMSE_Sodankyla_CLM(5) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLM_eval_all,...
    LW_sub_val_eval_all);
MBD_Sodankyla_CLM(5) = MBD(length(LW_sub_val_eval_all),LW_sub_CLM_eval_all,...
    LW_sub_val_eval_all);
RMSE_Sodankyla_CLMHM(5) = RMSE(length(LW_sub_val_eval_all),LW_sub_CLMHM_eval_all,...
    LW_sub_val_eval_all);
MBD_Sodankyla_CLMHM(5) = MBD(length(LW_sub_val_eval_all),LW_sub_CLMHM_eval_all,...
    LW_sub_val_eval_all);
RMSE_Sodankyla_SP(6) = RMSE(length(LW_atm_eval),LW_sub_SP_eval_avg,...
    LW_sub_val_eval_avg);
MBD_Sodankyla_SP(6) = MBD(length(LW_sub_val_eval_avg),LW_sub_SP_eval_avg,...
    LW_sub_val_eval_avg);
RMSE_Sodankyla_CLM(6) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg,...
    LW_sub_val_eval_avg);
MBD_Sodankyla_CLM(6) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLM_eval_avg,...
    LW_sub_val_eval_avg);
RMSE_Sodankyla_CLMHM(6) = RMSE(length(LW_sub_val_eval_avg),LW_sub_CLMHM_eval_avg,...
    LW_sub_val_eval_avg);
MBD_Sodankyla_CLMHM(6) = MBD(length(LW_sub_val_eval_avg),LW_sub_CLMHM_eval_avg,...
    LW_sub_val_eval_avg);
