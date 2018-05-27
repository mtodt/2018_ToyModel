% evaluation period
EP_C_1 = 12:443;
EP_C_2 = 468:923;

% cut out evaluation data
LW_sub_val_eval = vertcat(LW_in_bc_C_1h(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_C_1h(EP_C_2(1):EP_C_2(end),:));
LW_sub_val_eval_Sod = nan(length(LW_sub_val_eval(:,1)),1);
for l=1:length(LW_sub_val_eval_Sod)
    LW_sub_val_eval_Sod(l) = mean(LW_sub_val_eval(l,:));
end
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_CLM(EP_C_2(1):EP_C_2(end),:));
LW_sub_CLM_eval_Sod = nan(length(LW_sub_CLM_eval(:,1)),1);
for l=1:length(LW_sub_CLM_eval_Sod)
    LW_sub_CLM_eval_Sod(l) = mean(LW_sub_CLM_eval(l,:));
end
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_C_1(1):EP_C_1(end),:),...
    LW_in_bc_SP(EP_C_2(1):EP_C_2(end),:));
LW_sub_SP_eval_Sod = nan(length(LW_sub_SP_eval(:,1)),1);
for l=1:length(LW_sub_SP_eval_Sod)
    LW_sub_SP_eval_Sod(l) = mean(LW_sub_SP_eval(l,:));
end

% save data
save('ModelComp_LWsub_Sodankyla.mat','LW_sub_val_eval_Sod',...
    'LW_sub_CLM_eval_Sod','LW_sub_SP_eval_Sod')