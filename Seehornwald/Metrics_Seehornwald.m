%{
1) Run DataPrep_Seehornwald.m.
2) Run ToyModelAtSeehornwald.m with HM_in_CLM set to 'yes'.
3) Change LW_in_bc_CLM to LW_in_bc_CLM_HM and T_veg_CLM to T_veg_CLM_HM,
   and save both in CLM_Seehornwald_HM.mat.
        LW_in_bc_CLM_HM = LW_in_bc_CLM;
        T_veg_CLM_HM = T_veg_CLM;
        save('CLM_Seehornwald_HM.mat','LW_in_bc_CLM_HM','T_veg_CLM_HM')
4) Run ToyModelAtSeehornwald.m again with HM_in_CLM set to 'no'.
%}
load CLM_Seehornwald_HM.mat

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
%{
Full days are selected for evaluation, so that diurnal cycles can be
calculated and diurnal biases prevented generally. Meltout can be determined
from outgoing sub-canopy longwave radiation and is used to limit evaluation
period. Limits for each year are:
27 April 2008 0:00 (by eye: 26 April 2008 8:00)
1 April 2009 0:00
20 April 2010 0:00 (this one is hard)
29 March 2011 0:00 (this one is hard as well)
26 April 2012 0:00
Snow depth gives different, later results. However, since this is scaled
down from the open, we use earlier dates to prevent introducing
uncertainty. Results from snow depth are:
05 May 2008 0:00 (4920)
24 April 2009 0:00 (13416)
21 April 2010 0:00 (22104)
05 April 2011 0:00 (30480)
02 May 2012 0:00 (39912)
The difference for 2009 likely stems from the rapid, almost linear meltout,
so that measured LWR should be a combination of patches of soil and snow.
Start of each evaluation period is 1 January 1:00 (that's arbitrary).
%}
EP_2008 = 1921:4728;
EP_2009 = 10705:12864;
EP_2010 = 19465:22080;
EP_2011 = 28225:30312;
EP_2012 = 36985:39768;

%-------------------------------------------------------------------------%
%--------------  find & skip gaps and uncertain time steps  --------------%
% due to met forcing
find(isnan(EvalSeehornwald(EP_2008(1):EP_2008(end))));
find(isnan(EvalSeehornwald(EP_2009(1):EP_2009(end))));
find(isnan(EvalSeehornwald(EP_2010(1):EP_2010(end))));
gap_2010 = 1200:1224;   % 20 Feb 0:00 - 21 Feb 0:00
find(isnan(EvalSeehornwald(EP_2011(1):EP_2011(end))));
gap_2011 = 408:456;     % 18 Jan 0:00 - 20 Jan 0:00
find(isnan(EvalSeehornwald(EP_2012(1):EP_2012(end))));
%{
Due to the gaps lasting from 0:00 to 0:00, we need to exclude the entire
previous day as well. Effectively, gaps are therefore 19 February 2010 1:00
to 21 February 2010 0:00 and 17 January 2011 1:00 to 20 January 2011 0:00.
%}
% due to evaluation data
find(isnan(LW_in_bc_1h(EP_2008(1):EP_2008(end))));
gap_LWsub_2008 = 159;   % 7 Jan 15:00
find(isnan(LW_in_bc_1h(EP_2009(1):EP_2009(end))));
find(isnan(LW_in_bc_1h(EP_2010(1):EP_2010(end))));
gap_LWsub_2010 = 134:137;   % 6 Jan 14:00 - 17:00
gap_LWsub_2010_2 = 1200;    % 20 Feb 0:00
find(isnan(LW_in_bc_1h(EP_2011(1):EP_2011(end))));
gap_LWsub_2011 = 346:349;       % 15 Jan 10:00 - 13:00
gap_LWsub_2011_2 = 427:442;     % 18 Jan 19:00 - 19 Jan 10:00
find(isnan(LW_in_bc_1h(EP_2012(1):EP_2012(end))));
%{
Second gap for 2010 and 2011 each already covered by forcing data.
%}

%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%

%------------------------------  year 2008  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_1h(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_CLM(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_CLM_HM(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_SP(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));

RMSE_CLM_2008 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2008 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2008 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2008 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2008 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2008 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

LW_sub_val_eval_all = LW_sub_val_eval;
LW_sub_CLM_eval_all = LW_sub_CLM_eval;
LW_sub_CLMHM_eval_all = LW_sub_CLMHM_eval;
LW_sub_SP_eval_all = LW_sub_SP_eval;

%------------------------------  year 2009  ------------------------------%
LW_sub_val_eval = LW_in_bc_1h(EP_2009(1):EP_2009(end));
LW_sub_CLM_eval = LW_in_bc_CLM(EP_2009(1):EP_2009(end));
LW_sub_CLMHM_eval = LW_in_bc_CLM_HM(EP_2009(1):EP_2009(end));
LW_sub_SP_eval = LW_in_bc_SP(EP_2009(1):EP_2009(end));

RMSE_CLM_2009 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2009 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2009 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2009 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2009 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2009 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
temp = LW_sub_CLM_eval_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_all = vertcat(temp,LW_sub_CLM_eval);
temp = LW_sub_CLMHM_eval_all; clear LW_sub_CLMHM_eval_all
    LW_sub_CLMHM_eval_all = vertcat(temp,LW_sub_CLMHM_eval);
temp = LW_sub_SP_eval_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_all = vertcat(temp,LW_sub_SP_eval);

%------------------------------  year 2010  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_1h(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_1h(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_CLM(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_CLM(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_CLM_HM(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_CLM_HM(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_SP(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_SP(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));

RMSE_CLM_2010 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2010 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2010 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2010 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2010 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2010 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
temp = LW_sub_CLM_eval_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_all = vertcat(temp,LW_sub_CLM_eval);
temp = LW_sub_CLMHM_eval_all; clear LW_sub_CLMHM_eval_all
    LW_sub_CLMHM_eval_all = vertcat(temp,LW_sub_CLMHM_eval);
temp = LW_sub_SP_eval_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_all = vertcat(temp,LW_sub_SP_eval);

%------------------------------  year 2011  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_1h(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_1h(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_CLM(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_CLM(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_CLM_HM(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_CLM_HM(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_SP(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_SP(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));

RMSE_CLM_2011 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2011 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2011 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2011 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2011 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2011 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
temp = LW_sub_CLM_eval_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_all = vertcat(temp,LW_sub_CLM_eval);
temp = LW_sub_CLMHM_eval_all; clear LW_sub_CLMHM_eval_all
    LW_sub_CLMHM_eval_all = vertcat(temp,LW_sub_CLMHM_eval);
temp = LW_sub_SP_eval_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_all = vertcat(temp,LW_sub_SP_eval);

%------------------------------  year 2012  ------------------------------%
LW_sub_val_eval = LW_in_bc_1h(EP_2012(1):EP_2012(end));
LW_sub_CLM_eval = LW_in_bc_CLM(EP_2012(1):EP_2012(end));
LW_sub_CLMHM_eval = LW_in_bc_CLM_HM(EP_2012(1):EP_2012(end));
LW_sub_SP_eval = LW_in_bc_SP(EP_2012(1):EP_2012(end));

RMSE_CLM_2012 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2012 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2012 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2012 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2012 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2012 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
temp = LW_sub_CLM_eval_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_all = vertcat(temp,LW_sub_CLM_eval);
temp = LW_sub_CLMHM_eval_all; clear LW_sub_CLMHM_eval_all
    LW_sub_CLMHM_eval_all = vertcat(temp,LW_sub_CLMHM_eval);
temp = LW_sub_SP_eval_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_all = vertcat(temp,LW_sub_SP_eval);

%------------------------------  all years  ------------------------------%
RMSE_CLM_all = RMSE(length(LW_sub_val_eval_all),LW_sub_CLM_eval_all,LW_sub_val_eval_all);
MBD_CLM_all = MBD(length(LW_sub_val_eval_all),LW_sub_CLM_eval_all,LW_sub_val_eval_all);
RMSE_CLMHM_all = RMSE(length(LW_sub_val_eval_all),LW_sub_CLMHM_eval_all,LW_sub_val_eval_all);
MBD_CLMHM_all = MBD(length(LW_sub_val_eval_all),LW_sub_CLMHM_eval_all,LW_sub_val_eval_all);
RMSE_SP_all = RMSE(length(LW_sub_val_eval_all),LW_sub_SP_eval_all,LW_sub_val_eval_all);
MBD_SP_all = MBD(length(LW_sub_val_eval_all),LW_sub_SP_eval_all,LW_sub_val_eval_all);
