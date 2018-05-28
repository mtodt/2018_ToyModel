%{
1) Run DataPrep_Alptal.m.
2) Run ToyModelAtAlptal.m with HM_in_CLM set to 'yes'.
3) Change LW_in_bc_CLM to LW_in_bc_CLM_HM and T_veg_CLM to T_veg_CLM_HM,
   and save both in CLM_Alptal_HM.mat.
        LW_in_bc_CLM_HM = LW_in_bc_CLM;
        T_veg_CLM_HM = T_veg_CLM;
        save('CLM_Alptal_HM.mat','LW_in_bc_CLM_HM','T_veg_CLM_HM')
4) Run ToyModelAtAlptal.m again with HM_in_CLM set to 'no'.
%}
load CLM_Alptal_HM.mat

%-------------------------------------------------------------------------%
%----------------------  select evaluation periods  ----------------------%
%{
Full days are selected for evaluation, so that diurnal cycles can be
calculated and diurnal biases prevented generally. Meltout can be determined
from outgoing sub-canopy longwave radiation and soil temperature (at 20cm
depth) and is used to limit evaluation period. Limits for each year are:
- 2004
    17 April 17:00, although there seems to be an earlier snow-free period
    from 19 March 16:00 to 2 April 6:00.
    Manual snow measurements still have more than 10cm for 18 March
    although LWR features clear diurnal cycles with peaks of >340 W m^{-2},
    which might be tied to understory vegetation or patchy surface cover.
Limit evaluation period to 12 March 0:00 as period afterwards would
introduce uncertainty (surface temperature, soil moisture, snow depth,...)
and still gives us 48 days.
- 2005
    This is as easy as it gets, there is a clear change in LWR for 14 March
    12:00 to distinct diurnal cycles of 320 to 340 W m^{-2}. Again, manual
    snow depth measurements still give snow coverage, actually about 50cm
    while one week later there's only 15cm left, and there is another
    change in LWR for 30 March 2005 - but this would introduce uncertainty.
Limit evaluation period to 14 March 0:00.
- 2006
    There's a clear change in LWR from 18 March to 19 March, so we set the
end of evaluation to 19 March 0:00. Soil temperature also displays a
    decrease of about 0.5°C in the days after 19 March indicating re-coupling
    to the air and/or drainage of melt water.
    
- 2007
    This is a hard one again, as LWR display almost consitently diurnal
    cycles and regularly reaches 330-340 W m^{-2}. The only clear
    indication is an increase in the range of diurnal cycles from 5 April
    to 6 April (and the following handful of days). Soil temperature
    strengthens the case as there are extensive warming and cooling phases
    indicating periods of little to no snow cover. Manual measurements also
    show no snow depth larger than 12 cm until late March, although there
    is a gap of 3 weeks. This season can be used as a test.
Limit evaluation period to 6 April 0:00.
Update: The last manual snow measurement took place on 4 April 2007, so
evaluation period is limited to 4 April 2007 0:00.

Generally, it seems that LWR is distinctly larger than 320 W m^{-2} for the
final days of snowmelt indicating snow cover is already patchy and/or
understory vegetation plays a major role. To reduce uncertainty, these
periods are not included in evaluation.
Start of each evaluation period is 1 January 1:00 (that's arbitrary),
except for 2004 when radiation measurements start 23 January 7:00 and
evaluation period consequently starts 24 January 1:00.
%}
EP_2004 = start_2004+18:1728;   % 48 days
EP_2005 = 25:1752;              % 72 days
EP_2006 = 25:1872;              % 77 days
EP_2007 = 25:2256;              % 93 days

% consider gaps of evaluation data (snow on radiometer)
% 2004
gap_2004_1 = 601:624;           % 25 January 2004
gap_2004_2 = 1009:1032;         % 11 February 2004
gap_2004_3 = 1201:1224;         % 19 February 2004
gap_2004_4 = 1633:1728;         % 8 - 11 March 2004
% 2005
gap_2005_1 = 553:624;           % 23 - 25 January 2005
gap_2005_2 = 769:792;           % 1 February 2005
gap_2005_3 = 1057:1272;         % 13 - 21 February 2005
gap_2005_4 = 1609:1656;         % 8 - 9 March 2005
% 2006
gap_2006 = 937:1032;            % 8 - 11 February 2006           
% 2007
gap_2007_1 = 49:96;             % 2 - 3 January 2007
gap_2007_2 = 1873:2016;         % 19 - 24 March 2007


%-------------------------------------------------------------------------%
%---------------------  analyse: figures, RMSE, MBD  ---------------------%

%------------------------------  year 2004  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h_all(EP_2004(1):gap_2004_1(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
    LW_in_bc_1h_all(gap_2004_4(end)+1:EP_2004(end),1));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2004(1):gap_2004_1(1)-1,1),...
    LW_in_bc_CLM(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
    LW_in_bc_CLM(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
    LW_in_bc_CLM(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
    LW_in_bc_CLM(gap_2004_4(end)+1:EP_2004(end),1));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2004(1):gap_2004_1(1)-1,1),...
    LW_in_bc_CLM_HM(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
    LW_in_bc_CLM_HM(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
    LW_in_bc_CLM_HM(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
    LW_in_bc_CLM_HM(gap_2004_4(end)+1:EP_2004(end),1));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2004(1):gap_2004_1(1)-1,1),...
    LW_in_bc_SP(gap_2004_1(end)+1:gap_2004_2(1)-1,1),...
    LW_in_bc_SP(gap_2004_2(end)+1:gap_2004_3(1)-1,1),...
    LW_in_bc_SP(gap_2004_3(end)+1:gap_2004_4(1)-1,1),...
    LW_in_bc_SP(gap_2004_4(end)+1:EP_2004(end),1));

RMSE_CLM_2004 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2004 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2004 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2004 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2004 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2004 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

LW_sub_val_eval_all = LW_sub_val_eval;
LW_sub_CLM_eval_all = LW_sub_CLM_eval;
LW_sub_CLMHM_eval_all = LW_sub_CLMHM_eval;
LW_sub_SP_eval_all = LW_sub_SP_eval;

%------------------------------  year 2005  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h_all(EP_2005(1):gap_2005_1(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_1(end)+1:gap_2005_2(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_2(end)+1:gap_2005_3(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_3(end)+1:gap_2005_4(1)-1,2),...
    LW_in_bc_1h_all(gap_2005_4(end)+1:EP_2005(end),2));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2005(1):gap_2005_1(1)-1,2),...
    LW_in_bc_CLM(gap_2005_1(end)+1:gap_2005_2(1)-1,2),...
    LW_in_bc_CLM(gap_2005_2(end)+1:gap_2005_3(1)-1,2),...
    LW_in_bc_CLM(gap_2005_3(end)+1:gap_2005_4(1)-1,2),...
    LW_in_bc_CLM(gap_2005_4(end)+1:EP_2005(end),2));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2005(1):gap_2005_1(1)-1,2),...
    LW_in_bc_CLM_HM(gap_2005_1(end)+1:gap_2005_2(1)-1,2),...
    LW_in_bc_CLM_HM(gap_2005_2(end)+1:gap_2005_3(1)-1,2),...
    LW_in_bc_CLM_HM(gap_2005_3(end)+1:gap_2005_4(1)-1,2),...
    LW_in_bc_CLM_HM(gap_2005_4(end)+1:EP_2005(end),2));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2005(1):gap_2005_1(1)-1,2),...
    LW_in_bc_SP(gap_2005_1(end)+1:gap_2005_2(1)-1,2),...
    LW_in_bc_SP(gap_2005_2(end)+1:gap_2005_3(1)-1,2),...
    LW_in_bc_SP(gap_2005_3(end)+1:gap_2005_4(1)-1,2),...
    LW_in_bc_SP(gap_2005_4(end)+1:EP_2005(end),2));

RMSE_CLM_2005 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2005 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2005 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2005 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2005 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2005 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

temp = LW_sub_CLM_eval_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_all = vertcat(temp,LW_sub_CLM_eval);
temp = LW_sub_CLMHM_eval_all; clear LW_sub_CLMHM_eval_all
    LW_sub_CLMHM_eval_all = vertcat(temp,LW_sub_CLMHM_eval);
temp = LW_sub_SP_eval_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_all = vertcat(temp,LW_sub_SP_eval);

%------------------------------  year 2006  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h_all(EP_2006(1):gap_2006(1)-1,3),...
    LW_in_bc_1h_all(gap_2006(end)+1:EP_2006(end),3));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2006(1):gap_2006(1)-1,3),...
    LW_in_bc_CLM(gap_2006(end)+1:EP_2006(end),3));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2006(1):gap_2006(1)-1,3),...
    LW_in_bc_CLM_HM(gap_2006(end)+1:EP_2006(end),3));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2006(1):gap_2006(1)-1,3),...
    LW_in_bc_SP(gap_2006(end)+1:EP_2006(end),3));

RMSE_CLM_2006 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2006 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2006 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2006 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2006 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2006 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

temp = LW_sub_val_eval_all; clear LW_sub_val_eval_all
    LW_sub_val_eval_all = vertcat(temp,LW_sub_val_eval);
temp = LW_sub_CLM_eval_all; clear LW_sub_CLM_eval_all
    LW_sub_CLM_eval_all = vertcat(temp,LW_sub_CLM_eval);
temp = LW_sub_CLMHM_eval_all; clear LW_sub_CLMHM_eval_all
    LW_sub_CLMHM_eval_all = vertcat(temp,LW_sub_CLMHM_eval);
temp = LW_sub_SP_eval_all; clear LW_sub_SP_eval_all
    LW_sub_SP_eval_all = vertcat(temp,LW_sub_SP_eval);

%------------------------------  year 2007  ------------------------------%
LW_sub_val_eval = vertcat(LW_in_bc_1h_all(EP_2007(1):gap_2007_1(1)-1,4),...
    LW_in_bc_1h_all(gap_2007_1(end)+1:gap_2007_2(1)-1,4),...
    LW_in_bc_1h_all(gap_2007_2(end)+1:EP_2007(end),4));
LW_sub_CLM_eval = vertcat(LW_in_bc_CLM(EP_2007(1):gap_2007_1(1)-1,4),...
    LW_in_bc_CLM(gap_2007_1(end)+1:gap_2007_2(1)-1,4),...
    LW_in_bc_CLM(gap_2007_2(end)+1:EP_2007(end),4));
LW_sub_CLMHM_eval = vertcat(LW_in_bc_CLM_HM(EP_2007(1):gap_2007_1(1)-1,4),...
    LW_in_bc_CLM_HM(gap_2007_1(end)+1:gap_2007_2(1)-1,4),...
    LW_in_bc_CLM_HM(gap_2007_2(end)+1:EP_2007(end),4));
LW_sub_SP_eval = vertcat(LW_in_bc_SP(EP_2007(1):gap_2007_1(1)-1,4),...
    LW_in_bc_SP(gap_2007_1(end)+1:gap_2007_2(1)-1,4),...
    LW_in_bc_SP(gap_2007_2(end)+1:EP_2007(end),4));

RMSE_CLM_2007 = RMSE(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
MBD_CLM_2007 = MBD(length(LW_sub_val_eval),LW_sub_CLM_eval,LW_sub_val_eval);
RMSE_CLMHM_2007 = RMSE(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
MBD_CLMHM_2007 = MBD(length(LW_sub_val_eval),LW_sub_CLMHM_eval,LW_sub_val_eval);
RMSE_SP_2007 = RMSE(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);
MBD_SP_2007 = MBD(length(LW_sub_val_eval),LW_sub_SP_eval,LW_sub_val_eval);

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

