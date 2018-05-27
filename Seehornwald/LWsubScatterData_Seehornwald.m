% evaluation periods
EP_2008 = 1921:4728;
EP_2009 = 10705:12864;
EP_2010 = 19465:22080;
EP_2011 = 28225:30312;
EP_2012 = 36985:39768;

% gaps
gap_2010 = 1200:1224;   % 20 Feb 0:00 - 21 Feb 0:00
gap_2011 = 408:456;     % 18 Jan 0:00 - 20 Jan 0:00
gap_LWsub_2008 = 159;   % 7 Jan 15:00
gap_LWsub_2010 = 134:137;   % 6 Jan 14:00 - 17:00
gap_LWsub_2010_2 = 1200;    % 20 Feb 0:00
gap_LWsub_2011 = 346:349;       % 15 Jan 10:00 - 13:00
gap_LWsub_2011_2 = 427:442;     % 18 Jan 19:00 - 19 Jan 10:00

% cut out evaluation data
    % 2008
LW_sub_val_eval_SHW08 = vertcat(LW_in_bc_1h(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_1h(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_CLM_eval_SHW08 = vertcat(LW_in_bc_CLM(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_CLM(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
LW_sub_SP_eval_SHW08 = vertcat(LW_in_bc_SP(EP_2008(1):EP_2008(1)+gap_LWsub_2008-1-15),...
    LW_in_bc_SP(EP_2008(1)+gap_LWsub_2008-1+10:EP_2008(end)));
    % 2009
LW_sub_val_eval_SHW09 = LW_in_bc_1h(EP_2009(1):EP_2009(end));
LW_sub_CLM_eval_SHW09 = LW_in_bc_CLM(EP_2009(1):EP_2009(end));
LW_sub_SP_eval_SHW09 = LW_in_bc_SP(EP_2009(1):EP_2009(end));
    % 2010
LW_sub_val_eval_SHW10 = vertcat(LW_in_bc_1h(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_1h(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_1h(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_CLM_eval_SHW10 = vertcat(LW_in_bc_CLM(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_CLM(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_CLM(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
LW_sub_SP_eval_SHW10 = vertcat(LW_in_bc_SP(EP_2010(1):EP_2010(1)+gap_LWsub_2010(1)-1-14),...
    LW_in_bc_SP(EP_2010(1)+gap_LWsub_2010(end)-1+8:EP_2010(1)+gap_2010(1)-1-24),...
    LW_in_bc_SP(EP_2010(1)+gap_2010(end)-1+1:EP_2010(end)));
    % 2011
LW_sub_val_eval_SHW11 = vertcat(LW_in_bc_1h(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_1h(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_1h(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_CLM_eval_SHW11 = vertcat(LW_in_bc_CLM(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_CLM(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_CLM(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
LW_sub_SP_eval_SHW11 = vertcat(LW_in_bc_SP(EP_2011(1):EP_2011(1)+gap_LWsub_2011(1)-1-10),...
    LW_in_bc_SP(EP_2011(1)+gap_LWsub_2011(end)-1+12:EP_2011(1)+gap_2011(1)-1-24),...
    LW_in_bc_SP(EP_2011(1)+gap_2011(end)-1+1:EP_2011(end)));
    % 2012
LW_sub_val_eval_SHW12 = LW_in_bc_1h(EP_2012(1):EP_2012(end));
LW_sub_CLM_eval_SHW12 = LW_in_bc_CLM(EP_2012(1):EP_2012(end));
LW_sub_SP_eval_SHW12 = LW_in_bc_SP(EP_2012(1):EP_2012(end));

% concatenate years
LW_sub_val_eval_SHW = vertcat(LW_sub_val_eval_SHW08,LW_sub_val_eval_SHW09,...
    LW_sub_val_eval_SHW10,LW_sub_val_eval_SHW11,LW_sub_val_eval_SHW12);
LW_sub_CLM_eval_SHW = vertcat(LW_sub_CLM_eval_SHW08,LW_sub_CLM_eval_SHW09,...
    LW_sub_CLM_eval_SHW10,LW_sub_CLM_eval_SHW11,LW_sub_CLM_eval_SHW12);
LW_sub_SP_eval_SHW = vertcat(LW_sub_SP_eval_SHW08,LW_sub_SP_eval_SHW09,...
    LW_sub_SP_eval_SHW10,LW_sub_SP_eval_SHW11,LW_sub_SP_eval_SHW12);

% save data
save('ModelComp_LWsub_Seehornwald.mat','LW_sub_val_eval_SHW','LW_sub_CLM_eval_SHW','LW_sub_SP_eval_SHW')