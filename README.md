# 2018_ToyModel
Repository containing Matlab code used to run and analyze Toy Model simulations and create figures shown in Todt et al. (2018).
Instructions for each site with site name xxxxxxxare the following. To run Toy Model simulations and save output:
1) Run DataPrep_xxxxxxx.m.
2) Run ToyModelAtxxxxx.m with HM_in_CLM set to 'yes'.
3) Run MLRdata_xxxxxxx_HM.m.
4) Run ToyModelAtxxxxx.m again with HM_in_CLM set to 'no'.
5) Run MLRdata_xxxxxxx.m and LWsubScatterData_xxxxxxx.m for Alptal, Seehornwald, and Sodankylä.
6) Run Metrics_xxxxxxx.m to calculate metrics in Table 3.

To create figures:
7) Run MetOverviewPlots_Paper.m to get Figures 3 and 4.
8) Run LWsubCompSubplot.m to create Figure 5.
9) Run PaperPlots_AllSites.m for Figures 6-9.

For sensitivity tests in Supplementary Material run SensitivityTest_xxxxxxx.m.
