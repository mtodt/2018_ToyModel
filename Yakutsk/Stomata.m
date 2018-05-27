function [rs,psn] = Stomata(PFT,forc_pbot,thm,btran,apar,sla,ei,ea,o2,...
    co2,rb,dayl_factor,tl)
%{
Leaf stomatal resistance and leaf photosynthesis. Modifications for CN code.
%}

tgcm = thm; % air temperature at agcm reference height [K]

% constants
mpe = 10^(-6);
kc25 = 30;
akc = 2.1;
ko25 = 30000;
ako = 1.2;
avcmx = 2.4;
bp = 2000;
act25 = 3.6;
q10act = 2.4;
fnr = 7.16;

% PFT-dependent constants
    % quantum efficiency at 25°C [umol CO2 / umol photon]
qe25 = [0 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 0.06 ...
    0.06 0.04 0.06 0.06 0.04 0.06 0.06 0.06];
    % leaf C:N (gC/gN)
leafcn = [1 35 40 25 30 30 25 25 25 30 25 25 25 25 25 25 25 25 25 25 25];
    % fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
    % values in Tech Notes and pft-physiology.c110425.nc don't match !!!
flnr = [0 0.05 0.04 0.08 0.06 0.06 0.09 0.09 0.09 0.06 0.09 0.09 0.09 ...
    0.09 0.09 0.1 0.1 0.1 0.2 0.2 0.1];
    % foliage nitrogen limitation factor
fnitr = [0 0.72 0.78 0.79 0.83 0.71 0.66 0.64 0.7 0.62 0.6 0.76 0.68 ...
    0.61 0.64 0.61 0.61 0.55 0.55 0.55 0.55 ];
    % photosynthetic pathway: 0 = c4, 1 = c3
c3psn = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 0 1 1 1];
    % slope of conductance-to-photosynthesis relationship
mp = [9 6 6 6 9 9 9 9 9 9 9 9 9 9 5 9 9 4 9 9 9];

% convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
act25 = act25*1000/60;

% initialsie rs=rsmax and psn=0 because calculations are performed only
% when apar>0, in which case rs<=rsmax and psn>=0
rsmax0 = 20000;
cf = forc_pbot/(8.3145*10^3 * 0.001*tgcm) * 10^6;
if apar <= 0    % night time
    rs = min(rsmax0,cf/bp);
    psn = 0;
    lnc = 0;
    vcmx = 0;
else
    tc = tl - 273.15;
    ppf = 4.6*apar;
    j = ppf*qe25(PFT);
    kc = kc25*(akc^((tc-25)/10));
    ko = ko25*(ako^((tc-25)/10));
    awc = kc*(1+o2/ko);
    cp = 0.5*kc/ko*o2*0.21;
    
    lnc = 1/(sla*leafcn(PFT));
    act = act25*(q10act^((tc-25)/10));
    f2 = 1 + exp((-2.2*10^5 + 710*(tc+273.15))/(8.3145*10^3 * 0.001*(tc+273.15)));
    vcmx = lnc*flnr(PFT)*fnr*act/f2*btran*dayl_factor*fnitr(PFT);
    
    % first guess ci
    ci = 0.7*co2*c3psn(PFT) + 0.4*co2*(1-c3psn(PFT));
    
    % rb: s m^{-1} -> s m^2 mol^{-1}
    rb = rb/cf;
    
    % constrain ea
    cea = max(0.25*ei*c3psn(PFT) + 0.4*ei*(1-c3psn(PFT)),min(ea,ei));
    
    % ci iteration for 'actual' photosynthesis !!! not NEC_SX !!!
    for iter=1:3
       wj = max(ci-cp,0)*j/(ci+2*cp)*c3psn(PFT) + j*(1-c3psn(PFT)); 
       wc = max(ci-cp,0)*vcmx/(ci+awc)*c3psn(PFT) + vcmx*(1-c3psn(PFT)); 
       we = 0.5*vcmx*c3psn(PFT) + 4000*vcmx*ci/forc_pbot * (1-c3psn(PFT)); 
       psn = min(min(wj,wc),we);
       cs = max(co2 - 1.37*rb*forc_pbot*psn,mpe);
       atmp = mp(PFT)*psn*forc_pbot*cea/(cs*ei) + bp;
       btmp = (mp(PFT)*psn*forc_pbot/cs + bp)*rb - 1;
       ctmp = -rb;
       if btmp >= 0
           q = -0.5*(btmp + sqrt(btmp^2 - 4*atmp*ctmp));
       else
           q = -0.5*(btmp - sqrt(btmp^2 - 4*atmp*ctmp));
       end
       r1 = q/atmp;
       r2 = ctmp/q;
       rs = max(r1,r2);
       ci = max(cs - psn*forc_pbot*1.65*rs,0);
    end

    % rs, rb: s m^2 umol^{-1} -> s m^{-1}
    rs = min(rsmax0,rs*cf);
    rb = rb*cf;
end
    
end