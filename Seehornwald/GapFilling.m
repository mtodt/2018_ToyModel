function [tg] = GapFilling(tg,tp)
%{
Function to repair time series that contain gaps (NaNs). Gaps are filled
via interpolation from proxy time series. For example at Seehornwald, gaps
in time series of tower measurements are filled with data from FluxNet
(which should theoretically be identical but aren't for whatever reason).
Since both time series should be similar as they represent the same
variable for the same location, the filled data points are constrained to
represent roughly the same variations and not diverge too much from the
proxy.
Variables and parameters:
tg                      time series with gaps
tp                      proxy time series
r                       scaling ratio between gap and proxy
%}

if length(tg) ~= length(tp)
   disp('Input time series must be the same length.')
end
tl = length(tg);

for t=1:tl
    if isnan(tg(t))
        
% 1) determine gap length
        gap_start = t;
        gap = find(~isnan(tg(t:end)),1)-1;
        gap_end = t+gap-1;  % -1 because gap includes t
        
% 2) determine scaling ratio between gap and proxy
        r = (tg(gap_end+1)-tg(gap_start-1))/(tp(gap_end+1)-tp(gap_start-1));
        r = abs(r);
        if r < 2/3j
            r = 2/3;
        elseif r > 1.5
            r = 1.5;
        end
        
% 3) fill by incrementally applying scaled gradients from proxy
        for g=gap_start:gap_end
            if ~isnan(tp(g))
                tg(g) = tg(g-1) + r*(tp(g)-tp(g-1));
            else
                tg(g:gap_end) = nan;
                break
            end
        end
        
% 4) limit differences between filled data and proxy data
        maxdiff = max(tg(gap_start-1)-tp(gap_start-1),tg(gap_end+1)-tp(gap_end+1));
        mindiff = min(tg(gap_start-1)-tp(gap_start-1),tg(gap_end+1)-tp(gap_end+1));
        for g=gap_start:gap_end
            if ~isnan(tg(g))
                tg(g) = min(tg(g),tp(g)+maxdiff);
                tg(g) = max(tg(g),tp(g)+mindiff);
            end
        end
        
    end
end









end