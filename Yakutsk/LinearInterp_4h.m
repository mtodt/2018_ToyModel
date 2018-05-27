function bla = LinearInterp_4h(bla)
%{
Linear interpolation of hourly values for gaps of 4 hours or less.
%}
for b=2:length(bla)
    if isnan(bla(b)) && ~isnan(bla(b-1))
        
% 1) determine gap length
        gap_start = b;
        gap = find(~isnan(bla(b:end)),1)-1;
        gap_end = b+gap-1;  % -1 because gap includes b
        
% 2) fill by incrementally applying gradient if less than 4 hours
        if gap <= 4
            for g=gap_start:gap_end
                bla(g) = bla(g-1) + (bla(gap_end+1)-bla(gap_start-1))/(gap+1);
            end
        end
        
    end
end
end