function [day,month] = DoY_to_DayNMonth(year,doy)
%{
Pretty self-explanatory.
%}

if mod(year,400) == 0
    leapyear = 1;
elseif mod(year,4) == 0 && mod(year,100) ~= 0
    leapyear = 1;
else
    leapyear = 0;
end

if leapyear == 1
    months = [31 29 31 30 31 30 31 31 30 31 30 31];
elseif leapyear == 0
    months = [31 28 31 30 31 30 31 31 30 31 30 31];
end

for i=1:12
    bla = doy - sum(months(1:i));
    if bla < 0
        month = i;
        if i > 1
            day = doy - sum(months(1:i-1));
        else
            day = doy;
        end
        break
    elseif bla == 0
        month = i;
        day = months(i);
        break
    end
end

end