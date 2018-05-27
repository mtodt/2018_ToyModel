function PAI = InterpolateMonthlyPAI(kyr,kmo,kda,PAImonthly)

if mod(kyr,4) == 0 && mod(kyr,100) ~= 0
    leapyear = 'yes';
elseif mod(kyr,4) == 0 && mod(kyr,100) == 0 && mod(kyr,400) == 0
    leapyear = 'yes';
elseif mod(kyr,4) == 0 && mod(kyr,100) == 0 && mod(kyr,400) ~= 0
    leapyear = 'no';
else
    leapyear = 'no';
end

if strcmp(leapyear,'no') == 1
    ndaypm = [31 28 31 30 31 30 31 31 30 31 30 31];
elseif strcmp(leapyear,'yes') == 1
    ndaypm = [31 29 31 30 31 30 31 31 30 31 30 31];
end

t = (kda-0.5)/ndaypm(kmo);
it(1) = floor(t + 0.5);
it(2) = floor(it(1) + 1);
months(1) = kmo + it(1) - 1;
months(2) = kmo + it(2) - 1;
if months(1) <  1
    months(1) = 12;
end
if months(2) > 12
    months(2) = 1;
end
timwt(1) = (it(1)+0.5) - t;
timwt(2) = 1-timwt(1);

mpai2t(1) = PAImonthly(months(1));
mpai2t(2) = PAImonthly(months(2));

PAI = timwt(1)*mpai2t(1) + timwt(2)*mpai2t(2);
end