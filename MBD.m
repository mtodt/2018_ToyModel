function  result = MBD(t_end,sim,obs)

step = nan(t_end,1);

for t=1:t_end
    step(t) = sim(t) - obs(t);
end

result = mean(step);