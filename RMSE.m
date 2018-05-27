function  result = RMSE(t_end,sim,obs)

step = nan(t_end,1);

for t=1:t_end
    step(t) = (sim(t) - obs(t))^2;
end

result = sqrt(mean(step));