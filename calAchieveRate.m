function r = calAchieveRate(h, p, omega, sigma)
%% 计算achievable rate
r = omega*log2(1+h*p/sigma);

end