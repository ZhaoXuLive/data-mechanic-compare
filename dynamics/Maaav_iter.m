 function out = Maaav_iter(tspan, y, avp, u, n)
% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

dt = tspan(2) -  tspan(1);
x0 = y(1, :);  

out = zeros(length(tspan), n);
out(1, :) = x0;

for i = 2:length(tspan)
    out(i, 6:10) = out(i-1, 6:10) + Maaav(out(i - 1, :), avp, u(i - 1, :))' * dt;
    out(i, 1:5) = out(i-1, 1:5) + out(i-1, 6:10) * dt;
end