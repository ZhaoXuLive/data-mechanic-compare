function out = Maaav_SINDy_discrep_multi(Y, U, T, Xi, Xintercept, avp, polyorder, sep, poolMode, fitMode, Sigvars)

% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

dt = T(2) - T(1);
out = zeros(length(T), 10);

% 用于记录差异数据和候选函数数据，辅助测试
eff = zeros(length(T), 5);
poolff = zeros(length(T), length(Xi(:, 1)));

for iter = 1:length(sep(:, 1))
    out(sep(iter, 1), :) = Y(sep(iter, 1), :);
    
    for i = sep(iter, 1) + 1:sep(iter,2)
        
        if fitMode == 1
        	qdd = Maaav(out(i - 1, :), avp, U(i - 1, :))';
        else
            qdd = zeros(1, 5);
        end
        
        if poolMode == 1
            yPool = poolData_withU(out(i - 1, 3:10), U(i - 1, :), polyorder);
        end
        if poolMode == 2
            yPool = poolData_new(out(i - 1, 3:10), U(i - 1, :), polyorder);
        end
        if poolMode == 3
            augData = zeros(1, 13);
            augData(1, 1:8) = out(i - 1, 3:10);
            if i-1 > sep(iter, 1)
                augData(1, 9:13) = out(i - 1, 6:10) - out(i - 2, 6:10);
            end
            yPool = poolData_augNew(augData, U(i - 1, :), polyorder);
        end
        
        disc = zeros(1,5);
        for j = 1:5
            yPool_sig = yPool(:, Sigvars(2:Sigvars(1,j)+1, j));
            disc(j) = yPool_sig * Xi(2:Xi(1,j)+1, j) + Xintercept(j);
        end
%         disc = yPool * Xi + Xintercept'; 
        eff(i,:) = disc;
        poolff(i,:) = yPool;
        qdd = disc + qdd;

        out(i, 6:10) = out(i-1, 6:10) + qdd * dt;
        out(i, 1:5) = out(i-1, 1:5) + out(i-1, 6:10) * dt;

        % when predict over and over, do some limit
        if(hypot(out(i, 1) - out(i-1, 1), out(i, 2) - out(i-1, 2)) > 0.1)
            break;
%             disp(i);
%             out(i, 6:10) = out(i-1, 6:10) + (qdd - disc) * dt;
%             out(i, 1:5) = out(i-1, 1:5) + out(i-1, 6:10) * dt;
        end
    end
end
end

