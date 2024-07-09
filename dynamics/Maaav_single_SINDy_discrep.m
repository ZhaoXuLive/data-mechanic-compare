function out = Maaav_single_SINDy_discrep(Y, U, T, Xi, Xintercept, avp, polyorder, sep, poolMode, fitMode)

% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

dt = T(2) - T(1);
out = zeros(length(T), 10);

for iter = 1:length(sep(:, 1))
    out(sep(iter, 1), :) = Y(sep(iter, 1), :);
    
    for i = sep(iter, 1)+1:sep(iter, 2)
        if fitMode == 1
            qdd = Maaav(Y(i - 1, :), avp, U(i - 1, :))';
        else
            qdd = zeros(1, 5);
        end
        if poolMode == 1
            yPool = poolData_withU(Y(i - 1, 3:10), U(i - 1, :), polyorder);
        end
        if poolMode == 2
            yPool = poolData_new(Y(i - 1, 3:10), U(i - 1, :), polyorder);
        end
        if poolMode == 3
            augData = zeros(1, 13);
            augData(1, 1:8) = Y(i - 1, 3:10);
            if i-1 > sep(iter, 1)
                augData(1, 9:13) = Y(i - 1, 6:10) - Y(i - 2, 6:10);
            end
            yPool = poolData_augNew(augData, U(i - 1, :), polyorder);
        end
        disc = yPool * Xi(:, 1:5) + Xintercept'; 
        qdd = disc + qdd;

        out(i, 6:10) = Y(i-1, 6:10) + qdd * dt;
        out(i, 1:5) = Y(i-1, 1:5) + out(i-1, 6:10) * dt;
    end
end

end
