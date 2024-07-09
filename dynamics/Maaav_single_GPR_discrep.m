function out = Maaav_single_GPR_discrep(Y, U, T, discrepGPR, avp, sep, fitMode)

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
        
        disc = zeros(1, 5);
        for j = 1:5
            Xpre = [Y(i - 1, :) U(i - 1, :)];
            disc(1, j) = predict(discrepGPR{j}, Xpre);
        end
        qdd = disc + qdd;

        out(i, 6:10) = Y(i-1, 6:10) + qdd * dt;
        out(i, 1:5) = Y(i-1, 1:5) + out(i-1, 6:10) * dt;
    end
end

end
