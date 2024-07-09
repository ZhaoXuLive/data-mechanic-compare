function out = Maaav_GPR_discrep_divide(Y, U, T, models, avp, sep, fitMode, numBlocks)

% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

dt = T(2) - T(1);
out = zeros(length(T), 10);

for iter = 1:length(sep(:, 1))
    out(sep(iter, 1), :) = Y(sep(iter, 1), :);
    disp(iter);
    for i = sep(iter, 1) + 1:sep(iter,2)
        
        if fitMode == 1
        	qdd = Maaav(out(i - 1, :), avp, U(i - 1, :))';
            out(i, 6:10) = out(i-1, 6:10) + qdd * dt;
            out(i, 1:5) = out(i-1, 1:5) + out(i-1, 6:10) * dt; 
            
            disc = zeros(numBlocks, 10);
            Ypre = zeros(1, 10);
            for j = 1:10
                Xpre = [out(i - 1, 3:10) U(i - 1, :)];
                for k = 1:numBlocks
                    disc(k, j) = predict(models{k}, Xpre);
                end
                Ypre(1, j) = mean(disc(:, j));
            end
            out(i, :) = out(i, :) + Ypre * dt;
        else
            disc = zeros(numBlocks, 10);
            Ypre = zeros(1, 10);
            for j = 1:10
                Xpre = [out(i - 1, 3:10) U(i - 1, :)];
                for k = 1:numBlocks
                    disc(k, j) = predict(models{k}, Xpre);
                end
                Ypre(1, j) = mean(disc(:, j));
            end
            out(i, :) = out(i-1, :) + Ypre;
        end
    end
end
end

