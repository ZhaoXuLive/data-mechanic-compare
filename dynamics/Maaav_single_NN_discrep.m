function out = Maaav_single_NN_discrep(Y, U, T, net, avp, sep, fitMode)

% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

dt = T(2) - T(1);
out = zeros(length(T), 10);

for iter = 1:length(sep(:, 1))
    out(sep(iter, 1), :) = Y(sep(iter, 1), :);
    
    for i = sep(iter, 1)+1:sep(iter, 2)
        
%         if fitMode == 1
%             qdd = Maaav(Y(i - 1, :), avp, U(i - 1, :))';
%         else
%             qdd = zeros(1, 5);
%         end
%         disc = net([Y(i - 1, :) U(i - 1, :)]')'; 
%         qdd = disc + qdd;
%         out(i, 6:10) = Y(i-1, 6:10) + qdd * dt;
%         out(i, 1:5) = Y(i-1, 1:5) + out(i-1, 6:10) * dt;

        if fitMode == 1
            qdd = Maaav(Y(i - 1, :), avp, U(i - 1, :))';
            out(i, 6:10) = Y(i-1, 6:10) + qdd * dt;
            out(i, 1:5) = Y(i-1, 1:5) + out(i-1, 6:10) * dt;
            disc = net([Y(i - 1, 3:10) U(i - 1, :)]')';
            out(i, :)  = out(i, :) + disc * dt;
        else
            disc = net([Y(i - 1, 3:10) U(i - 1, :)]')';
            out(i, :) = disc;
        end
        
    end
end

end
