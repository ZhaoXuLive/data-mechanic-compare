function out = Maaav_DMD_discrep_fixb(Phi, Lambda, b, Y, U, T, avp, sep, fitMode)

dt = T(2) - T(1);
out = zeros(length(T), 10);

for iter = 1:length(sep(:, 1))
    out(sep(iter, 1), :) = Y(sep(iter, 1), :);
            
    for i = sep(iter, 1) + 1:sep(iter,2)
        
        if fitMode == 1
%         	qdd = Maaav(out(i - 1, :), avp, U(i - 1, :))';
%             disc = real(Phi * (Lambda^(i-1)) * b);
%             qdd = disc(6:10, :)' + qdd;
%             out(i, 6:10) = out(i-1, 6:10) + qdd * dt;
%             out(i, 1:5) = out(i-1, 1:5) + out(i-1, 6:10) * dt;
            qdd = Maaav(out(i - 1, :), avp, U(i - 1, :))';
            out(i, 6:10) = out(i-1, 6:10) + qdd * dt;
            out(i, 1:5) = out(i-1, 1:5) + out(i-1, 6:10) * dt;
            disc = real(Phi * (Lambda^(i-1)) * b);
            out(i, :) = out(i, :) + disc' * dt;
        end
        if fitMode == 2
            disc = real(Phi * (Lambda^(i-1)) * b);
            out(i, :) = disc';
        end
    end
end

end