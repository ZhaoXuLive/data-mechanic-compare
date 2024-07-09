function out = Maaav_single_DMD_discrep(Phi, Lambda, Y, U, T, avp, sep, fitMode)

dt = T(2) - T(1);
out = zeros(length(T), 10);


for iter = 1:length(sep(:, 1))
    out(sep(iter, 1), :) = Y(sep(iter, 1), :);
    
    for i = sep(iter, 1) + 1:sep(iter,2)
        
%         if fitMode == 1
%         	qdd = Maaav(Y(i - 1, :), avp, U(i - 1, :))';
%         else
%             qdd = zeros(1, 5);
%         end
%         disc = zeros(1, 10);
%         if i - sep(iter, 1) > 1
%             lastDY = Y(i, :) - Y(i-1, :);
%             lastX1 = zeros(1, 10);
%             lastX1(6:10) = Y(i-1, 6:10) + Maaav(Y(i - 1, :), avp, U(i - 1, :))' * dt;
%             lastX1(1:5) = Y(i-1, 1:5) + lastX1(6:10) * dt;
%             lastX = zeros(1, 10);
%             lastX(6:10) = Y(i, 6:10) + Maaav(Y(i, :), avp, U(i, :))' * dt;
%             lastX(1:5) = Y(i, 1:5) + lastX(6:10) * dt;
%             lastDX = lastX - lastX1;
%             lastEf = (lastDY - lastDX) / dt;
%             b = Phi \ lastEf';
%             disc = real(Phi * Lambda * b);
%             qdd = disc(6:10, :)' + qdd;
%         end
%         out(i, 6:10) = Y(i-1, 6:10) + qdd * dt;
%         out(i, 1:5) = Y(i-1, 1:5) + out(i-1, 6:10) * dt;

        if fitMode == 1
            qdd = Maaav(Y(i - 1, :), avp, U(i - 1, :))';
            out(i, 6:10) = Y(i-1, 6:10) + qdd * dt;
            out(i, 1:5) = Y(i-1, 1:5) + out(i-1, 6:10) * dt;
            disc = zeros(1, 10);
            if i - sep(iter, 1) > 1
                lastDY = Y(i, :) - Y(i-1, :);
                lastX1 = zeros(1, 10);
                lastX1(6:10) = Y(i-1, 6:10) + Maaav(Y(i - 1, :), avp, U(i - 1, :))' * dt;
                lastX1(1:5) = Y(i-1, 1:5) + lastX1(6:10) * dt;
                lastX = zeros(1, 10);
                lastX(6:10) = Y(i, 6:10) + Maaav(Y(i, :), avp, U(i, :))' * dt;
                lastX(1:5) = Y(i, 1:5) + lastX(6:10) * dt;
                lastDX = lastX - lastX1;
                lastEf = (lastDY - lastDX) / dt;
                b = Phi \ lastEf';
                disc = real(Phi * Lambda * b)';
            end
            out(i, :) = out(i, :) + disc * dt;
        else
            disc = zeros(1, 10);
            if i - sep(iter, 1) > 1
                lastDY = Y(i, :) - Y(i-1, :);
                lastX1 = zeros(1, 10);
                lastX1(6:10) = Y(i-1, 6:10) + Maaav(Y(i - 1, :), avp, U(i - 1, :))' * dt;
                lastX1(1:5) = Y(i-1, 1:5) + lastX1(6:10) * dt;
                lastX = zeros(1, 10);
                lastX(6:10) = Y(i, 6:10) + Maaav(Y(i, :), avp, U(i, :))' * dt;
                lastX(1:5) = Y(i, 1:5) + lastX(6:10) * dt;
                lastDX = lastX - lastX1;
                lastEf = (lastDY - lastDX) / dt;
                b = Phi \ lastEf';
                disc = real(Phi * Lambda * b);
            end
            out(i, :) = disc';
        end
    end
end

end