function yout = poolData_withU(yin, u, polyorder)
% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

yn = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

%% initial
% 把y和u综合放下一起
% 引入cos和sin 3276? over and over
% fxlen = 6+6+5+8;
% uy = zeros(yn, fxlen);
% uy(:, 1) = sin(yin(:, 1));
% uy(:, 2) = sin(yin(:, 2));
% uy(:, 3) = sin(yin(:, 3));
% uy(:, 4) = cos(yin(:, 1));
% uy(:, 5) = cos(yin(:, 2));
% uy(:, 6) = cos(yin(:, 3));
% uy(:, 7) = sin(yin(:, 2) - yin(:, 1));
% uy(:, 8) = sin(yin(:, 3) - yin(:, 2));
% uy(:, 9) = sin(yin(:, 3) - yin(:, 1));
% uy(:, 10) = cos(yin(:, 2) - yin(:, 1));
% uy(:, 11) = cos(yin(:, 3) - yin(:, 2));
% uy(:, 12) = cos(yin(:, 3) - yin(:, 1));
% uy(:, 12+1:12+5) = yin(:, 4:8);
% uy(:, 12+5+1:end) = u;

% 1540
fxlen = 6+5+8;
uy = zeros(yn, fxlen);
uy(:, 1) = sin(yin(:, 1));
uy(:, 2) = sin(yin(:, 2));
uy(:, 3) = sin(yin(:, 3));
uy(:, 4) = cos(yin(:, 1));
uy(:, 5) = cos(yin(:, 2));
uy(:, 6) = cos(yin(:, 3));
uy(:, 6+1:6+5) = yin(:, 4:8);
uy(:, 6+5+1:end) = u;


%% add
ind = 1;
% poly order 0
yout(:,ind) = ones(yn,1);
ind = ind+1;

% poly order 1
for i=1:fxlen
    yout(:,ind) = uy(:,i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:fxlen
        for j=i:fxlen
            yout(:,ind) = uy(:,i).*uy(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:fxlen
        for j=i:fxlen
            for k=j:fxlen
                yout(:,ind) = uy(:,i).*uy(:,j).*uy(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:fxlen
        for j=i:fxlen
            for k=j:fxlen
                for l=k:fxlen
                    yout(:,ind) = uy(:,i).*uy(:,j).*uy(:,k).*uy(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:fxlen
        for j=i:fxlen
            for k=j:fxlen
                for l=k:fxlen
                    for m=l:fxlen
                        yout(:,ind) = uy(:,i).*uy(:,j).*uy(:,k).*uy(:,l).*uy(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

end