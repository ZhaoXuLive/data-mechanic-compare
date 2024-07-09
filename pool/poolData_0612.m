function yout = poolData_new(yin, u, polyorder)
% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

yn = size(yin,1);
% yout = zeros(n,1+nVars+(nVars*(nVars+1)/2)+(nVars*(nVars+1)*(nVars+2)/(2*3))+11);

%% initial
% 自主设计候选函数 474

% 三角函数
tri = zeros(yn, 12);
tri(:, 1) = sin(yin(:, 1));
tri(:, 2) = sin(yin(:, 2));
tri(:, 3) = sin(yin(:, 3));
tri(:, 4) = cos(yin(:, 1));
tri(:, 5) = cos(yin(:, 2));
tri(:, 6) = cos(yin(:, 3));
tri(:, 7) = sin(yin(:, 2) - yin(:, 1));
tri(:, 8) = sin(yin(:, 3) - yin(:, 2));
tri(:, 9) = sin(yin(:, 3) - yin(:, 1));
tri(:, 10) = cos(yin(:, 2) - yin(:, 1));
tri(:, 11) = cos(yin(:, 3) - yin(:, 2));
tri(:, 12) = cos(yin(:, 3) - yin(:, 1));

% 速度
vl = zeros(yn, 2);
vl(:, 1) = yin(:, 4);
vl(:, 2) = yin(:, 5);

% 转向角
del = zeros(yn, 8);
del(:, :) = u(:, :);

% 角速度
wd = zeros(yn, 3);
wd(:, :) = yin(:, 6:8);

%% mix
allmix = [tri vl del wd];

ind = 1;
% poly order 0
yout(:,ind) = ones(yn,1);
ind = ind+1;

% poly order 1
for i=1:length(allmix(1, :))
    yout(:,ind) = allmix(:,i);
    ind = ind+1;
end

% poly order 2
for i=1:length(allmix(1, :))
    for j=1:length(allmix(1, :))
        yout(:,ind) = allmix(:,i).*allmix(:,j);
        ind = ind+1;
    end
end

% poly order 3
if polyorder >= 3
    for i=1:length(tri(1, :))
        for j=1:length(allmix(1, :))
            for k=1:length(allmix(1, :))
                yout(:,ind) = tri(:,i).*allmix(:,j).*allmix(:,k);
                ind = ind+1;
            end
        end
    end
    for i=1:length(wd(1, :))
        for j=1:length(wd(1, :))
            for k=1:length(del(1, :))
                yout(:,ind) = wd(:,i).*wd(:,j).*del(:,k);
                ind = ind+1;
            end
        end
    end
    for i=1:length(wd(1, :))
        for j=1:length(wd(1, :))
            for k=1:length(vl(1, :))
                yout(:,ind) = wd(:,i).*wd(:,j).*vl(:,k);
                ind = ind+1;
            end
        end
    end
    for i=1:length(wd(1, :))
        for j=1:length(del(1, :))
            for k=1:length(vl(1, :))
                yout(:,ind) = wd(:,i).*del(:,j).*vl(:,k);
                ind = ind+1;
            end
        end
    end
    for i=1:length(wd(1, :))
        for j=1:length(vl(1, :))
            for k=1:length(vl(1, :))
                yout(:,ind) = wd(:,i).*vl(:,j).*vl(:,k);
                ind = ind+1;
            end
        end
    end
end

if polyorder >= 4
    % poly order 4
    for i=1:length(tri(1, :))
        for j=1:length(tri(1, :))
            for k=1:length(del(1, :))
                for l=1:length(wd(1, :))
                    yout(:,ind) = tri(:,i).*tri(:,j).*del(:,k).*wd(:,l);
                    ind = ind+1;
                end
            end
        end
    end
    for i=1:length(tri(1, :))
        for j=1:length(tri(1, :))
            for k=1:length(del(1, :))
                for l=1:length(vl(1, :))
                    yout(:,ind) = tri(:,i).*tri(:,j).*del(:,k).*vl(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

end