function yout = poolDataLIST_withU(uyin, ahat, nVars, un, polyorder)
% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

% poolDataLIST_withU(signy, Xi, n, un, polyorder);

ind = 1;
% poly order 0
yout{ind,1} = ['1'];
ind = ind+1;

% poly order 1
for i=1:nVars+un
    yout(ind,1) = uyin(i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars+un
        for j=i:nVars+un
            yout{ind,1} = [uyin{i},uyin{j}];
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars+un
        for j=i:nVars+un
            for k=j:nVars+un
                yout{ind,1} = [uyin{i},uyin{j},uyin{k}];
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars+un
        for j=i:nVars+un
            for k=j:nVars+un
                for l=k:nVars+un
                    yout{ind,1} = [uyin{i},uyin{j},uyin{k},uyin{l}];
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars+un
        for j=i:nVars+un
            for k=j:nVars+un
                for l=k:nVars+un
                    for m=l:nVars+un
                        yout{ind,1} = [uyin{i},uyin{j},uyin{k},uyin{l},uyin{m}];
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

yin = uyin(1:10);

output = yout;
newout(1) = {''};
title = {'ax', 'ay', 'psi1dd', 'psi2dd','psi3dd'};
for k=1:5
    newout{1,1+k} = [title{k}];
end
% newout = {'','xdot','ydot','udot'};
for k=1:size(ahat,1)
    newout(k+1,1) = output(k);
    for j=1:5
        newout{k+1,1+j} = ahat(k,j);
    end
end
newout