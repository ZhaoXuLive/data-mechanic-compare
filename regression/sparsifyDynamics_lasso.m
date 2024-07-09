function Xi = sparsifyDynamics_lasso(Theta, dXdt)

% % 数据标准化（可选）
% Theta_mean = mean(Theta);
% Theta_std = std(Theta);
% Theta_std(Theta_std == 0) = 1; % 避免除以零
% Theta = (Theta - Theta_mean) ./ Theta_std;
% dXdt_mean = mean(dXdt);
% dXdt = dXdt - dXdt_mean;

Xi = zeros(size(Theta,2), 5);
for i = 1:5
    [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'CV', 10);
    % 选择一个合适的lambda值，例如使用交叉验证的最优值
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    MaxLambda = FitInfo.Lambda1SE
    Xi(:, i) = B(:, idxLambdaMinMSE);
end

end