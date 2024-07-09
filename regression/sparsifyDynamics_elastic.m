function [Xi, Xintercept] = sparsifyDynamics_elastic(Theta, dXdt)

Xi = zeros(size(Theta,2), 5);
Xintercept = zeros(5, 1);
for i = 1:5
    
    % 设置 alpha 值
    alpha = 0.5;  % Elastic Net 回归中的 alpha 参数
    
    % 使用 lasso 函数进行交叉验证以选择最佳的 lambda
    [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'Alpha', alpha, 'CV', 10);  % 10 折交叉验证
    
    % 找到具有最小均方误差的最佳 lambda
    minMSEIndex = FitInfo.IndexMinMSE;
%     bestLambda = FitInfo.Lambda(minMSEIndex);
    
    % 获取模型系数和截距
    Xi(:, i) = B(:, minMSEIndex);
    Xintercept(i) = FitInfo.Intercept(minMSEIndex);
    
%     % 绘制 MSE 和 Lambda 的关系图
%     lassoPlot(B, FitInfo, 'PlotType', 'CV');

end

end