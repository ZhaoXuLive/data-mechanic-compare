function [Xi, Xintercept] = sparsifyDynamics_lars(Theta, dXdt)

Xi = zeros(size(Theta,2), 5);
Xintercept = zeros(5, 1);
for i = 1:5
    
    % 使用 lasso 函数进行最小角回归
    [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'CV', 10, 'Lambda', 'auto', 'DFmax', size(Theta, 2));
    
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