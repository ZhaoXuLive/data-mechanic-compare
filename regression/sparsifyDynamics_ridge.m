function [Xi, Xintercept] = sparsifyDynamics_ridge(Theta, dXdt)

Xi = zeros(size(Theta,2), 5);
Xintercept = zeros(5, 1);
lassoLambda = [9.6128e-08 2.5671e-07 9.7338e-09 8.6769e-09 1.1747e-08];
for i = 1:5
    
%     % 设置 alpha 值
%     alpha = 0;  % 0 ridge 1 lasso
%     
%     % 使用 lasso 函数进行交叉验证以选择最佳的 lambda
%     [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'Alpha', alpha, 'CV', 10);  % 10 折交叉验证
%     
%     % 找到具有最小均方误差的最佳 lambda
%     minMSEIndex = FitInfo.IndexMinMSE;
% %     bestLambda = FitInfo.Lambda(minMSEIndex);
%     
%     % 获取模型系数和截距
%     Xi(:, i) = B(:, minMSEIndex);
%     Xintercept(i) = FitInfo.Intercept(minMSEIndex);
%     
% %     % 绘制 MSE 和 Lambda 的关系图
% %     lassoPlot(B, FitInfo, 'PlotType', 'CV');

%     lambdas = logspace(-4, 1, 50);  % 从 10^-4 到 10^1 之间的20个值
% 
%     % 使用交叉验证选择最佳参数
%     Mdl = fitrlinear(Theta, dXdt(:, i+5), 'Learner', 'leastsquares', 'Regularization', 'ridge', ...
%                      'Lambda', lambdas, 'KFold', 5);
% 
%     % 找到最好的Lambda值
%     [~, idx] = min(kfoldLoss(Mdl));
%     bestLambda = lambdas(idx);
%     fprintf('最佳Lambda: %f\n', bestLambda);

    bestLambda = lassoLambda(i);
    % 使用最佳Lambda重新训练模型
    ridgeModelBest = fitrlinear(Theta, dXdt(:, i+5), 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);

    % 由于存在系数和截距，不能简单将系数存储起来，需要进一步记录截距
    Xi(:, i) = ridgeModelBest.Beta;
    Xintercept(i) = ridgeModelBest.Bias;

end

end