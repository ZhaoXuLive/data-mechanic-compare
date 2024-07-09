function [Xi, Xintercept] = sparsifyDynamics_ridge(Theta, dXdt)

Xi = zeros(size(Theta,2), 5);
Xintercept = zeros(5, 1);
lassoLambda = [9.6128e-08 2.5671e-07 9.7338e-09 8.6769e-09 1.1747e-08];
for i = 1:5
    
%     % ���� alpha ֵ
%     alpha = 0;  % 0 ridge 1 lasso
%     
%     % ʹ�� lasso �������н�����֤��ѡ����ѵ� lambda
%     [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'Alpha', alpha, 'CV', 10);  % 10 �۽�����֤
%     
%     % �ҵ�������С����������� lambda
%     minMSEIndex = FitInfo.IndexMinMSE;
% %     bestLambda = FitInfo.Lambda(minMSEIndex);
%     
%     % ��ȡģ��ϵ���ͽؾ�
%     Xi(:, i) = B(:, minMSEIndex);
%     Xintercept(i) = FitInfo.Intercept(minMSEIndex);
%     
% %     % ���� MSE �� Lambda �Ĺ�ϵͼ
% %     lassoPlot(B, FitInfo, 'PlotType', 'CV');

%     lambdas = logspace(-4, 1, 50);  % �� 10^-4 �� 10^1 ֮���20��ֵ
% 
%     % ʹ�ý�����֤ѡ����Ѳ���
%     Mdl = fitrlinear(Theta, dXdt(:, i+5), 'Learner', 'leastsquares', 'Regularization', 'ridge', ...
%                      'Lambda', lambdas, 'KFold', 5);
% 
%     % �ҵ���õ�Lambdaֵ
%     [~, idx] = min(kfoldLoss(Mdl));
%     bestLambda = lambdas(idx);
%     fprintf('���Lambda: %f\n', bestLambda);

    bestLambda = lassoLambda(i);
    % ʹ�����Lambda����ѵ��ģ��
    ridgeModelBest = fitrlinear(Theta, dXdt(:, i+5), 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);

    % ���ڴ���ϵ���ͽؾ࣬���ܼ򵥽�ϵ���洢��������Ҫ��һ����¼�ؾ�
    Xi(:, i) = ridgeModelBest.Beta;
    Xintercept(i) = ridgeModelBest.Bias;

end

end