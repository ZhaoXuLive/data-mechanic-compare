function [Xi, Xintercept] = sparsifyDynamics_lars(Theta, dXdt)

Xi = zeros(size(Theta,2), 5);
Xintercept = zeros(5, 1);
for i = 1:5
    
    % ʹ�� lasso ����������С�ǻع�
    [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'CV', 10, 'Lambda', 'auto', 'DFmax', size(Theta, 2));
    
    % �ҵ�������С����������� lambda
    minMSEIndex = FitInfo.IndexMinMSE;
%     bestLambda = FitInfo.Lambda(minMSEIndex);
    
    % ��ȡģ��ϵ���ͽؾ�
    Xi(:, i) = B(:, minMSEIndex);
    Xintercept(i) = FitInfo.Intercept(minMSEIndex);
    
%     % ���� MSE �� Lambda �Ĺ�ϵͼ
%     lassoPlot(B, FitInfo, 'PlotType', 'CV');

end

end