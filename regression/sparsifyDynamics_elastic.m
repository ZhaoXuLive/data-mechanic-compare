function [Xi, Xintercept] = sparsifyDynamics_elastic(Theta, dXdt)

Xi = zeros(size(Theta,2), 5);
Xintercept = zeros(5, 1);
for i = 1:5
    
    % ���� alpha ֵ
    alpha = 0.5;  % Elastic Net �ع��е� alpha ����
    
    % ʹ�� lasso �������н�����֤��ѡ����ѵ� lambda
    [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'Alpha', alpha, 'CV', 10);  % 10 �۽�����֤
    
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