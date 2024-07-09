function Xi = sparsifyDynamics_lasso(Theta, dXdt)

% % ���ݱ�׼������ѡ��
% Theta_mean = mean(Theta);
% Theta_std = std(Theta);
% Theta_std(Theta_std == 0) = 1; % ���������
% Theta = (Theta - Theta_mean) ./ Theta_std;
% dXdt_mean = mean(dXdt);
% dXdt = dXdt - dXdt_mean;

Xi = zeros(size(Theta,2), 5);
for i = 1:5
    [B, FitInfo] = lasso(Theta, dXdt(:, i+5), 'CV', 10);
    % ѡ��һ�����ʵ�lambdaֵ������ʹ�ý�����֤������ֵ
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    MaxLambda = FitInfo.Lambda1SE
    Xi(:, i) = B(:, idxLambdaMinMSE);
end

end