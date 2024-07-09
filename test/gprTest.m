%% 自回归GPR

% % 生成一些样本数据
% X = (1:10)';
% y = sin(X) + 0.1*randn(10,1);
% 
% % 使用默认设置拟合高斯过程回归模型
% gprMdl = fitrgp(X, y);
% 
% % 使用拟合的模型进行预测
% Xnew = (1:0.1:10)';
% [ypred, ysd, yint] = predict(gprMdl, Xnew);
% 
% % 绘制结果
% figure;
% hold on;
% plot(X, y, 'ko', 'MarkerFaceColor', 'r'); % 样本点
% plot(Xnew, ypred, 'b-'); % 预测值
% plot(Xnew, yint, 'm--'); % 预测区间
% hold off;
% legend('Data', 'Prediction', 'Prediction Interval');
% xlabel('X');
% ylabel('Y');
% title('Gaussian Process Regression');

%% 带输入的GPR

% 系统维度
n = 1; % 状态维度
m = 1; % 输入维度

% 初始状态
x0 = 0.5;

% 输入信号
N = 100; % 时间步数
U = randn(m, N); % 随机输入

% 生成状态数据
X = zeros(n, N+1);
X(:,1) = x0;
Y = zeros(n, N);

for k = 1:N
    X(:,k+1) = sin(X(:,k)) + U(:,k) + 0.1*randn;
    Y(:,k) = X(:,k+1); % 目标是下一个状态
end

% 状态和输入数据矩阵
X_data = X(:,1:N);
U_data = U;

% 将状态和输入数据合并为训练数据
trainData = [X_data; U_data]';

% GPR模型
gprMdl = fitrgp(trainData, Y');

% 测试数据
x_test = linspace(-2, 2, 100)';
u_test = zeros(size(x_test));
testData = [x_test u_test];

% 预测
[y_pred, y_sd] = predict(gprMdl, testData);

% 绘制结果
figure;
hold on;
plot(X(1,1:end-1), Y, 'r.', 'MarkerSize', 15);
plot(x_test, y_pred, 'b-', 'LineWidth', 2);
fill([x_test; flipud(x_test)], [y_pred-y_sd; flipud(y_pred+y_sd)], [7 7 7]/8, 'EdgeColor', 'none');
xlabel('State x_k');
ylabel('State x_{k+1}');
legend('Training data', 'GPR prediction', '95% prediction interval');
title('Gaussian Process Regression for Dynamic System');
hold off;
