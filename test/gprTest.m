%% �Իع�GPR

% % ����һЩ��������
% X = (1:10)';
% y = sin(X) + 0.1*randn(10,1);
% 
% % ʹ��Ĭ��������ϸ�˹���̻ع�ģ��
% gprMdl = fitrgp(X, y);
% 
% % ʹ����ϵ�ģ�ͽ���Ԥ��
% Xnew = (1:0.1:10)';
% [ypred, ysd, yint] = predict(gprMdl, Xnew);
% 
% % ���ƽ��
% figure;
% hold on;
% plot(X, y, 'ko', 'MarkerFaceColor', 'r'); % ������
% plot(Xnew, ypred, 'b-'); % Ԥ��ֵ
% plot(Xnew, yint, 'm--'); % Ԥ������
% hold off;
% legend('Data', 'Prediction', 'Prediction Interval');
% xlabel('X');
% ylabel('Y');
% title('Gaussian Process Regression');

%% �������GPR

% ϵͳά��
n = 1; % ״̬ά��
m = 1; % ����ά��

% ��ʼ״̬
x0 = 0.5;

% �����ź�
N = 100; % ʱ�䲽��
U = randn(m, N); % �������

% ����״̬����
X = zeros(n, N+1);
X(:,1) = x0;
Y = zeros(n, N);

for k = 1:N
    X(:,k+1) = sin(X(:,k)) + U(:,k) + 0.1*randn;
    Y(:,k) = X(:,k+1); % Ŀ������һ��״̬
end

% ״̬���������ݾ���
X_data = X(:,1:N);
U_data = U;

% ��״̬���������ݺϲ�Ϊѵ������
trainData = [X_data; U_data]';

% GPRģ��
gprMdl = fitrgp(trainData, Y');

% ��������
x_test = linspace(-2, 2, 100)';
u_test = zeros(size(x_test));
testData = [x_test u_test];

% Ԥ��
[y_pred, y_sd] = predict(gprMdl, testData);

% ���ƽ��
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
