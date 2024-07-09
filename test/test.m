% 扩展DMD示例代码

% 生成示例数据
t = 0:1:9; % 时间向量，10个时间点
n = length(t); % 数据长度

% 系统矩阵和输入矩阵
A = [1.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];

% 生成输入信号
u = sin(0.1 * t); % 输入信号，可以根据需要调整

% 初始化状态向量
x0 = [1; 1]; % 初始状态
X = zeros(2, n); % 状态矩阵
X(:, 1) = x0;

% 生成状态矩阵
for i = 2:n
    X(:, i) = A * X(:, i-1) + B * u(i-1);
end

% 构建扩展状态矩阵，包括输入信号
X_ext = [X(:, 1:end-1); u(1:end-1)];
X_prime_ext = X(:, 2:end);

% DMD计算
[U, S, V] = svd(X_ext, 'econ');
Atilde = U' * X_prime_ext * V * diag(1./diag(S));
[W, D] = eig(Atilde);
Phi = X_prime_ext * V * diag(1./diag(S)) * W;

% 动态模式
lambda = diag(D);
omega = log(lambda); % 频率

% 初始状态
b = Phi \ X_prime_ext(:, 1); % 初始状态

% 时间演化
time_dynamics = zeros(length(b), n - 1);
for i = 1:n-1
    time_dynamics(:,i) = b .* exp(omega * t(i));
end

% 重构信号
x_reconstructed = Phi * time_dynamics;

% 绘图
figure;
plot(t, X(1, :), 'b', 'DisplayName', '原始信号');
hold on;
plot(t(2:end), real(x_reconstructed(1, :)), 'r--', 'DisplayName', '重构信号');
legend;
title('信号重构');
xlabel('时间');
ylabel('状态变量');

% 验证集上的预测
t_pred = 10:1:14; % 新的时间向量
u_pred = sin(0.1 * t_pred); % 验证集输入信号
X_pred = zeros(2, length(t_pred));
X_pred(:, 1) = X(:, end);

for i = 2:length(t_pred)
    X_pred(:, i) = A * X_pred(:, i-1) + B * u_pred(i-1);
end

X_pred_ext = [X_pred(:, 1:end-1); u_pred(1:end-1)];

% 重构验证集信号
time_dynamics_pred = zeros(length(b), length(t_pred)-1);
for i = 1:length(t_pred)-1
    time_dynamics_pred(:,i) = b .* exp(omega * (t_pred(i) - t(end)));
end

x_reconstructed_pred = Phi * time_dynamics_pred;

% 绘图
figure;
plot(t_pred, X_pred(1, :), 'b', 'DisplayName', '原始信号');
hold on;
plot(t_pred(2:end), real(x_reconstructed_pred(1, :)), 'r--', 'DisplayName', '重构信号');
legend;
title('验证集信号重构');
xlabel('时间');
ylabel('状态变量');
