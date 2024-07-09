% 延迟DMD处理带输入的动态系统

% 参数设置
n = 2;  % 状态维度
m = 1;  % 输入维度
M = 10; % 数据长度
d = 3; % 延迟嵌入维度

% 生成状态和输入数据
A = [1.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];
x = zeros(n, M);
% 输入
u = ones(1, M);

x(:,1) = [1; 1];

for k = 1:M-1
    x(:,k+1) = A*x(:,k) + B*u(k);
end

% 构建延迟嵌入矩阵
Z = [];
Z_prime = [];
for k = 1:M-d
    Z = [Z, reshape(x(:,k:k+d-1), [], 1)];
    Z_prime = [Z_prime, reshape(x(:,k+1:k+d), [], 1)];
end

% 构建扩展的输入矩阵
U = [];
for k = 1:M-d
    U = [U, reshape(u(k:k+d-1), [], 1)];
end

W = [Z; U];

% 求解动态矩阵 K
K = Z_prime * pinv(W);

% 分离矩阵 A_d 和 B_d
A_d = K(1:n*d, 1:n*d);
B_d = K(1:n*d, n*d+1:end);

% 模态分析
[Phi, Lambda] = eig(A_d);

% 初始条件
Z0 = Z(:,1);
b = Phi \ Z0;

% 重建状态
Z_reconstructed = zeros(size(Z));
for k = 0:M-d-1
    Z_reconstructed(:,k+1) = Phi * (Lambda^k) * b;
end

% 绘制结果
figure;
subplot(2,1,1);
plot(1:M, x(1,:), 'b', 'DisplayName', '真实数据');
hold on;
plot(d+1:M, Z_reconstructed(1,:), 'r--', 'DisplayName', '重建数据');
xlabel('时间步');
ylabel('状态x_1');
legend;

subplot(2,1,2);
plot(1:M, x(2,:), 'b', 'DisplayName', '真实数据');
hold on;
plot(d+1:M, Z_reconstructed(2,:), 'r--', 'DisplayName', '重建数据');
xlabel('时间步');
ylabel('状态x_2');
legend;

sgtitle('带输入的动态系统延迟DMD重建');

% 延迟DMD方法如何处理带输入的动态系统
