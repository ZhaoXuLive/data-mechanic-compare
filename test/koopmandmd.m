
% 定义系统参数
A = [0.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];
x0 = [1; 1];
u = 1;
m = 10;  % 时间点数量

% 生成状态矩阵 X 和 X_prime
X = zeros(2, m);
X_prime = zeros(2, m);
U = u * ones(1, m);
X(:, 1) = x0;

for k = 1:m-1
    X_prime(:, k) = A * X(:, k) + B * U(:, k);
    X(:, k+1) = X_prime(:, k);
end
X_prime(:, m) = A * X(:, m) + B * U(:, m);

% 使用多项式基函数进行扩展
phi = @(x) [x(1); x(2); x(1)^2; x(2)^2; x(1)*x(2)];
Psi_X = zeros(5, m);
Psi_X_prime = zeros(5, m);

for k = 1:m
    Psi_X(:, k) = phi(X(:, k));
    Psi_X_prime(:, k) = phi(X_prime(:, k));
end

%%
% % Koopman DMD 分析
% [Phi, Lambda] = koopmanDMD(Psi_X, Psi_X_prime);

% 矩阵大小
[q, m] = size(Psi_X);

% 计算 K 矩阵
K = Psi_X_prime * pinv(Psi_X);

% 特征分解 K
[W, Lambda] = eig(K);

% 计算 Koopman DMD 模态
Phi = X_prime * W * inv(Lambda);
      
%%
% 重建训练数据
X_reconstructed = reconstructData(Phi, Lambda, Psi_X(:,1), m);

% 打印结果
disp('Koopman DMD Modes (Phi):');
disp(Phi);
disp('Eigenvalues (Lambda):');
disp(diag(Lambda));
disp('Reconstructed Training Data:');
disp(X_reconstructed);

% 测试集数据验证
n_test = 5;  % 测试集时间点数量
X_test = zeros(2, n_test);
X_test(:, 1) = X(:, end);  % 测试集初始状态为训练集末状态

for k = 1:n_test-1
    X_test(:, k+1) = A * X_test(:, k) + B * u;
end

% 扩展测试数据
Psi_X_test = zeros(5, n_test);
for k = 1:n_test
    Psi_X_test(:, k) = phi(X_test(:, k));
end

% 预测测试集数据
X_test_predicted = reconstructData(Phi, Lambda, Psi_X_test(:,1), n_test);

% 计算误差
error = norm(X_test - X_test_predicted(1:2, :), 'fro');
disp('Test Data:');
disp(X_test);
disp('Predicted Test Data:');
disp(X_test_predicted(1:2, :));
disp(['Reconstruction Error: ', num2str(error)]);

function X_reconstructed = reconstructData(Phi, Lambda, Psi_x0, num_steps)
    % 使用 DMD 模态和特征值重建数据
    % 输入:
    % Phi         - DMD 模态矩阵
    % Lambda      - DMD 特征值矩阵
    % Psi_x0      - 初始状态在扩展空间的表示
    % num_steps   - 时间步数量
    % 输出:
    % X_reconstructed - 重建的数据矩阵

    q = size(Phi, 1);  % 扩展空间维度
    X_reconstructed = zeros(q, num_steps);
    X_reconstructed(:, 1) = Psi_x0;

    for k = 2:num_steps
        X_reconstructed(:, k) = Phi * (Lambda^(k-1)) * (Phi \ Psi_x0);
    end
end
