
% 系统矩阵
A = [1.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];

% 初始状态
x0 = [1; 1];

% 输入
u = 1;

% 时间点数量
m = 10;

% 生成状态矩阵 X 和 X_prime
X = zeros(2, m);
X_prime = zeros(2, m);
U = u * ones(1, m);

% 初始状态赋值
X(:, 1) = x0;

% 生成数据
for k = 1:m-1
    X_prime(:, k) = A * X(:, k) + B * U(:, k);
    X(:, k+1) = X_prime(:, k);
end
X_prime(:, m) = A * X(:, m) + B * U(:, m);

% 调用优化DMD函数
optimizedDMDwithInputs(X, X_prime, U);

function optimizedDMDwithInputs(X, X_prime, U)
    % 输入:
    % X       - 状态矩阵，大小为 n x m
    % X_prime - 下一个时间步的状态矩阵，大小为 n x m
    % U       - 输入矩阵，大小为 p x m

    % 矩阵大小
    [n, m] = size(X);
    p = size(U, 1);

    % 构建扩展矩阵 Omega
    Omega = [X; U];

    % 计算 Omega 的伪逆
    Omega_pseudo_inverse = pinv(Omega);

    % 计算 K 矩阵
    K = X_prime * Omega_pseudo_inverse;

    % 从 K 中分解出 A 和 B
    A = K(:, 1:n);
    B = K(:, n+1:end);

    % 奇异值分解 X
    [U_X, Sigma_X, V_X] = svd(X, 'econ');

    % 计算降维后的 A_tilde
    A_tilde = U_X' * A * U_X;

    % 特征分解 A_tilde
    [W, Lambda] = eig(A_tilde);

    % 计算 DMD 模态
    Phi = X_prime * V_X / Sigma_X * W;

    % 初始条件的DMD坐标
    b = Phi \ X(:,1);

    % 优化DMD步（这里以最小二乘法优化为例）
    Phi_optimized = optimizeModes(Phi, Lambda, b, X, m);

    % 重建训练集数据
    X_dmd = zeros(n, m);
    for i = 1:m
        X_dmd(:,i) = Phi_optimized * (Lambda^(i-1)) * b;
    end

    % 输出结果
    disp('A matrix:');
    disp(A);
    disp('B matrix:');
    disp(B);
    
    t = 1:1:10;
    figure,
    plot(t, X(1, :));
    hold on;
    plot(t, X_dmd(1, :), 'r');
    figure,
    plot(t, X(2, :));
    hold on;
    plot(t, X_dmd(2, :), 'r');

end

function Phi_optimized = optimizeModes(Phi, Lambda, b, X, m)
    % 优化DMD模态的函数，采用最小二乘法优化
    % 这里我们只是展示一个简单的优化示例
    % 实际应用中可以使用更复杂的优化方法

    % 初始化优化后的Phi
    Phi_optimized = Phi;

    % 优化步骤
    for i = 1:10 % 优化迭代次数
        % 计算重建误差
        X_reconstructed = zeros(size(X));
        for t = 1:m
            X_reconstructed(:,t) = Phi_optimized * (Lambda^(t-1)) * b;
        end
        error = X - X_reconstructed;
        
        avgError = zeros(2, 1);
        avgError(1) = mean(error(1, :));
        avgError(2) = mean(error(2, :));
        
        % 更新Phi_optimized
        Phi_optimized = Phi_optimized + 0.01 * avgError * b'; % 学习率为0.01
    end
end
