function Phi_optimized = optimizeModes(Phi, Lambda, b, X, m, itern, learnv)
    % 优化DMD模态的函数，采用最小二乘法优化
    % 这里我们只是展示一个简单的优化示例
    % 实际应用中可以使用更复杂的优化方法

    % 初始化优化后的Phi
    Phi_optimized = Phi;

    % 优化步骤
    for i = 1:itern % 优化迭代次数
        % 计算重建误差
        X_reconstructed = zeros(size(X));
        for t = 1:m
            X_reconstructed(:,t) = Phi_optimized * (Lambda^(t-1)) * b;
        end
        error = X - X_reconstructed;
        
        avgError = zeros(10, 1);
        for j = 1:10
            avgError(j) = mean(error(j, :));
        end
        
        % 更新Phi_optimized
        Phi_optimized = Phi_optimized + learnv * avgError * b'; 
    end
end