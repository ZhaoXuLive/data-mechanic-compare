%% basic parameter

% 获取可用的 CPU 核心数量
numCores = feature('numcores');

% 启动并行池
parpool(numCores);

% 生成示例数据
n = 10000; % 样本数量
p = 100; % 特征数量
X = randn(n, p);
y = X(:,1) - 2*X(:,2) + 0.5*randn(n, 1); % 假设的线性关系

% 将数据划分为训练集和测试集
trainRatio = 0.8;
trainIdx = 1:floor(trainRatio * n);
testIdx = (floor(trainRatio * n) + 1):n;

X_train = X(trainIdx, :);
y_train = y(trainIdx);
X_test = X(testIdx, :);
y_test = y(testIdx);

tStart = tic;

%% PCA降维

% 进行PCA降维
[coeff, score_train, ~, ~, explained] = pca(X_train);

% 选择解释90%方差的主成分
numComponents = find(cumsum(explained) >= 90, 1);
X_train_reduced = score_train(:, 1:numComponents);

% 将测试数据投影到主成分空间
X_test_reduced = X_test * coeff(:, 1:numComponents);

%% 使用parfor分块进行逐步回归

% 定义分块数
numBlocks = numCores;
blockSize = floor(length(trainIdx) / numBlocks);

% 初始化模型存储单元和显著变量索引
mdl_blocks = cell(numBlocks, 1);
sig_vars_blocks = cell(numBlocks, 1);

% 使用 parfor 循环进行分块回归
parfor i = 1:numBlocks
    startIdx = (i-1) * blockSize + 1;
    if i == numBlocks
        endIdx = length(trainIdx);
    else
        endIdx = i * blockSize;
    end
    X_block = X_train_reduced(startIdx:endIdx, :);
    y_block = y_train(startIdx:endIdx);
    mdl_blocks{i} = stepwiselm(X_block, y_block, 'linear', 'Upper', 'linear', 'Lower', 'constant', ...
                   'PEnter', 0.05, 'PRemove', 0.10, 'Verbose', 0);
    % 获取显著变量索引
    pValues = mdl_blocks{i}.anova.pValue;
    sig_vars_blocks{i} = find(pValues < 0.05); % 选择显著性水平为0.05
end

% 合并显著变量索引
sig_vars = unique(cell2mat(sig_vars_blocks));

% 使用显著变量选择特征
X_train_sig = X_train_reduced(:, sig_vars);
X_test_sig = X_test_reduced(:, sig_vars);

%% LASSO回归

% 使用LASSO进行回归
[B, FitInfo] = lasso(X_train_sig, y_train, 'CV', 10);

% 获取最佳LASSO模型的系数
idxLambdaMinMSE = FitInfo.IndexMinMSE;
B_optimal = B(:, idxLambdaMinMSE);
lasso_intercept = FitInfo.Intercept(idxLambdaMinMSE);

tEnd = toc(tStart)

% 对训练数据进行预测
y_train_pred_lasso = X_train_sig * B_optimal + lasso_intercept;

% 对测试数据进行预测
y_test_pred_lasso = X_test_sig * B_optimal + lasso_intercept;

% 计算LASSO回归的训练集和测试集的MSE
mse_train_lasso = mean((y_train - y_train_pred_lasso).^2);
mse_test_lasso = mean((y_test - y_test_pred_lasso).^2);

% 关闭并行池
delete(gcp);

%% 输出结果

disp(['LASSO回归的训练集MSE: ', num2str(mse_train_lasso)]);
disp(['LASSO回归的测试集MSE: ', num2str(mse_test_lasso)]);
