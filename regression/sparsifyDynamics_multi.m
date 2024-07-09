function [Xi, Xintercept, Sigvars] = sparsifyDynamics_multi(Theta, dXdt)

%% 基础操作

% 获取可用的 CPU 核心数量
numCores = feature('numcores');

% 启动并行池
parpool(numCores);

% 开始计时
tStart = tic;

%% 核心内容

Xi = zeros(size(Theta,2), 5);
Sigvars = zeros(size(Theta,2), 5);
X_train = Theta;

for i = 1:5
    y_train = dXdt(:, i+5);
    
    disp(i);
    
    % 使用parfor分块进行逐步回归
    % 定义分块数
    numBlocks = numCores;
    trainIdx = size(Theta,1);
    blockSize = floor(trainIdx / numBlocks);

    % 初始化模型存储单元和显著变量索引
    mdl_blocks = cell(numBlocks, 1);
    sig_vars_blocks = cell(numBlocks, 1);

    % 使用 parfor 循环进行分块回归
    parfor j = 1:numBlocks
        startIdx = (j-1) * blockSize + 1;
        if j == numBlocks
            endIdx = trainIdx;
        else
            endIdx = j * blockSize;
        end
        X_block = X_train(startIdx:endIdx, :);
        y_block = y_train(startIdx:endIdx);
        mdl_blocks{j} = stepwiselm(X_block, y_block, 'linear', 'Upper', 'linear', 'Lower', 'constant', ...
                       'PEnter', 0.05, 'PRemove', 0.10, 'Verbose', 0);
        % 获取显著变量索引
        pValues = mdl_blocks{j}.anova.pValue;
        sig_vars_blocks{j} = find(pValues < 0.05); % 选择显著性水平为0.05
    end

    % 合并显著变量索引
    sig_vars = unique(cell2mat(sig_vars_blocks));
    
    % 记录变量索引，用于测试集
    Sigvars(1, i) = length(sig_vars);
    Sigvars(2:length(sig_vars)+1, i) = sig_vars;

    % 使用显著变量选择特征
    X_train_sig = X_train(:, sig_vars);
    
    % 使用LASSO进行回归
    [B, FitInfo] = lasso(X_train_sig, y_train, 'CV', 10);

    % 获取最佳LASSO模型的系数
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    B_optimal = B(:, idxLambdaMinMSE);
    lasso_intercept = FitInfo.Intercept(idxLambdaMinMSE);
    
    Xi(1, i) = length(B_optimal);
    Xi(2:length(B_optimal)+1, i) = B_optimal;
    Xintercept(i) = lasso_intercept;
end

tEnd = toc(tStart)
delete(gcp);

end