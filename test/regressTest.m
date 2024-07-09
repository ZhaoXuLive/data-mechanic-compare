% 1. 生成数据
n = 1000; % 样本数量
p = 100; % 特征数量
X = randn(n, p); % 随机生成特征矩阵
true_beta = [3; -2; 0; 0; 1; zeros(p-5, 1)]; % 真实回归系数，部分为0
y = X * true_beta + randn(n, 1); % 添加噪声的响应变量

% 2. 数据标准化（逐步回归前处理）
X_mean = mean(X);
X_std = std(X);
X_std(X_std == 0) = 1; % 避免除以零
X = (X - X_mean) ./ X_std;
y_mean = mean(y);
y = y - y_mean;

%% lasso

% % 3. 使用Lasso回归进行变量选择
% [B, FitInfo] = lasso(X, y, 'CV', 10);
% lassoLambda = FitInfo.Lambda1SE; % 使用1-SE法则选择Lambda
% lassoBeta = B(:, FitInfo.Index1SE);
% 
% % 打印结果
% selectedLassoVars = find(lassoBeta ~= 0);
% fprintf('Lasso选择的变量索引: %s\n', mat2str(selectedLassoVars));
% fprintf('Lasso选择的系数: \n');
% disp(lassoBeta(lassoBeta ~= 0));
% 
% % 使用 LassoBeta 预测响应变量
% y_pred = X * lassoBeta;
% 
% % 计算均方误差 (MSE)
% MSE = mean((y - y_pred).^2);
% 
% % 输出结果
% fprintf('Lasso 回归的均方误差 (MSE): %.4f\n', MSE);

%% lasso + stepwise

% 获取可用的 CPU 核心数量
numCores = feature('numcores');

% 启动并行池
parpool(numCores);

tStart = tic;

% 定义分块数（不要超过可用核心数）
numBlocks = numCores;
blockSize = floor(n / numBlocks);

% 初始化模型存储单元
mdl = cell(numBlocks, 1);

% 进行PCA降维
[coeff, score, ~, ~, explained] = pca(X);

% 选择解释90%方差的主成分
numComponents = find(cumsum(explained) >= 90, 1);
X_reduced = score(:, 1:numComponents);

% 使用 parfor 循环进行分块回归
parfor i = 1:numBlocks
    startIdx = (i-1) * blockSize + 1;
    if i == numBlocks
        endIdx = n; % 确保最后一个块包含所有剩余数据
    else
        endIdx = i * blockSize;
    end
    X_block = X_reduced(startIdx:endIdx, :);
    y_block = y(startIdx:endIdx);
    mdl{i} = stepwiselm(X_block, y_block, 'linear', 'Upper', 'linear', 'Lower', 'constant', ...
                   'PEnter', 0.05, 'PRemove', 0.10, 'Verbose', 0);
end

% % 3. 使用逐步回归进行初步变量选择，确保只保留一阶自变量
% mdl = stepwiselm(X, y, 'linear', 'Upper', 'linear', 'Lower', 'constant', ...
%                    'PEnter', 0.05, 'PRemove', 0.10, 'Verbose', 0);

tEnd = toc(tStart)

% 关闭并行池
delete(gcp);
                      
% 4. 选择逐步回归中的显著变量
coeffs = mdl.Coefficients;
p_values = coeffs.pValue(2:end); % 排除截距项的p值
selectedVars = find(p_values < 0.10); % 选择显著变量

% 检查是否有显著变量
if ~isempty(selectedVars)
    % 5. 使用Lasso回归在逐步回归选择的变量基础上进行进一步优化
    X_selected = X(:, selectedVars);
    [B, FitInfo] = lasso(X_selected, y, 'CV', 10);
    lassoLambda = FitInfo.Lambda1SE; % 使用1-SE法则选择Lambda
    lassoBeta = B(:, FitInfo.Index1SE);

    % 打印结果
    selectedLassoVars = selectedVars(lassoBeta ~= 0);
    fprintf('Lasso选择的变量索引: %s\n', mat2str(selectedLassoVars));
    fprintf('Lasso选择的系数: \n');
    disp(lassoBeta(lassoBeta ~= 0));
    
    % 使用 LassoBeta 预测响应变量
    allBeta = zeros(p, 1);
    allBeta(1:size(lassoBeta(:, 1), 1), :)  = lassoBeta(:, :);
    y_pred = X * allBeta;

    % 计算均方误差 (MSE)
    MSE = mean((y - y_pred).^2);

    % 输出结果
    fprintf('StepWise + Lasso 回归的均方误差 (MSE): %.4f\n', MSE);
    
%     % 绘制Lasso路径图
%     figure;
%     lassoPlot(B, FitInfo, 'PlotType', 'Lambda', 'XScale', 'log');
%     xlabel('Lambda');
%     ylabel('回归系数');
%     title('Lasso路径图');
%     grid on;
else
    disp('逐步回归未选择任何变量。');
end

% tEnd = toc(tStart)
