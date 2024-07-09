%% basic parameter

% ��ȡ���õ� CPU ��������
numCores = feature('numcores');

% �������г�
parpool(numCores);

% ����ʾ������
n = 10000; % ��������
p = 100; % ��������
X = randn(n, p);
y = X(:,1) - 2*X(:,2) + 0.5*randn(n, 1); % ��������Թ�ϵ

% �����ݻ���Ϊѵ�����Ͳ��Լ�
trainRatio = 0.8;
trainIdx = 1:floor(trainRatio * n);
testIdx = (floor(trainRatio * n) + 1):n;

X_train = X(trainIdx, :);
y_train = y(trainIdx);
X_test = X(testIdx, :);
y_test = y(testIdx);

tStart = tic;

%% PCA��ά

% ����PCA��ά
[coeff, score_train, ~, ~, explained] = pca(X_train);

% ѡ�����90%��������ɷ�
numComponents = find(cumsum(explained) >= 90, 1);
X_train_reduced = score_train(:, 1:numComponents);

% ����������ͶӰ�����ɷֿռ�
X_test_reduced = X_test * coeff(:, 1:numComponents);

%% ʹ��parfor�ֿ�����𲽻ع�

% ����ֿ���
numBlocks = numCores;
blockSize = floor(length(trainIdx) / numBlocks);

% ��ʼ��ģ�ʹ洢��Ԫ��������������
mdl_blocks = cell(numBlocks, 1);
sig_vars_blocks = cell(numBlocks, 1);

% ʹ�� parfor ѭ�����зֿ�ع�
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
    % ��ȡ������������
    pValues = mdl_blocks{i}.anova.pValue;
    sig_vars_blocks{i} = find(pValues < 0.05); % ѡ��������ˮƽΪ0.05
end

% �ϲ�������������
sig_vars = unique(cell2mat(sig_vars_blocks));

% ʹ����������ѡ������
X_train_sig = X_train_reduced(:, sig_vars);
X_test_sig = X_test_reduced(:, sig_vars);

%% LASSO�ع�

% ʹ��LASSO���лع�
[B, FitInfo] = lasso(X_train_sig, y_train, 'CV', 10);

% ��ȡ���LASSOģ�͵�ϵ��
idxLambdaMinMSE = FitInfo.IndexMinMSE;
B_optimal = B(:, idxLambdaMinMSE);
lasso_intercept = FitInfo.Intercept(idxLambdaMinMSE);

tEnd = toc(tStart)

% ��ѵ�����ݽ���Ԥ��
y_train_pred_lasso = X_train_sig * B_optimal + lasso_intercept;

% �Բ������ݽ���Ԥ��
y_test_pred_lasso = X_test_sig * B_optimal + lasso_intercept;

% ����LASSO�ع��ѵ�����Ͳ��Լ���MSE
mse_train_lasso = mean((y_train - y_train_pred_lasso).^2);
mse_test_lasso = mean((y_test - y_test_pred_lasso).^2);

% �رղ��г�
delete(gcp);

%% ������

disp(['LASSO�ع��ѵ����MSE: ', num2str(mse_train_lasso)]);
disp(['LASSO�ع�Ĳ��Լ�MSE: ', num2str(mse_test_lasso)]);
