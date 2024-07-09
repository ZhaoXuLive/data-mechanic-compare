function [Xi, Xintercept, Sigvars] = sparsifyDynamics_multi(Theta, dXdt)

%% ��������

% ��ȡ���õ� CPU ��������
numCores = feature('numcores');

% �������г�
parpool(numCores);

% ��ʼ��ʱ
tStart = tic;

%% ��������

Xi = zeros(size(Theta,2), 5);
Sigvars = zeros(size(Theta,2), 5);
X_train = Theta;

for i = 1:5
    y_train = dXdt(:, i+5);
    
    disp(i);
    
    % ʹ��parfor�ֿ�����𲽻ع�
    % ����ֿ���
    numBlocks = numCores;
    trainIdx = size(Theta,1);
    blockSize = floor(trainIdx / numBlocks);

    % ��ʼ��ģ�ʹ洢��Ԫ��������������
    mdl_blocks = cell(numBlocks, 1);
    sig_vars_blocks = cell(numBlocks, 1);

    % ʹ�� parfor ѭ�����зֿ�ع�
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
        % ��ȡ������������
        pValues = mdl_blocks{j}.anova.pValue;
        sig_vars_blocks{j} = find(pValues < 0.05); % ѡ��������ˮƽΪ0.05
    end

    % �ϲ�������������
    sig_vars = unique(cell2mat(sig_vars_blocks));
    
    % ��¼�������������ڲ��Լ�
    Sigvars(1, i) = length(sig_vars);
    Sigvars(2:length(sig_vars)+1, i) = sig_vars;

    % ʹ����������ѡ������
    X_train_sig = X_train(:, sig_vars);
    
    % ʹ��LASSO���лع�
    [B, FitInfo] = lasso(X_train_sig, y_train, 'CV', 10);

    % ��ȡ���LASSOģ�͵�ϵ��
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