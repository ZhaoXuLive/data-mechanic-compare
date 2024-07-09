% 1. ��������
n = 1000; % ��������
p = 100; % ��������
X = randn(n, p); % ���������������
true_beta = [3; -2; 0; 0; 1; zeros(p-5, 1)]; % ��ʵ�ع�ϵ��������Ϊ0
y = X * true_beta + randn(n, 1); % �����������Ӧ����

% 2. ���ݱ�׼�����𲽻ع�ǰ����
X_mean = mean(X);
X_std = std(X);
X_std(X_std == 0) = 1; % ���������
X = (X - X_mean) ./ X_std;
y_mean = mean(y);
y = y - y_mean;

%% lasso

% % 3. ʹ��Lasso�ع���б���ѡ��
% [B, FitInfo] = lasso(X, y, 'CV', 10);
% lassoLambda = FitInfo.Lambda1SE; % ʹ��1-SE����ѡ��Lambda
% lassoBeta = B(:, FitInfo.Index1SE);
% 
% % ��ӡ���
% selectedLassoVars = find(lassoBeta ~= 0);
% fprintf('Lassoѡ��ı�������: %s\n', mat2str(selectedLassoVars));
% fprintf('Lassoѡ���ϵ��: \n');
% disp(lassoBeta(lassoBeta ~= 0));
% 
% % ʹ�� LassoBeta Ԥ����Ӧ����
% y_pred = X * lassoBeta;
% 
% % ���������� (MSE)
% MSE = mean((y - y_pred).^2);
% 
% % ������
% fprintf('Lasso �ع�ľ������ (MSE): %.4f\n', MSE);

%% lasso + stepwise

% ��ȡ���õ� CPU ��������
numCores = feature('numcores');

% �������г�
parpool(numCores);

tStart = tic;

% ����ֿ�������Ҫ�������ú�������
numBlocks = numCores;
blockSize = floor(n / numBlocks);

% ��ʼ��ģ�ʹ洢��Ԫ
mdl = cell(numBlocks, 1);

% ����PCA��ά
[coeff, score, ~, ~, explained] = pca(X);

% ѡ�����90%��������ɷ�
numComponents = find(cumsum(explained) >= 90, 1);
X_reduced = score(:, 1:numComponents);

% ʹ�� parfor ѭ�����зֿ�ع�
parfor i = 1:numBlocks
    startIdx = (i-1) * blockSize + 1;
    if i == numBlocks
        endIdx = n; % ȷ�����һ�����������ʣ������
    else
        endIdx = i * blockSize;
    end
    X_block = X_reduced(startIdx:endIdx, :);
    y_block = y(startIdx:endIdx);
    mdl{i} = stepwiselm(X_block, y_block, 'linear', 'Upper', 'linear', 'Lower', 'constant', ...
                   'PEnter', 0.05, 'PRemove', 0.10, 'Verbose', 0);
end

% % 3. ʹ���𲽻ع���г�������ѡ��ȷ��ֻ����һ���Ա���
% mdl = stepwiselm(X, y, 'linear', 'Upper', 'linear', 'Lower', 'constant', ...
%                    'PEnter', 0.05, 'PRemove', 0.10, 'Verbose', 0);

tEnd = toc(tStart)

% �رղ��г�
delete(gcp);
                      
% 4. ѡ���𲽻ع��е���������
coeffs = mdl.Coefficients;
p_values = coeffs.pValue(2:end); % �ų��ؾ����pֵ
selectedVars = find(p_values < 0.10); % ѡ����������

% ����Ƿ�����������
if ~isempty(selectedVars)
    % 5. ʹ��Lasso�ع����𲽻ع�ѡ��ı��������Ͻ��н�һ���Ż�
    X_selected = X(:, selectedVars);
    [B, FitInfo] = lasso(X_selected, y, 'CV', 10);
    lassoLambda = FitInfo.Lambda1SE; % ʹ��1-SE����ѡ��Lambda
    lassoBeta = B(:, FitInfo.Index1SE);

    % ��ӡ���
    selectedLassoVars = selectedVars(lassoBeta ~= 0);
    fprintf('Lassoѡ��ı�������: %s\n', mat2str(selectedLassoVars));
    fprintf('Lassoѡ���ϵ��: \n');
    disp(lassoBeta(lassoBeta ~= 0));
    
    % ʹ�� LassoBeta Ԥ����Ӧ����
    allBeta = zeros(p, 1);
    allBeta(1:size(lassoBeta(:, 1), 1), :)  = lassoBeta(:, :);
    y_pred = X * allBeta;

    % ���������� (MSE)
    MSE = mean((y - y_pred).^2);

    % ������
    fprintf('StepWise + Lasso �ع�ľ������ (MSE): %.4f\n', MSE);
    
%     % ����Lasso·��ͼ
%     figure;
%     lassoPlot(B, FitInfo, 'PlotType', 'Lambda', 'XScale', 'log');
%     xlabel('Lambda');
%     ylabel('�ع�ϵ��');
%     title('Lasso·��ͼ');
%     grid on;
else
    disp('�𲽻ع�δѡ���κα�����');
end

% tEnd = toc(tStart)
