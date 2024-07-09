%% GPR model
% Perform Gaussian Process regression on discrepancy signal

divideMode = 2;

Xtrain = [Y(:, 3:10) U];
% 'KernelFunction': squaredexponential 平方指数核（默认） matern52 Matérn 5/2核 rationalquadratic 有理二次核
% 'BasisFunction': 指定基函数, 可以使用 'constant'，'linear' 或自定义基函数
% 'exact': 精确方法计算预测结果时，不进行任何近似。适用于数据规模较小的情况
% 'bcd': Bayesian Committee Machine 贝叶斯委员会机将数据分成多个子集，对每个子集分别训练模型并进行预测，最后将结果进行整合
% 'fic': Fully Independent Conditional 完全独立条件方法基于稀疏高斯过程模型，通过选择一组代表性样本（活动集）来近似整个数据集
% 'sd': Subspace Design 子空间设计方法将数据投影到低维子空间，通过在低维子空间中进行计算来提高效率
% 'sr': Sparse Representation 稀疏表示方法选择一组代表性样本（活动集）来近似整个数据集，从而减少计算复杂度
% 'PredictMethod'：指定预测方法，例如 'exact'，'bcd'，'fic'，'sd'。
% 'FitMethod'：指定拟合方法，例如 'exact'，'sr'，'fic'，'sd'
% 'ActiveSetSize' 用于控制稀疏表示方法中的活动集的大小, 从训练数据中选择的一组样本，这些样本用于近似完整数据集的影响
% 注意，ActiveSetSize仅适用于以下拟合方法：sr、fic、sd

% Get the number of available cores
numCores = feature('numcores');

% Determine the number of blocks
numBlocks = min(numCores, size(Xtrain, 1));
blockSize = ceil(size(X, 1) / numBlocks);

for i = 1:10
    if fitMode == 1
        Ytrain = Ef(:,i);
        if divideMode == 1
            gprT = [Xtrain Ytrain];
            tbl = array2table(gprT);
            tbl.Properties.VariableNames = {'Xj1','Yj1','psi1','psi2', 'psi3', 'Vxj1', 'Vyj1', 'psi1d', 'psi2d', 'psi3d', ...
                                          'delta1', 'delta2', 'delta3', 'delta4', 'delta5', 'delta6', 'delta7', 'delta8', ...
                                            'Ef'};                        
    %         discrepGPR{i} = fitrgp(tbl, 'Ef', 'KernelFunction', 'squaredexponential');
            discrepGPR{i} = fitrgp(tbl, 'Ef', 'KernelFunction', 'squaredexponential', ...
                                        'ActiveSetSize',100,'FitMethod','sr','PredictMethod','fic');
        end
        if divideMode == 2
            % Split data into blocks
            blocks = cell(numBlocks, 1);
            for j = 1:numBlocks
                startIdx = (j-1) * blockSize + 1;
                endIdx = min(j * blockSize, size(X, 1));
                blocks{j}.X = Xtrain(startIdx:endIdx, :);
                blocks{j}.Y = Ytrain(startIdx:endIdx);
            end
            % Store models
            models = cell(numBlocks, 1);
            % Train Gaussian Process Regression models in parallel
            parfor k = 1:numBlocks
                gprMdl = fitrgp(blocks{k}.X, blocks{k}.Y);
                models{k} = gprMdl;
            end
        end
        disp(i);
    end
    if fitMode == 2
        Ytrain = dY(:,i);
        if divideMode == 1
            gprT = [Xtrain Ytrain];
            tbl = array2table(gprT);
            tbl.Properties.VariableNames = {'Xj1','Yj1','psi1','psi2', 'psi3', 'Vxj1', 'Vyj1', 'psi1d', 'psi2d', 'psi3d', ...
                                          'delta1', 'delta2', 'delta3', 'delta4', 'delta5', 'delta6', 'delta7', 'delta8', ...
                                            'dY'};                        
    %         discrepGPR{i} = fitrgp(tbl, 'dY', 'KernelFunction', 'squaredexponential');
            discrepGPR{i} = fitrgp(tbl, 'dY', 'KernelFunction', 'squaredexponential', ...
                                        'ActiveSetSize',100,'FitMethod','sr','PredictMethod','fic');
        end
        if divideMode == 2
            % Split data into blocks
            blocks = cell(numBlocks, 1);
            for j = 1:numBlocks
                startIdx = (j-1) * blockSize + 1;
                endIdx = min(j * blockSize, size(X, 1));
                blocks{j}.X = Xtrain(startIdx:endIdx, :);
                blocks{j}.Y = Ytrain(startIdx:endIdx);
            end
            % Store models
            models = cell(numBlocks, 1);
            % Train Gaussian Process Regression models in parallel
            parfor k = 1:numBlocks
                gprMdl = fitrgp(blocks{k}.X, blocks{k}.Y);
                models{k} = gprMdl;
            end
        end
        disp(i);
    end
end

%% 重构训练数据

if dataGenerateMode == 1
    if divideMode == 1
        
        XGS = Maaav_single_GPR_discrep(Y, U, T, discrepGPR, avp, sep, fitMode);
    end
    if divideMode == 2
        XGS = Maaav_single_GPR_discrep_divide(Y, U, T, models, avp, sep, fitMode, numBlocks);
    end
    XO = Maaav_no_discrep(Y, U, T, avp, sep);
end
if divideMode == 1
    XGM = Maaav_GPR_discrep(Y, U, T, discrepGPR, avp, sep, fitMode);
end
if divideMode == 2
    XGM = Maaav_GPR_discrep_divide(Y, U, T, models, avp, sep, fitMode, numBlocks);
end

% Plotting and Calculating Reconstruction of Training Data

if dataGenerateMode == 1
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XGS(sep(i,1):sep(i,2),1), XGS(sep(i,1):sep(i,2),2), 'b--','LineWidth',1), hold on
        plot(XO(sep(i,1):sep(i,2),1), XO(sep(i,1):sep(i,2),2), 'r','LineWidth',1), hold on
        plot(XGM(sep(i,1):sep(i,2),1), XGM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    end

    % data compare
    MSE_NS = 0;
    MSE_NM = 0;
    MSE_GS = 0;
    MSE_GM = 0;
    RMSE_NS = 0;
    RMSE_NM = 0;
    RMSE_GS = 0;
    RMSE_GM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_NS = MSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_NM = MSE_NM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2));
        MSE_GS = MSE_GS + hypot(XGS(i, 1) - Y(i, 1), XGS(i, 2) - Y(i, 2));
        MSE_GM = MSE_GM + hypot(XGM(i, 1) - Y(i, 1), XGM(i, 2) - Y(i, 2));
        RMSE_NS = RMSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_NM = RMSE_NM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2))^2;
        RMSE_GS = RMSE_GS + hypot(XGS(i, 1) - Y(i, 1), XGS(i, 2) - Y(i, 2))^2;
        RMSE_GM = RMSE_GM + hypot(XGM(i, 1) - Y(i, 1), XGM(i, 2) - Y(i, 2))^2;
    end
    MSE_NS = MSE_NS / data_len;
    MSE_NM = MSE_NM / data_len
    MSE_GS = MSE_GS / data_len;
    MSE_GM = MSE_GM / data_len
    RMSE_NS = RMSE_NS / data_len;
    RMSE_NM = RMSE_NM / data_len;
    RMSE_GS = RMSE_GS / data_len;
    RMSE_GM = RMSE_GM / data_len;
end
if dataGenerateMode == 2
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XGM(sep(i,1):sep(i,2),1), XGM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','multi Ideal','multi Learned');
    end

    % data compare
    MSE_NS = 0;
    MSE_GM = 0;
    RMSE_NS = 0;
    RMSE_GM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_NS = MSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_GM = MSE_GM + hypot(XGM(i, 1) - Y(i, 1), XGM(i, 2) - Y(i, 2));
        RMSE_NS = RMSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_GM = RMSE_GM + hypot(XGM(i, 1) - Y(i, 1), XGM(i, 2) - Y(i, 2))^2;
    end
    MSE_NS = MSE_NS / data_len
    MSE_GM = MSE_GM / data_len
    RMSE_NS = RMSE_NS / data_len;
    RMSE_GM = RMSE_GM / data_len;
end

%% 验证集

getValidData;

if dataGenerateMode == 1
    if divideMode == 1
        XGS_v = Maaav_single_GPR_discrep(ValidY, ValidU, ValidT, discrepGPR, avp, sepV, fitMode);
    end
    if divideMode == 2
        XGS_v = Maaav_single_GPR_discrep_divide(ValidY, ValidU, ValidT, models, avp, sepV, fitMode, numBlocks);
    end
    XO_v = Maaav_no_discrep(ValidY, ValidU, ValidT, avp, sepV);
end
if divideMode == 1
    XGM_v = Maaav_GPR_discrep(ValidY, ValidU, ValidT, discrepGPR, avp, sepV, fitMode);
end
if divideMode == 2
    XGM_v = Maaav_GPR_discrep_divide(ValidY, ValidU, ValidT, models, avp, sepV, fitMode, numBlocks);
end

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XGS_v(1:sepV(1,2),1), XGS_v(1:sepV(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(1:sepV(1,2),1), XO_v(1:sepV(1,2),2), 'r','LineWidth',1), hold on
    plot(XGM_v(1:sepV(1,2),1), XGM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XGS_v(sepV(2,1):sepV(2,2),1), XGS_v(sepV(2,1):sepV(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(sepV(2,1):sepV(2,2),1), XO_v(sepV(2,1):sepV(2,2),2), 'r','LineWidth',1), hold on
    plot(XGM_v(sepV(2,1):sepV(2,2),1), XGM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);
    
    % data compare
    MSE_NS_v = 0;
    MSE_NM_v = 0;
    MSE_GS_v = 0;
    MSE_GM_v = 0;
    RMSE_NS_v = 0;
    RMSE_NM_v = 0;
    RMSE_GS_v = 0;
    RMSE_GM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_NS_v = MSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_NM_v = MSE_NM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2));
        MSE_GS_v = MSE_GS_v + hypot(XGS_v(i, 1) - ValidY(i, 1), XGS_v(i, 2) - ValidY(i, 2));
        MSE_GM_v = MSE_GM_v + hypot(XGM_v(i, 1) - ValidY(i, 1), XGM_v(i, 2) - ValidY(i, 2));
        RMSE_NS_v = RMSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_NM_v = RMSE_NM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2))^2;
        RMSE_GS_v = RMSE_GS_v + hypot(XGS_v(i, 1) - ValidY(i, 1), XGS_v(i, 2) - ValidY(i, 2))^2;
        RMSE_GM_v = RMSE_GM_v + hypot(XGM_v(i, 1) - ValidY(i, 1), XGM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_NS_v = MSE_NS_v / data_len_v;
    MSE_NM_v = MSE_NM_v / data_len_v
    MSE_GS_v = MSE_GS_v / data_len_v;
    MSE_GM_v = MSE_GM_v / data_len_v
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_NM_v = RMSE_NM_v / data_len_v;
    RMSE_GS_v = RMSE_GS_v / data_len_v;
    RMSE_GM_v = RMSE_GM_v / data_len_v;
end

if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XGM_v(1:sepV(1,2),1), XGM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XGM_v(sepV(2,1):sepV(2,2),1), XGM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);
    
    % data compare
    MSE_NS_v = 0;
    MSE_GM_v = 0;
    RMSE_NS_v = 0;
    RMSE_GM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_NS_v = MSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_GM_v = MSE_GM_v + hypot(XGM_v(i, 1) - ValidY(i, 1), XGM_v(i, 2) - ValidY(i, 2));
        RMSE_NS_v = RMSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_GM_v = RMSE_GM_v + hypot(XGM_v(i, 1) - ValidY(i, 1), XGM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_NS_v = MSE_NS_v / data_len_v
    MSE_GM_v = MSE_GM_v / data_len_v
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_GM_v = RMSE_GM_v / data_len_v;
end

%% 测试集

getTestData;
 
if dataGenerateMode == 1
    if divideMode == 1
        XGS_t = Maaav_single_GPR_discrep(TestY, TestU, TestT, discrepGPR, avp, sepT, fitMode);
    end
    if divideMode == 2
        XGS_t = Maaav_single_GPR_discrep_divide(TestY, TestU, TestT, models, avp, sepT, fitMode, numBlocks);
    end
    XO_t = Maaav_no_discrep(TestY, TestU, TestT, avp, sepT);
end
if divideMode == 1
    XGM_t = Maaav_GPR_discrep(TestY, TestU, TestT, discrepGPR, avp, sepT, fitMode);
end
if divideMode == 2
    XGM_t = Maaav_GPR_discrep_divide(TestY, TestU, TestT, models, avp, sepT, fitMode, numBlocks);
end

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XGS_t(1:sepT(1,2),1), XGS_t(1:sepT(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(1:sepT(1,2),1), XO_t(1:sepT(1,2),2), 'r','LineWidth',1), hold on
    plot(XGM_t(1:sepT(1,2),1), XGM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XGS_t(sepT(2,1):sepT(2,2),1), XGS_t(sepT(2,1):sepT(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(sepT(2,1):sepT(2,2),1), XO_t(sepT(2,1):sepT(2,2),2), 'r','LineWidth',1), hold on
    plot(XGM_t(sepT(2,1):sepT(2,2),1), XGM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);
    
    % data compare
    MSE_NS_t = 0;
    MSE_NM_t = 0;
    MSE_GS_t = 0;
    MSE_GM_t = 0;
    RMSE_NS_t = 0;
    RMSE_NM_t = 0;
    RMSE_GS_t = 0;
    RMSE_GM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_NS_t = MSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_NM_t = MSE_NM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2));
        MSE_GS_t = MSE_GS_t + hypot(XGS_t(i, 1) - TestY(i, 1), XGS_t(i, 2) - TestY(i, 2));
        MSE_GM_t = MSE_GM_t + hypot(XGM_t(i, 1) - TestY(i, 1), XGM_t(i, 2) - TestY(i, 2));
        RMSE_NS_t = RMSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_NM_t = RMSE_NM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2))^2;
        RMSE_GS_t = RMSE_GS_t + hypot(XGS_t(i, 1) - TestY(i, 1), XGS_t(i, 2) - TestY(i, 2))^2;
        RMSE_GM_t = RMSE_GM_t + hypot(XGM_t(i, 1) - TestY(i, 1), XGM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_NS_t = MSE_NS_t / data_len_t;
    MSE_NM_t = MSE_NM_t / data_len_t
    MSE_GS_t = MSE_GS_t / data_len_t;
    MSE_GM_t = MSE_GM_t / data_len_t
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_NM_t = RMSE_NM_t / data_len_t;
    RMSE_GS_t = RMSE_GS_t / data_len_t;
    RMSE_GM_t = RMSE_GM_t / data_len_t;
end
    
if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XGM_t(1:sepT(1,2),1), XGM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XGM_t(sepT(2,1):sepT(2,2),1), XGM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);

    % data compare
    MSE_NS_t = 0;
    MSE_GM_t = 0;
    RMSE_NS_t = 0;
    RMSE_GM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_NS_t = MSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_GM_t = MSE_GM_t + hypot(XGM_t(i, 1) - TestY(i, 1), XGM_t(i, 2) - TestY(i, 2));
        RMSE_NS_t = RMSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_GM_t = RMSE_GM_t + hypot(XGM_t(i, 1) - TestY(i, 1), XGM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_NS_t = MSE_NS_t / data_len_t
    MSE_GM_t = MSE_GM_t / data_len_t
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_GM_t = RMSE_GM_t / data_len_t;
end
