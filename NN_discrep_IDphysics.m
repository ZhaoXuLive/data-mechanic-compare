%% NN Learn Physics

%% Train Neural Network Parameter

% 神经网络隐藏层数，和各层神经元数目
numNodes = [16 20 20 20 20 20 20 20 20 10];
% feedforwardnet，生成前馈神经网络
% trainbr为训练函数，这里指贝叶斯正则化
% traingd - 梯度下降
% trainrp - 弹性反向传播
net = feedforwardnet(numNodes, 'trainscg');
% 训练轮次
net.trainParam.epochs = 50000;
% 各层传递函数
% Sigmod
net.layers{1}.transferFcn = 'logsig';
% 径向基
net.layers{2}.transferFcn = 'radbas';
% 线性传递函数
net.layers{3}.transferFcn = 'satlins';
% tan-sigmoid函数
net.layers{4}.transferFcn = 'tansig';
% softmax函数
net.layers{5}.transferFcn = 'softmax';
% Elliot Sigmoid函数
net.layers{6}.transferFcn = 'elliotsig';
% Hard Limit函数
net.layers{7}.transferFcn = 'hardlim';
% Saturating Linear函数
net.layers{8}.transferFcn = 'satlin';
% Symmetric Saturating Linear函数
net.layers{9}.transferFcn = 'tribas';
% Triangular Basis Function函数
net.layers{10}.transferFcn = 'purelin';
% % Sigmod
% net.layers{1}.transferFcn = 'logsig';
% % 径向基
% net.layers{2}.transferFcn = 'radbas';
% % 线性传递函数
% net.layers{3}.transferFcn = 'purelin';
% % tan-sigmoid函数
% net.layers{4}.transferFcn = 'tansig';
% % softmax函数
% net.layers{5}.transferFcn = 'softmax';
% % Elliot Sigmoid函数
% net.layers{6}.transferFcn = 'elliotsig';
% % Hard Limit函数
% net.layers{7}.transferFcn = 'hardlim';
% % Saturating Linear函数
% net.layers{8}.transferFcn = 'satlin';
% % Symmetric Saturating Linear函数
% net.layers{9}.transferFcn = 'satlins';
% % Triangular Basis Function函数
% net.layers{10}.transferFcn = 'tribas';
% % Sigmod
% net.layers{11}.transferFcn = 'logsig';
% % 径向基
% net.layers{12}.transferFcn = 'radbas';
% % 线性传递函数
% net.layers{13}.transferFcn = 'purelin';
% % tan-sigmoid函数
% net.layers{14}.transferFcn = 'tansig';
% % softmax函数
% net.layers{15}.transferFcn = 'softmax';

%% Train process

if fitMode == 1
%     TrainX = Y(:, 3:10)';
    TrainX = [Y(:, 3:10) U]';
    TrainT = Ef';
end
if fitMode == 2
%     TrainX = Y(:, 3:10)';
    TrainX = [Y(:, 3:10) U]';
    TrainT = dY';
end

% net = train(net,x,t); x 输入值 t 相关联的目标输出值，这里应为n*size
% net = fitnet(numNodes,'trainbr'); 函数拟合神经网络
% net = train(net, TrainX, TrainT);
net = train(net, TrainX, TrainT, 'useParallel', 'yes', 'useGPU', 'only');

% % 训练集网络估计结果
% result = net(TrainX);
% % 可视化
% view(net)
% % 评估经过训练的网络的性能。默认性能函数是均方误差
% perf = perform(net, result, TrainT);

%% 重建训练数据

if dataGenerateMode == 1
    XNS = Maaav_single_NN_discrep(Y, U, T, net, avp, sep, fitMode);
end
XO = Maaav_no_discrep(Y, U, T, avp, sep);
XNM = Maaav_NN_discrep(Y, U, T, net, avp, sep, fitMode);

% Plotting and Calculating Reconstruction of Training Data

if dataGenerateMode == 1
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XNS(sep(i,1):sep(i,2),1), XNS(sep(i,1):sep(i,2),2), 'b--','LineWidth',1), hold on
        plot(XO(sep(i,1):sep(i,2),1), XO(sep(i,1):sep(i,2),2), 'r','LineWidth',1), hold on
        plot(XNM(sep(i,1):sep(i,2),1), XNM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    end

    % data compare
    MSE_LS = 0;
    MSE_LM = 0;
    MSE_NS = 0;
    MSE_NM = 0;
    RMSE_LS = 0;
    RMSE_LM = 0;
    RMSE_NS = 0;
    RMSE_NM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_LS = MSE_LS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_LM = MSE_LM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2));
        MSE_NS = MSE_NS + hypot(XNS(i, 1) - Y(i, 1), XNS(i, 2) - Y(i, 2));
        MSE_NM = MSE_NM + hypot(XNM(i, 1) - Y(i, 1), XNM(i, 2) - Y(i, 2));
        RMSE_LS = RMSE_LS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_LM = RMSE_LM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2))^2;
        RMSE_NS = RMSE_NS + hypot(XNS(i, 1) - Y(i, 1), XNS(i, 2) - Y(i, 2))^2;
        RMSE_NM = RMSE_NM + hypot(XNM(i, 1) - Y(i, 1), XNM(i, 2) - Y(i, 2))^2;
    end
    MSE_LS = MSE_LS / data_len;
    MSE_LM = MSE_LM / data_len
    MSE_NS = MSE_NS / data_len;
    MSE_NM = MSE_NM / data_len
    RMSE_LS = RMSE_LS / data_len;
    RMSE_LM = RMSE_LM / data_len;
    RMSE_NS = RMSE_NS / data_len;
    RMSE_NM = RMSE_NM / data_len;
end
if dataGenerateMode == 2
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XNM(sep(i,1):sep(i,2),1), XNM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    end

    % data compare
    MSE_LS = 0;
    MSE_NM = 0;
    RMSE_LS = 0;
    RMSE_NM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_LS = MSE_LS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_NM = MSE_NM + hypot(XNM(i, 1) - Y(i, 1), XNM(i, 2) - Y(i, 2));
        RMSE_LS = RMSE_LS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_NM = RMSE_NM + hypot(XNM(i, 1) - Y(i, 1), XNM(i, 2) - Y(i, 2))^2;
    end
    MSE_LS = MSE_LS / data_len
    MSE_NM = MSE_NM / data_len
    RMSE_LS = RMSE_LS / data_len;
    RMSE_NM = RMSE_NM / data_len;
end

%% 验证集

getValidData;

if dataGenerateMode == 1
    XNS_v = Maaav_single_NN_discrep(ValidY, ValidU, ValidT, net, avp, sepV, fitMode);
end
XO_v = Maaav_no_discrep(ValidY, ValidU, ValidT, avp, sepV);
XNM_v = Maaav_NN_discrep(ValidY, ValidU, ValidT, net, avp, sepV, fitMode);

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XNS_v(1:sepV(1,2),1), XNS_v(1:sepV(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(1:sepV(1,2),1), XO_v(1:sepV(1,2),2), 'r','LineWidth',1), hold on
    plot(XNM_v(1:sepV(1,2),1), XNM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XNS_v(sepV(2,1):sepV(2,2),1), XNS_v(sepV(2,1):sepV(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(sepV(2,1):sepV(2,2),1), XO_v(sepV(2,1):sepV(2,2),2), 'r','LineWidth',1), hold on
    plot(XNM_v(sepV(2,1):sepV(2,2),1), XNM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);
    
    % data compare
    MSE_LS_v = 0;
    MSE_LM_v = 0;
    MSE_NS_v = 0;
    MSE_NM_v = 0;
    RMSE_LS_v = 0;
    RMSE_LM_v = 0;
    RMSE_NS_v = 0;
    RMSE_NM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_LS_v = MSE_LS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_LM_v = MSE_LM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2));
        MSE_NS_v = MSE_NS_v + hypot(XNS_v(i, 1) - ValidY(i, 1), XNS_v(i, 2) - ValidY(i, 2));
        MSE_NM_v = MSE_NM_v + hypot(XNM_v(i, 1) - ValidY(i, 1), XNM_v(i, 2) - ValidY(i, 2));
        RMSE_LS_v = RMSE_LS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_LM_v = RMSE_LM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2))^2;
        RMSE_NS_v = RMSE_NS_v + hypot(XNS_v(i, 1) - ValidY(i, 1), XNS_v(i, 2) - ValidY(i, 2))^2;
        RMSE_NM_v = RMSE_NM_v + hypot(XNM_v(i, 1) - ValidY(i, 1), XNM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_LS_v = MSE_LS_v / data_len_v;
    MSE_LM_v = MSE_LM_v / data_len_v
    MSE_NS_v = MSE_NS_v / data_len_v;
    MSE_NM_v = MSE_NM_v / data_len_v
    RMSE_LS_v = RMSE_LS_v / data_len_v;
    RMSE_LM_v = RMSE_LM_v / data_len_v;
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_NM_v = RMSE_NM_v / data_len_v;
end

if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XNM_v(1:sepV(1,2),1), XNM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XNM_v(sepV(2,1):sepV(2,2),1), XNM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);
    
    % data compare
    MSE_LS_v = 0;
    MSE_NM_v = 0;
    RMSE_LS_v = 0;
    RMSE_NM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_LS_v = MSE_LS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_NM_v = MSE_NM_v + hypot(XNM_v(i, 1) - ValidY(i, 1), XNM_v(i, 2) - ValidY(i, 2));
        RMSE_LS_v = RMSE_LS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_NM_v = RMSE_NM_v + hypot(XNM_v(i, 1) - ValidY(i, 1), XNM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_LS_v = MSE_LS_v / data_len_v
    MSE_NM_v = MSE_NM_v / data_len_v
    RMSE_LS_v = RMSE_LS_v / data_len_v;
    RMSE_NM_v = RMSE_NM_v / data_len_v;
end

%% 测试集

getTestData;

if dataGenerateMode == 1
    XNS_t = Maaav_single_NN_discrep(TestY, TestU, TestT, net, avp, sepT, fitMode);
end
XO_t = Maaav_no_discrep(TestY, TestU, TestT, avp, sepT);
XNM_t = Maaav_NN_discrep(TestY, TestU, TestT, net, avp, sepT, fitMode);

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XNS_t(1:sepT(1,2),1), XNS_t(1:sepT(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(1:sepT(1,2),1), XO_t(1:sepT(1,2),2), 'r','LineWidth',1), hold on
    plot(XNM_t(1:sepT(1,2),1), XNM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XNS_t(sepT(2,1):sepT(2,2),1), XNS_t(sepT(2,1):sepT(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(sepT(2,1):sepT(2,2),1), XO_t(sepT(2,1):sepT(2,2),2), 'r','LineWidth',1), hold on
    plot(XNM_t(sepT(2,1):sepT(2,2),1), XNM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);
    
    % data compare
    MSE_LS_t = 0;
    MSE_LM_t = 0;
    MSE_NS_t = 0;
    MSE_NM_t = 0;
    RMSE_LS_t = 0;
    RMSE_LM_t = 0;
    RMSE_NS_t = 0;
    RMSE_NM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_LS_t = MSE_LS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_LM_t = MSE_LM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2));
        MSE_NS_t = MSE_NS_t + hypot(XNS_t(i, 1) - TestY(i, 1), XNS_t(i, 2) - TestY(i, 2));
        MSE_NM_t = MSE_NM_t + hypot(XNM_t(i, 1) - TestY(i, 1), XNM_t(i, 2) - TestY(i, 2));
        RMSE_LS_t = RMSE_LS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_LM_t = RMSE_LM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2))^2;
        RMSE_NS_t = RMSE_NS_t + hypot(XNS_t(i, 1) - TestY(i, 1), XNS_t(i, 2) - TestY(i, 2))^2;
        RMSE_NM_t = RMSE_NM_t + hypot(XNM_t(i, 1) - TestY(i, 1), XNM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_LS_t = MSE_LS_t / data_len_t;
    MSE_LM_t = MSE_LM_t / data_len_t
    MSE_NS_t = MSE_NS_t / data_len_t;
    MSE_NM_t = MSE_NM_t / data_len_t
    RMSE_LS_t = RMSE_LS_t / data_len_t;
    RMSE_LM_t = RMSE_LM_t / data_len_t;
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_NM_t = RMSE_NM_t / data_len_t;
end
    
if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XNM_t(1:sepT(1,2),1), XNM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XNM_t(sepT(2,1):sepT(2,2),1), XNM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);

    % data compare
    MSE_LS_t = 0;
    MSE_NM_t = 0;
    RMSE_LS_t = 0;
    RMSE_NM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_LS_t = MSE_LS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_NM_t = MSE_NM_t + hypot(XNM_t(i, 1) - TestY(i, 1), XNM_t(i, 2) - TestY(i, 2));
        RMSE_LS_t = RMSE_LS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_NM_t = RMSE_NM_t + hypot(XNM_t(i, 1) - TestY(i, 1), XNM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_LS_t = MSE_LS_t / data_len_t
    MSE_NM_t = MSE_NM_t / data_len_t
    RMSE_LS_t = RMSE_LS_t / data_len_t;
    RMSE_NM_t = RMSE_NM_t / data_len_t;
end
