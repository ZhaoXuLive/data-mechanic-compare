%% DMD model

optMode = 1;
% 优化迭代次数
itern = 10;
% 学习率
learnv = 0.001;

if fitMode == 1
    X_train = Ef(1:end-1, :)';
    X_train_prime = Ef(2:end, :)';
    U_train = U(1:end-1, :)';
    
    % 矩阵大小
    [Xtn, Xtm] = size(X_train);
    Xtp = size(U_train, 1);

    % 构建扩展矩阵 Omega
    Omega = [X_train; U_train];
    % 计算 Omega 的伪逆
    Omega_pseudo_inverse = pinv(Omega);
    % 计算 K 矩阵
    K = X_train_prime * Omega_pseudo_inverse;

    % 从 K 中分解出 A 和 B
    A = K(:, 1:Xtn);
    B = K(:, Xtn+1:end);

    % 奇异值分解 X
    [U_X, Sigma_X, V_X] = svd(X_train, 'econ');
    % 计算降维后的 A_tilde
    A_tilde = U_X' * A * U_X;
    % 特征分解 A_tilde
    [W, Lambda] = eig(A_tilde);
    % 计算 DMD 模态
    Phi = X_train_prime * V_X / Sigma_X * W;
    b = Phi \ X_train(:,1);
    
    if optMode == 2
        Phi_optimized = optimizeModes(Phi, Lambda, b, X_train, Xtm, itern, learnv);
        Phi = Phi_optimized;
    end
end
if fitMode == 2
    X_train = Y(1:end-1, :)';
    X_train_prime = Y(2:end, :)';
    U_train = U(1:end-1, :)';
    
    % 矩阵大小
    [Xtn, Xtm] = size(X_train);
    Xtp = size(U_train, 1);

    % 构建扩展矩阵 Omega
    Omega = [X_train; U_train];
    % 计算 Omega 的伪逆
    Omega_pseudo_inverse = pinv(Omega);
    % 计算 K 矩阵
    K = X_train_prime * Omega_pseudo_inverse;

    % 从 K 中分解出 A 和 B
    A = K(:, 1:Xtn);
    B = K(:, Xtn+1:end);

    % 奇异值分解 X
    [U_X, Sigma_X, V_X] = svd(X_train, 'econ');
    % 计算降维后的 A_tilde
    A_tilde = U_X' * A * U_X;
    % 特征分解 A_tilde
    [W, Lambda] = eig(A_tilde);
    % 计算 DMD 模态
    Phi = X_train_prime * V_X / Sigma_X * W;
    b = Phi \ X_train(:,1);
    
    if optMode == 2
        Phi_optimized = optimizeModes(Phi, Lambda, b, X_train, Xtm, itern, learnv);
        Phi = Phi_optimized;
    end
end

%% 重构训练数据

if dataGenerateMode == 1
    XDS = Maaav_single_DMD_discrep(Phi, Lambda, Y, U, T, avp, sep, fitMode);
    XO = Maaav_no_discrep(Y, U, T, avp, sep);
end
% XDM = Maaav_DMD_discrep(Phi, Lambda, Y, U, T, avp, sep, fitMode);
XDM = Maaav_DMD_discrep_fixb(Phi, Lambda, b, Y, U, T, avp, sep, fitMode);


%% Plotting and Calculating Reconstruction of Training Data

if dataGenerateMode == 1
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XDS(sep(i,1):sep(i,2),1), XDS(sep(i,1):sep(i,2),2), 'b--','LineWidth',1), hold on
        plot(XO(sep(i,1):sep(i,2),1), XO(sep(i,1):sep(i,2),2), 'r','LineWidth',1), hold on
        plot(XDM(sep(i,1):sep(i,2),1), XDM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    end

    % data compare
    MSE_NS = 0;
    MSE_NM = 0;
    MSE_DS = 0;
    MSE_DM = 0;
    RMSE_NS = 0;
    RMSE_NM = 0;
    RMSE_DS = 0;
    RMSE_DM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_NS = MSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_NM = MSE_NM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2));
        MSE_DS = MSE_DS + hypot(XDS(i, 1) - Y(i, 1), XDS(i, 2) - Y(i, 2));
        MSE_DM = MSE_DM + hypot(XDM(i, 1) - Y(i, 1), XDM(i, 2) - Y(i, 2));
        RMSE_NS = RMSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_NM = RMSE_NM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2))^2;
        RMSE_DS = RMSE_DS + hypot(XDS(i, 1) - Y(i, 1), XDS(i, 2) - Y(i, 2))^2;
        RMSE_DM = RMSE_DM + hypot(XDM(i, 1) - Y(i, 1), XDM(i, 2) - Y(i, 2))^2;
    end
    MSE_NS = MSE_NS / data_len;
    MSE_NM = MSE_NM / data_len
    MSE_DS = MSE_DS / data_len;
    MSE_DM = MSE_DM / data_len
    RMSE_NS = RMSE_NS / data_len;
    RMSE_NM = RMSE_NM / data_len;
    RMSE_DS = RMSE_DS / data_len;
    RMSE_DM = RMSE_DM / data_len;
end
if dataGenerateMode == 2
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XDM(sep(i,1):sep(i,2),1), XDM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','multi Ideal','multi Learned');
    end

    % data compare
    MSE_NS = 0;
    MSE_DM = 0;
    RMSE_NS = 0;
    RMSE_DM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_NS = MSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_DM = MSE_DM + hypot(XDM(i, 1) - Y(i, 1), XDM(i, 2) - Y(i, 2));
        RMSE_NS = RMSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_DM = RMSE_DM + hypot(XDM(i, 1) - Y(i, 1), XDM(i, 2) - Y(i, 2))^2;
    end
    MSE_NS = MSE_NS / data_len
    MSE_DM = MSE_DM / data_len
    RMSE_NS = RMSE_NS / data_len;
    RMSE_DM = RMSE_DM / data_len;
end

%% 验证集

getValidData;

if dataGenerateMode == 1
    XDS_v = Maaav_single_DMD_discrep(Phi, Lambda, ValidY, ValidU, ValidT, avp, sepV, fitMode);
    XO_v = Maaav_no_discrep(ValidY, ValidU, ValidT, avp, sepV);
end
% XDM_v = Maaav_DMD_discrep(Phi, Lambda, ValidY, ValidU, ValidT, avp, sepV, fitMode);
XDM_v = Maaav_DMD_discrep_fixb(Phi, Lambda, b, ValidY, ValidU, ValidT, avp, sepV, fitMode);

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XDS_v(1:sepV(1,2),1), XDS_v(1:sepV(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(1:sepV(1,2),1), XO_v(1:sepV(1,2),2), 'r','LineWidth',1), hold on
    plot(XDM_v(1:sepV(1,2),1), XDM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XDS_v(sepV(2,1):sepV(2,2),1), XDS_v(sepV(2,1):sepV(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(sepV(2,1):sepV(2,2),1), XO_v(sepV(2,1):sepV(2,2),2), 'r','LineWidth',1), hold on
    plot(XDM_v(sepV(2,1):sepV(2,2),1), XDM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);
    
    % data compare
    MSE_NS_v = 0;
    MSE_NM_v = 0;
    MSE_DS_v = 0;
    MSE_DM_v = 0;
    RMSE_NS_v = 0;
    RMSE_NM_v = 0;
    RMSE_DS_v = 0;
    RMSE_DM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_NS_v = MSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_NM_v = MSE_NM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2));
        MSE_DS_v = MSE_DS_v + hypot(XDS_v(i, 1) - ValidY(i, 1), XDS_v(i, 2) - ValidY(i, 2));
        MSE_DM_v = MSE_DM_v + hypot(XDM_v(i, 1) - ValidY(i, 1), XDM_v(i, 2) - ValidY(i, 2));
        RMSE_NS_v = RMSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_NM_v = RMSE_NM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2))^2;
        RMSE_DS_v = RMSE_DS_v + hypot(XDS_v(i, 1) - ValidY(i, 1), XDS_v(i, 2) - ValidY(i, 2))^2;
        RMSE_DM_v = RMSE_DM_v + hypot(XDM_v(i, 1) - ValidY(i, 1), XDM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_NS_v = MSE_NS_v / data_len_v;
    MSE_NM_v = MSE_NM_v / data_len_v
    MSE_DS_v = MSE_DS_v / data_len_v;
    MSE_DM_v = MSE_DM_v / data_len_v
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_NM_v = RMSE_NM_v / data_len_v;
    RMSE_DS_v = RMSE_DS_v / data_len_v;
    RMSE_DM_v = RMSE_DM_v / data_len_v;
end

if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XDM_v(1:sepV(1,2),1), XDM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XDM_v(sepV(2,1):sepV(2,2),1), XDM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);
    
    % data compare
    MSE_NS_v = 0;
    MSE_DM_v = 0;
    RMSE_NS_v = 0;
    RMSE_DM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_NS_v = MSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_DM_v = MSE_DM_v + hypot(XDM_v(i, 1) - ValidY(i, 1), XDM_v(i, 2) - ValidY(i, 2));
        RMSE_NS_v = RMSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_DM_v = RMSE_DM_v + hypot(XDM_v(i, 1) - ValidY(i, 1), XDM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_NS_v = MSE_NS_v / data_len_v
    MSE_DM_v = MSE_DM_v / data_len_v
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_DM_v = RMSE_DM_v / data_len_v;
end

%% 测试集

getTestData;

if dataGenerateMode == 1
    XDS_t = Maaav_single_DMD_discrep(Phi, Lambda, TestY, TestU, TestT, avp, sepT, fitMode);
    XO_t = Maaav_no_discrep(TestY, TestU, TestT, avp, sepT);
end
% XDM_t = Maaav_DMD_discrep(Phi, Lambda, TestY, TestU, TestT, avp, sepT, fitMode);
XDM_t = Maaav_DMD_discrep_fixb(Phi, Lambda, b, TestY, TestU, TestT, avp, sepT, fitMode);

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XDS_t(1:sepT(1,2),1), XDS_t(1:sepT(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(1:sepT(1,2),1), XO_t(1:sepT(1,2),2), 'r','LineWidth',1), hold on
    plot(XDM_t(1:sepT(1,2),1), XDM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XDS_t(sepT(2,1):sepT(2,2),1), XDS_t(sepT(2,1):sepT(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(sepT(2,1):sepT(2,2),1), XO_t(sepT(2,1):sepT(2,2),2), 'r','LineWidth',1), hold on
    plot(XDM_t(sepT(2,1):sepT(2,2),1), XDM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);
    
    % data compare
    MSE_NS_t = 0;
    MSE_NM_t = 0;
    MSE_DS_t = 0;
    MSE_DM_t = 0;
    RMSE_NS_t = 0;
    RMSE_NM_t = 0;
    RMSE_DS_t = 0;
    RMSE_DM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_NS_t = MSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_NM_t = MSE_NM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2));
        MSE_DS_t = MSE_DS_t + hypot(XDS_t(i, 1) - TestY(i, 1), XDS_t(i, 2) - TestY(i, 2));
        MSE_DM_t = MSE_DM_t + hypot(XDM_t(i, 1) - TestY(i, 1), XDM_t(i, 2) - TestY(i, 2));
        RMSE_NS_t = RMSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_NM_t = RMSE_NM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2))^2;
        RMSE_DS_t = RMSE_DS_t + hypot(XDS_t(i, 1) - TestY(i, 1), XDS_t(i, 2) - TestY(i, 2))^2;
        RMSE_DM_t = RMSE_DM_t + hypot(XDM_t(i, 1) - TestY(i, 1), XDM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_NS_t = MSE_NS_t / data_len_t;
    MSE_NM_t = MSE_NM_t / data_len_t
    MSE_DS_t = MSE_DS_t / data_len_t;
    MSE_DM_t = MSE_DM_t / data_len_t
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_NM_t = RMSE_NM_t / data_len_t;
    RMSE_DS_t = RMSE_DS_t / data_len_t;
    RMSE_DM_t = RMSE_DM_t / data_len_t;
end
    
if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XDM_t(1:sepT(1,2),1), XDM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XDM_t(sepT(2,1):sepT(2,2),1), XDM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);

    % data compare
    MSE_NS_t = 0;
    MSE_DM_t = 0;
    RMSE_NS_t = 0;
    RMSE_DM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_NS_t = MSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_DM_t = MSE_DM_t + hypot(XDM_t(i, 1) - TestY(i, 1), XDM_t(i, 2) - TestY(i, 2));
        RMSE_NS_t = RMSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_DM_t = RMSE_DM_t + hypot(XDM_t(i, 1) - TestY(i, 1), XDM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_NS_t = MSE_NS_t / data_len_t
    MSE_DM_t = MSE_DM_t / data_len_t
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_DM_t = RMSE_DM_t / data_len_t;
end
