%% Sparse Identification of Nonlinear Dynamics for discrepancy modeling & learning missing physics

% pool Data  (i.e., build library of nonlinear time series)

% 去除无用状态
state = Y(:, 3:10);

% 在修改训练集之后无有效结果
if poolMode == 1
    Theta = poolData_withU(state, U, polyorder);
end
if poolMode == 2
    Theta = poolData_new(state, U, polyorder);
end
if poolMode == 3
    augState = getAugState(state, sep);
    Theta = poolData_augNew(augState, U, polyorder);
end

%% Compute sparse regression for discrepancy g(x): sequential least squares

Xintercept = zeros(5, 1);
if fitMode == 1
    if sparsifyMode == 1
        % sequential least squares
        Xi = sparsifyDynamics_sls(Theta, Ef, lambda);
    end
    if sparsifyMode == 2
        % sequential least squares with library normalization
        Xi = sparsifyDynamics_nsls(Theta, Ef, lambda);
    end
    if sparsifyMode == 3
        % lasso
        Xi = sparsifyDynamics_lasso(Theta, Ef);
    end
    if sparsifyMode == 4
        % sequential least squares with Stepwise Sparse Regressor
        Xi = sparsifyDynamics_ssr(Theta, Ef);
    end
    if sparsifyMode == 5
        % ridge regression
        [Xi, Xintercept] = sparsifyDynamics_ridge(Theta, Ef);
    end  
    if sparsifyMode == 6
        % Elastic regression
        [Xi, Xintercept] = sparsifyDynamics_elastic(Theta, Ef);
    end
    if sparsifyMode == 7
        % pca + parfor + stepwise + lasso
        [Xi, Xintercept, Sigvars] = sparsifyDynamics_multi(Theta, Ef);
    end
%     Error = Ef(:, 6:10) - (Theta * Xi + Xintercept');
end
if fitMode == 2
    if sparsifyMode == 1
        Xi = sparsifyDynamics_sls(Theta, dY, lambda);
    end
    if sparsifyMode == 2
        Xi = sparsifyDynamics_nsls(Theta, dY, lambda);
    end
    if sparsifyMode == 3
        Xi = sparsifyDynamics_lasso(Theta, dY);
    end
    if sparsifyMode == 4
        % sequential least squares with Stepwise Sparse Regressor
        Xi = sparsifyDynamics_ssr(Theta, dY);
    end
    if sparsifyMode == 5
        % ridge regression
        [Xi, Xintercept] = sparsifyDynamics_ridge(Theta, dY);
    end  
    if sparsifyMode == 6
        % Elastic regression
        [Xi, Xintercept] = sparsifyDynamics_elastic(Theta, dY);
    end  
    if sparsifyMode == 7
        % pca + parfor + stepwise + lasso
        [Xi, Xintercept, Sigvars] = sparsifyDynamics_multi(Theta, dY);
    end
%     Error = dY(:, 6:10) - (Theta * Xi + Xintercept');
end

% figure,
% title('Error of learned and discrepency Data')
% for i = 1:5
%     subplot(1,5,i)
%     plot(T(1:end), Error(1:end,i),'k','Linewidth',[2]);
% end

%% Reconstruction of Training Data

% signy = {'psi1', 'psi2', 'psi3', 'vx', 'vy', 'psi1d', 'psi2d', 'psi3d',...
%          'delta1', 'delta2', 'delta3', 'delta4', 'delta5', 'delta6', 'delta7', 'delta8'};
% poolDataLIST_withU(signy, Xi, n-2, un, polyorder);

% 如果为单步产生的数据集，
% 在差异模型基础下，那么XS、XO、XM都可以计算，
% XS用于和X、Y比较，比较的是修正后的模型在单步预测下的精度区别；
% XM用于和XO、Y比较，比较的是修正后的模型在多步预测下的精度区别；
% 在完整模型基础下，也应该都计算，因为比较的东西还是一样的呀！

% 如果为多步产生的数据集，
% 在差异模型基础下，只需要XM，比较的是Y、X、XM；
% 在完整模型基础下，只需要XM，比较的是Y、X、XM；

if sparsifyMode == 7
    if dataGenerateMode == 1
        XS = Maaav_single_SINDy_discrep_multi(Y, U, T, Xi, Xintercept, avp, polyorder, sep, poolMode, fitMode, Sigvars);
        XO = Maaav_no_discrep(Y, U, T, avp, sep);
    end
    XM = Maaav_SINDy_discrep_multi(Y, U, T, Xi, Xintercept, avp, polyorder, sep, poolMode, fitMode, Sigvars);
else
    if dataGenerateMode == 1
        XS = Maaav_single_SINDy_discrep(Y, U, T, Xi, Xintercept, avp, polyorder, sep, poolMode, fitMode);
        XO = Maaav_no_discrep(Y, U, T, avp, sep);
    end
    XM = Maaav_SINDy_discrep(Y, U, T, Xi, Xintercept, avp, polyorder, sep, poolMode, fitMode);
end

%% Plotting and Calculating Reconstruction of Training Data

if dataGenerateMode == 1
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XS(sep(i,1):sep(i,2),1), XS(sep(i,1):sep(i,2),2), 'b--','LineWidth',1), hold on
        plot(XO(sep(i,1):sep(i,2),1), XO(sep(i,1):sep(i,2),2), 'r','LineWidth',1), hold on
        plot(XM(sep(i,1):sep(i,2),1), XM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    end
    
    % data compare
    MSE_NS = 0;
    MSE_NM = 0;
    MSE_SS = 0;
    MSE_SM = 0;
    RMSE_NS = 0;
    RMSE_NM = 0;
    RMSE_SS = 0;
    RMSE_SM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_NS = MSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_NM = MSE_NM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2));
        MSE_SS = MSE_SS + hypot(XS(i, 1) - Y(i, 1), XS(i, 2) - Y(i, 2));
        MSE_SM = MSE_SM + hypot(XM(i, 1) - Y(i, 1), XM(i, 2) - Y(i, 2));
        RMSE_NS = RMSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_NM = RMSE_NM + hypot(XO(i, 1) - Y(i, 1), XO(i, 2) - Y(i, 2))^2;
        RMSE_SS = RMSE_SS + hypot(XS(i, 1) - Y(i, 1), XS(i, 2) - Y(i, 2))^2;
        RMSE_SM = RMSE_SM + hypot(XM(i, 1) - Y(i, 1), XM(i, 2) - Y(i, 2))^2;
    end
    MSE_NS = MSE_NS / data_len;
    MSE_NM = MSE_NM / data_len
    MSE_SS = MSE_SS / data_len;
    MSE_SM = MSE_SM / data_len
    RMSE_NS = RMSE_NS / data_len;
    RMSE_NM = RMSE_NM / data_len;
    RMSE_SS = RMSE_SS / data_len;
    RMSE_SM = RMSE_SM / data_len;
end
if dataGenerateMode == 2
    figure,
    for i = 1:length(sep(:, 1))
        subplot(1,6,i);
        plot(Y(sep(i,1):sep(i,2),1), Y(sep(i,1):sep(i,2),2), 'k','LineWidth',1), hold on
        plot(X(sep(i,1):sep(i,2),1), X(sep(i,1):sep(i,2),2), 'r--','LineWidth',1), hold on
        plot(XM(sep(i,1):sep(i,2),1), XM(sep(i,1):sep(i,2),2), 'b','LineWidth',1), hold on
        legend('True','multi Ideal','multi Learned');
    end

    % data compare
    MSE_NS = 0;
    MSE_SM = 0;
    RMSE_NS = 0;
    RMSE_SM = 0;
    data_len = length(X(:, 1));
    for i = 1:data_len
        MSE_NS = MSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2));
        MSE_SM = MSE_SM + hypot(XM(i, 1) - Y(i, 1), XM(i, 2) - Y(i, 2));
        RMSE_NS = RMSE_NS + hypot(X(i, 1) - Y(i, 1), X(i, 2) - Y(i, 2))^2;
        RMSE_SM = RMSE_SM + hypot(XM(i, 1) - Y(i, 1), XM(i, 2) - Y(i, 2))^2;
    end
    MSE_NS = MSE_NS / data_len
    MSE_SM = MSE_SM / data_len
    RMSE_NS = RMSE_NS / data_len;
    RMSE_SM = RMSE_SM / data_len;
end

% % figure compare
% n = length(Y(1,:));
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tspan,Y(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(tspan,XS(:,i),'r--','Linewidth',[2])
%     hold on;
%     plot(tspan,X(:,i),'b--','Linewidth',[2])
%     hold on;
%     plot(tspan,XM(:,i),'r---','Linewidth',[2])
%     hold on;
%     plot(tspan,XO(:,i),'b---','Linewidth',[2])
%     hold on;
%     legend('True','Learned single','Ideal single','Learned multi','Ideal multi')
%     title('Reconstruction of Training Data')
% end

%% Forecasting on Validation Data

getValidData;

if dataGenerateMode == 1
    XS_v = Maaav_single_SINDy_discrep(ValidY, ValidU, ValidT, Xi, Xintercept, avp, polyorder, sepV, poolMode, fitMode);
    XO_v = Maaav_no_discrep(ValidY, ValidU, ValidT, avp, sepV);
end
XM_v = Maaav_SINDy_discrep(ValidY, ValidU, ValidT, Xi, Xintercept, avp, polyorder, sepV, poolMode, fitMode);

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XS_v(1:sepV(1,2),1), XS_v(1:sepV(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(1:sepV(1,2),1), XO_v(1:sepV(1,2),2), 'r','LineWidth',1), hold on
    plot(XM_v(1:sepV(1,2),1), XM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XS_v(sepV(2,1):sepV(2,2),1), XS_v(sepV(2,1):sepV(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_v(sepV(2,1):sepV(2,2),1), XO_v(sepV(2,1):sepV(2,2),2), 'r','LineWidth',1), hold on
    plot(XM_v(sepV(2,1):sepV(2,2),1), XM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);

    % third picture
%     subplot(1,3,3);
%     plot(ValidY(sepV(3,1):sepV(3,2),1), ValidY(sepV(3,1):sepV(3,2),2), 'k','LineWidth',1), hold on
%     plot(ValidX(sepV(3,1):sepV(3,2),1), ValidX(sepV(3,1):sepV(3,2),2), 'r--','LineWidth',1), hold on
%     plot(XS_v(sepV(3,1):sepV(3,2),1), XS_v(sepV(3,1):sepV(3,2),2), 'b--','LineWidth',1), hold on
%     plot(XO_v(sepV(3,1):sepV(3,2),1), XO_v(sepV(3,1):sepV(3,2),2), 'r','LineWidth',1), hold on
%     plot(XM_v(sepV(3,1):sepV(3,2),1), XM_v(sepV(3,1):sepV(3,2),2), 'b','LineWidth',1), hold on
%     legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     sgtitle('Trajectory of Systems of Valid environment')
%     axis([325 330 -305 -280]);
    
    % data compare
    MSE_NS_v = 0;
    MSE_NM_v = 0;
    MSE_SS_v = 0;
    MSE_SM_v = 0;
    RMSE_NS_v = 0;
    RMSE_NM_v = 0;
    RMSE_SS_v = 0;
    RMSE_SM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_NS_v = MSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_NM_v = MSE_NM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2));
        MSE_SS_v = MSE_SS_v + hypot(XS_v(i, 1) - ValidY(i, 1), XS_v(i, 2) - ValidY(i, 2));
        MSE_SM_v = MSE_SM_v + hypot(XM_v(i, 1) - ValidY(i, 1), XM_v(i, 2) - ValidY(i, 2));
        RMSE_NS_v = RMSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_NM_v = RMSE_NM_v + hypot(XO_v(i, 1) - ValidY(i, 1), XO_v(i, 2) - ValidY(i, 2))^2;
        RMSE_SS_v = RMSE_SS_v + hypot(XS_v(i, 1) - ValidY(i, 1), XS_v(i, 2) - ValidY(i, 2))^2;
        RMSE_SM_v = RMSE_SM_v + hypot(XM_v(i, 1) - ValidY(i, 1), XM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_NS_v = MSE_NS_v / data_len_v;
    MSE_NM_v = MSE_NM_v / data_len_v
    MSE_SS_v = MSE_SS_v / data_len_v;
    MSE_SM_v = MSE_SM_v / data_len_v
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_NM_v = RMSE_NM_v / data_len_v;
    RMSE_SS_v = RMSE_SS_v / data_len_v;
    RMSE_SM_v = RMSE_SM_v / data_len_v;
end

if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(ValidY(1:sepV(1,2),1), ValidY(1:sepV(1,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(1:sepV(1,2),1), ValidX(1:sepV(1,2),2), 'r--','LineWidth',1), hold on
    plot(XM_v(1:sepV(1,2),1), XM_v(1:sepV(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(ValidY(sepV(2,1):sepV(2,2),1), ValidY(sepV(2,1):sepV(2,2),2), 'k','LineWidth',1), hold on
    plot(ValidX(sepV(2,1):sepV(2,2),1), ValidX(sepV(2,1):sepV(2,2),2), 'r--','LineWidth',1), hold on
    plot(XM_v(sepV(2,1):sepV(2,2),1), XM_v(sepV(2,1):sepV(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
%     axis([358 374 -175 -135]);

    % third picture
%     subplot(1,3,3);
%     plot(ValidY(sepV(3,1):sepV(3,2),1), ValidY(sepV(3,1):sepV(3,2),2), 'k','LineWidth',1), hold on
%     plot(ValidX(sepV(3,1):sepV(3,2),1), ValidX(sepV(3,1):sepV(3,2),2), 'r--','LineWidth',1), hold on
%     plot(XM_v(sepV(3,1):sepV(3,2),1), XM_v(sepV(3,1):sepV(3,2),2), 'b','LineWidth',1), hold on
%     legend('True','multi Ideal','multi Learned');
%     sgtitle('Trajectory of Systems of Valid environment')
%     axis([325 330 -305 -280]);
    
    % data compare
    MSE_NS_v = 0;
    MSE_SM_v = 0;
    RMSE_NS_v = 0;
    RMSE_SM_v = 0;

    data_len_v = length(ValidX(:, 1));
    for i = 1:data_len_v
        MSE_NS_v = MSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2));
        MSE_SM_v = MSE_SM_v + hypot(XM_v(i, 1) - ValidY(i, 1), XM_v(i, 2) - ValidY(i, 2));
        RMSE_NS_v = RMSE_NS_v + hypot(ValidX(i, 1) - ValidY(i, 1), ValidX(i, 2) - ValidY(i, 2))^2;
        RMSE_SM_v = RMSE_SM_v + hypot(XM_v(i, 1) - ValidY(i, 1), XM_v(i, 2) - ValidY(i, 2))^2;
    end

    MSE_NS_v = MSE_NS_v / data_len_v
    MSE_SM_v = MSE_SM_v / data_len_v
    RMSE_NS_v = RMSE_NS_v / data_len_v;
    RMSE_SM_v = RMSE_SM_v / data_len_v;
end

%% Forecasting on Test Data

getTestData;

if dataGenerateMode == 1
    XS_t = Maaav_single_SINDy_discrep(TestY, TestU, TestT, Xi, Xintercept, avp, polyorder, sepT, poolMode, fitMode);
    XO_t = Maaav_no_discrep(TestY, TestU, TestT, avp, sepT);
end
XM_t = Maaav_SINDy_discrep(TestY, TestU, TestT, Xi, Xintercept, avp, polyorder, sepT, poolMode, fitMode);

if dataGenerateMode == 1
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XS_t(1:sepT(1,2),1), XS_t(1:sepT(1,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(1:sepT(1,2),1), XO_t(1:sepT(1,2),2), 'r','LineWidth',1), hold on
    plot(XM_t(1:sepT(1,2),1), XM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XS_t(sepT(2,1):sepT(2,2),1), XS_t(sepT(2,1):sepT(2,2),2), 'b--','LineWidth',1), hold on
    plot(XO_t(sepT(2,1):sepT(2,2),1), XO_t(sepT(2,1):sepT(2,2),2), 'r','LineWidth',1), hold on
    plot(XM_t(sepT(2,1):sepT(2,2),1), XM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);

    % third picture
%     subplot(1,3,3);
%     plot(TestY(sepT(3,1):sepT(3,2),1), TestY(sepT(3,1):sepT(3,2),2), 'k','LineWidth',1), hold on
%     plot(TestX(sepT(3,1):sepT(3,2),1), TestX(sepT(3,1):sepT(3,2),2), 'r--','LineWidth',1), hold on
%     plot(XS_t(sepT(3,1):sepT(3,2),1), XS_t(sepT(3,1):sepT(3,2),2), 'b--','LineWidth',1), hold on
%     plot(XO_t(sepT(3,1):sepT(3,2),1), XO_t(sepT(3,1):sepT(3,2),2), 'r','LineWidth',1), hold on
%     plot(XM_t(sepT(3,1):sepT(3,2),1), XM_t(sepT(3,1):sepT(3,2),2), 'b','LineWidth',1), hold on
%     legend('True','single Ideal', 'single Learned', 'multi Ideal','multi Learned');
%     sgtitle('Trajectory of Systems of Test environment')
    % axis([325 330 -305 -280]);
    
    % data compare
    MSE_NS_t = 0;
    MSE_NM_t = 0;
    MSE_SS_t = 0;
    MSE_SM_t = 0;
    RMSE_NS_t = 0;
    RMSE_NM_t = 0;
    RMSE_SS_t = 0;
    RMSE_SM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_NS_t = MSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_NM_t = MSE_NM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2));
        MSE_SS_t = MSE_SS_t + hypot(XS_t(i, 1) - TestY(i, 1), XS_t(i, 2) - TestY(i, 2));
        MSE_SM_t = MSE_SM_t + hypot(XM_t(i, 1) - TestY(i, 1), XM_t(i, 2) - TestY(i, 2));
        RMSE_NS_t = RMSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_NM_t = RMSE_NM_t + hypot(XO_t(i, 1) - TestY(i, 1), XO_t(i, 2) - TestY(i, 2))^2;
        RMSE_SS_t = RMSE_SS_t + hypot(XS_t(i, 1) - TestY(i, 1), XS_t(i, 2) - TestY(i, 2))^2;
        RMSE_SM_t = RMSE_SM_t + hypot(XM_t(i, 1) - TestY(i, 1), XM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_NS_t = MSE_NS_t / data_len_t;
    MSE_NM_t = MSE_NM_t / data_len_t
    MSE_SS_t = MSE_SS_t / data_len_t;
    MSE_SM_t = MSE_SM_t / data_len_t
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_NM_t = RMSE_NM_t / data_len_t;
    RMSE_SS_t = RMSE_SS_t / data_len_t;
    RMSE_SM_t = RMSE_SM_t / data_len_t;
end
    
if dataGenerateMode == 2
    figure,
    % first picture
    subplot(1,2,1);
    plot(TestY(1:sepT(1,2),1), TestY(1:sepT(1,2),2), 'k','LineWidth',1), hold on
    plot(TestX(1:sepT(1,2),1), TestX(1:sepT(1,2),2), 'r--','LineWidth',1), hold on
    plot(XM_t(1:sepT(1,2),1), XM_t(1:sepT(1,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([360 365 -177 -162]);

    % second picture
    subplot(1,2,2);
    plot(TestY(sepT(2,1):sepT(2,2),1), TestY(sepT(2,1):sepT(2,2),2), 'k','LineWidth',1), hold on
    plot(TestX(sepT(2,1):sepT(2,2),1), TestX(sepT(2,1):sepT(2,2),2), 'r--','LineWidth',1), hold on
    plot(XM_t(sepT(2,1):sepT(2,2),1), XM_t(sepT(2,1):sepT(2,2),2), 'b','LineWidth',1), hold on
    legend('True','multi Ideal','multi Learned');
    % axis([358 374 -175 -135]);

    % third picture
%     subplot(1,3,3);
%     plot(TestY(sepT(3,1):sepT(3,2),1), TestY(sepT(3,1):sepT(3,2),2), 'k','LineWidth',1), hold on
%     plot(TestX(sepT(3,1):sepT(3,2),1), TestX(sepT(3,1):sepT(3,2),2), 'r--','LineWidth',1), hold on
%     plot(XM_t(sepT(3,1):sepT(3,2),1), XM_t(sepT(3,1):sepT(3,2),2), 'b','LineWidth',1), hold on
%     legend('True','multi Ideal','multi Learned');
%     sgtitle('Trajectory of Systems of Test environment')
    % axis([325 330 -305 -280]);
    
    % data compare
    MSE_NS_t = 0;
    MSE_SM_t = 0;
    RMSE_NS_t = 0;
    RMSE_SM_t = 0;

    data_len_t = length(TestX(:, 1));
    for i = 1:data_len_t
        MSE_NS_t = MSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2));
        MSE_SM_t = MSE_SM_t + hypot(XM_t(i, 1) - TestY(i, 1), XM_t(i, 2) - TestY(i, 2));
        RMSE_NS_t = RMSE_NS_t + hypot(TestX(i, 1) - TestY(i, 1), TestX(i, 2) - TestY(i, 2))^2;
        RMSE_SM_t = RMSE_SM_t + hypot(XM_t(i, 1) - TestY(i, 1), XM_t(i, 2) - TestY(i, 2))^2;
    end

    MSE_NS_t = MSE_NS_t / data_len_t
    MSE_SM_t = MSE_SM_t / data_len_t
    RMSE_NS_t = RMSE_NS_t / data_len_t;
    RMSE_SM_t = RMSE_SM_t / data_len_t;
end
