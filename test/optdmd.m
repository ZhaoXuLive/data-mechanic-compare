
% ϵͳ����
A = [1.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];

% ��ʼ״̬
x0 = [1; 1];

% ����
u = 1;

% ʱ�������
m = 10;

% ����״̬���� X �� X_prime
X = zeros(2, m);
X_prime = zeros(2, m);
U = u * ones(1, m);

% ��ʼ״̬��ֵ
X(:, 1) = x0;

% ��������
for k = 1:m-1
    X_prime(:, k) = A * X(:, k) + B * U(:, k);
    X(:, k+1) = X_prime(:, k);
end
X_prime(:, m) = A * X(:, m) + B * U(:, m);

% �����Ż�DMD����
optimizedDMDwithInputs(X, X_prime, U);

function optimizedDMDwithInputs(X, X_prime, U)
    % ����:
    % X       - ״̬���󣬴�СΪ n x m
    % X_prime - ��һ��ʱ�䲽��״̬���󣬴�СΪ n x m
    % U       - ������󣬴�СΪ p x m

    % �����С
    [n, m] = size(X);
    p = size(U, 1);

    % ������չ���� Omega
    Omega = [X; U];

    % ���� Omega ��α��
    Omega_pseudo_inverse = pinv(Omega);

    % ���� K ����
    K = X_prime * Omega_pseudo_inverse;

    % �� K �зֽ�� A �� B
    A = K(:, 1:n);
    B = K(:, n+1:end);

    % ����ֵ�ֽ� X
    [U_X, Sigma_X, V_X] = svd(X, 'econ');

    % ���㽵ά��� A_tilde
    A_tilde = U_X' * A * U_X;

    % �����ֽ� A_tilde
    [W, Lambda] = eig(A_tilde);

    % ���� DMD ģ̬
    Phi = X_prime * V_X / Sigma_X * W;

    % ��ʼ������DMD����
    b = Phi \ X(:,1);

    % �Ż�DMD������������С���˷��Ż�Ϊ����
    Phi_optimized = optimizeModes(Phi, Lambda, b, X, m);

    % �ؽ�ѵ��������
    X_dmd = zeros(n, m);
    for i = 1:m
        X_dmd(:,i) = Phi_optimized * (Lambda^(i-1)) * b;
    end

    % ������
    disp('A matrix:');
    disp(A);
    disp('B matrix:');
    disp(B);
    
    t = 1:1:10;
    figure,
    plot(t, X(1, :));
    hold on;
    plot(t, X_dmd(1, :), 'r');
    figure,
    plot(t, X(2, :));
    hold on;
    plot(t, X_dmd(2, :), 'r');

end

function Phi_optimized = optimizeModes(Phi, Lambda, b, X, m)
    % �Ż�DMDģ̬�ĺ�����������С���˷��Ż�
    % ��������ֻ��չʾһ���򵥵��Ż�ʾ��
    % ʵ��Ӧ���п���ʹ�ø����ӵ��Ż�����

    % ��ʼ���Ż����Phi
    Phi_optimized = Phi;

    % �Ż�����
    for i = 1:10 % �Ż���������
        % �����ؽ����
        X_reconstructed = zeros(size(X));
        for t = 1:m
            X_reconstructed(:,t) = Phi_optimized * (Lambda^(t-1)) * b;
        end
        error = X - X_reconstructed;
        
        avgError = zeros(2, 1);
        avgError(1) = mean(error(1, :));
        avgError(2) = mean(error(2, :));
        
        % ����Phi_optimized
        Phi_optimized = Phi_optimized + 0.01 * avgError * b'; % ѧϰ��Ϊ0.01
    end
end
