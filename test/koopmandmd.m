
% ����ϵͳ����
A = [0.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];
x0 = [1; 1];
u = 1;
m = 10;  % ʱ�������

% ����״̬���� X �� X_prime
X = zeros(2, m);
X_prime = zeros(2, m);
U = u * ones(1, m);
X(:, 1) = x0;

for k = 1:m-1
    X_prime(:, k) = A * X(:, k) + B * U(:, k);
    X(:, k+1) = X_prime(:, k);
end
X_prime(:, m) = A * X(:, m) + B * U(:, m);

% ʹ�ö���ʽ������������չ
phi = @(x) [x(1); x(2); x(1)^2; x(2)^2; x(1)*x(2)];
Psi_X = zeros(5, m);
Psi_X_prime = zeros(5, m);

for k = 1:m
    Psi_X(:, k) = phi(X(:, k));
    Psi_X_prime(:, k) = phi(X_prime(:, k));
end

%%
% % Koopman DMD ����
% [Phi, Lambda] = koopmanDMD(Psi_X, Psi_X_prime);

% �����С
[q, m] = size(Psi_X);

% ���� K ����
K = Psi_X_prime * pinv(Psi_X);

% �����ֽ� K
[W, Lambda] = eig(K);

% ���� Koopman DMD ģ̬
Phi = X_prime * W * inv(Lambda);
      
%%
% �ؽ�ѵ������
X_reconstructed = reconstructData(Phi, Lambda, Psi_X(:,1), m);

% ��ӡ���
disp('Koopman DMD Modes (Phi):');
disp(Phi);
disp('Eigenvalues (Lambda):');
disp(diag(Lambda));
disp('Reconstructed Training Data:');
disp(X_reconstructed);

% ���Լ�������֤
n_test = 5;  % ���Լ�ʱ�������
X_test = zeros(2, n_test);
X_test(:, 1) = X(:, end);  % ���Լ���ʼ״̬Ϊѵ����ĩ״̬

for k = 1:n_test-1
    X_test(:, k+1) = A * X_test(:, k) + B * u;
end

% ��չ��������
Psi_X_test = zeros(5, n_test);
for k = 1:n_test
    Psi_X_test(:, k) = phi(X_test(:, k));
end

% Ԥ����Լ�����
X_test_predicted = reconstructData(Phi, Lambda, Psi_X_test(:,1), n_test);

% �������
error = norm(X_test - X_test_predicted(1:2, :), 'fro');
disp('Test Data:');
disp(X_test);
disp('Predicted Test Data:');
disp(X_test_predicted(1:2, :));
disp(['Reconstruction Error: ', num2str(error)]);

function X_reconstructed = reconstructData(Phi, Lambda, Psi_x0, num_steps)
    % ʹ�� DMD ģ̬������ֵ�ؽ�����
    % ����:
    % Phi         - DMD ģ̬����
    % Lambda      - DMD ����ֵ����
    % Psi_x0      - ��ʼ״̬����չ�ռ�ı�ʾ
    % num_steps   - ʱ�䲽����
    % ���:
    % X_reconstructed - �ؽ������ݾ���

    q = size(Phi, 1);  % ��չ�ռ�ά��
    X_reconstructed = zeros(q, num_steps);
    X_reconstructed(:, 1) = Psi_x0;

    for k = 2:num_steps
        X_reconstructed(:, k) = Phi * (Lambda^(k-1)) * (Phi \ Psi_x0);
    end
end
