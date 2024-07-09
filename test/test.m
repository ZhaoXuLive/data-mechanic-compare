% ��չDMDʾ������

% ����ʾ������
t = 0:1:9; % ʱ��������10��ʱ���
n = length(t); % ���ݳ���

% ϵͳ������������
A = [1.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];

% ���������ź�
u = sin(0.1 * t); % �����źţ����Ը�����Ҫ����

% ��ʼ��״̬����
x0 = [1; 1]; % ��ʼ״̬
X = zeros(2, n); % ״̬����
X(:, 1) = x0;

% ����״̬����
for i = 2:n
    X(:, i) = A * X(:, i-1) + B * u(i-1);
end

% ������չ״̬���󣬰��������ź�
X_ext = [X(:, 1:end-1); u(1:end-1)];
X_prime_ext = X(:, 2:end);

% DMD����
[U, S, V] = svd(X_ext, 'econ');
Atilde = U' * X_prime_ext * V * diag(1./diag(S));
[W, D] = eig(Atilde);
Phi = X_prime_ext * V * diag(1./diag(S)) * W;

% ��̬ģʽ
lambda = diag(D);
omega = log(lambda); % Ƶ��

% ��ʼ״̬
b = Phi \ X_prime_ext(:, 1); % ��ʼ״̬

% ʱ���ݻ�
time_dynamics = zeros(length(b), n - 1);
for i = 1:n-1
    time_dynamics(:,i) = b .* exp(omega * t(i));
end

% �ع��ź�
x_reconstructed = Phi * time_dynamics;

% ��ͼ
figure;
plot(t, X(1, :), 'b', 'DisplayName', 'ԭʼ�ź�');
hold on;
plot(t(2:end), real(x_reconstructed(1, :)), 'r--', 'DisplayName', '�ع��ź�');
legend;
title('�ź��ع�');
xlabel('ʱ��');
ylabel('״̬����');

% ��֤���ϵ�Ԥ��
t_pred = 10:1:14; % �µ�ʱ������
u_pred = sin(0.1 * t_pred); % ��֤�������ź�
X_pred = zeros(2, length(t_pred));
X_pred(:, 1) = X(:, end);

for i = 2:length(t_pred)
    X_pred(:, i) = A * X_pred(:, i-1) + B * u_pred(i-1);
end

X_pred_ext = [X_pred(:, 1:end-1); u_pred(1:end-1)];

% �ع���֤���ź�
time_dynamics_pred = zeros(length(b), length(t_pred)-1);
for i = 1:length(t_pred)-1
    time_dynamics_pred(:,i) = b .* exp(omega * (t_pred(i) - t(end)));
end

x_reconstructed_pred = Phi * time_dynamics_pred;

% ��ͼ
figure;
plot(t_pred, X_pred(1, :), 'b', 'DisplayName', 'ԭʼ�ź�');
hold on;
plot(t_pred(2:end), real(x_reconstructed_pred(1, :)), 'r--', 'DisplayName', '�ع��ź�');
legend;
title('��֤���ź��ع�');
xlabel('ʱ��');
ylabel('״̬����');
