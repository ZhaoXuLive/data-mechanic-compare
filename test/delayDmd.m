% �ӳ�DMD���������Ķ�̬ϵͳ

% ��������
n = 2;  % ״̬ά��
m = 1;  % ����ά��
M = 10; % ���ݳ���
d = 3; % �ӳ�Ƕ��ά��

% ����״̬����������
A = [1.9, 0.1; -0.2, 0.8];
B = [0.1; 0.05];
x = zeros(n, M);
% ����
u = ones(1, M);

x(:,1) = [1; 1];

for k = 1:M-1
    x(:,k+1) = A*x(:,k) + B*u(k);
end

% �����ӳ�Ƕ�����
Z = [];
Z_prime = [];
for k = 1:M-d
    Z = [Z, reshape(x(:,k:k+d-1), [], 1)];
    Z_prime = [Z_prime, reshape(x(:,k+1:k+d), [], 1)];
end

% ������չ���������
U = [];
for k = 1:M-d
    U = [U, reshape(u(k:k+d-1), [], 1)];
end

W = [Z; U];

% ��⶯̬���� K
K = Z_prime * pinv(W);

% ������� A_d �� B_d
A_d = K(1:n*d, 1:n*d);
B_d = K(1:n*d, n*d+1:end);

% ģ̬����
[Phi, Lambda] = eig(A_d);

% ��ʼ����
Z0 = Z(:,1);
b = Phi \ Z0;

% �ؽ�״̬
Z_reconstructed = zeros(size(Z));
for k = 0:M-d-1
    Z_reconstructed(:,k+1) = Phi * (Lambda^k) * b;
end

% ���ƽ��
figure;
subplot(2,1,1);
plot(1:M, x(1,:), 'b', 'DisplayName', '��ʵ����');
hold on;
plot(d+1:M, Z_reconstructed(1,:), 'r--', 'DisplayName', '�ؽ�����');
xlabel('ʱ�䲽');
ylabel('״̬x_1');
legend;

subplot(2,1,2);
plot(1:M, x(2,:), 'b', 'DisplayName', '��ʵ����');
hold on;
plot(d+1:M, Z_reconstructed(2,:), 'r--', 'DisplayName', '�ؽ�����');
xlabel('ʱ�䲽');
ylabel('״̬x_2');
legend;

sgtitle('������Ķ�̬ϵͳ�ӳ�DMD�ؽ�');

% �ӳ�DMD������δ��������Ķ�̬ϵͳ
