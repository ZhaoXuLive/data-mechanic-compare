function Phi_optimized = optimizeModes(Phi, Lambda, b, X, m, itern, learnv)
    % �Ż�DMDģ̬�ĺ�����������С���˷��Ż�
    % ��������ֻ��չʾһ���򵥵��Ż�ʾ��
    % ʵ��Ӧ���п���ʹ�ø����ӵ��Ż�����

    % ��ʼ���Ż����Phi
    Phi_optimized = Phi;

    % �Ż�����
    for i = 1:itern % �Ż���������
        % �����ؽ����
        X_reconstructed = zeros(size(X));
        for t = 1:m
            X_reconstructed(:,t) = Phi_optimized * (Lambda^(t-1)) * b;
        end
        error = X - X_reconstructed;
        
        avgError = zeros(10, 1);
        for j = 1:10
            avgError(j) = mean(error(j, :));
        end
        
        % ����Phi_optimized
        Phi_optimized = Phi_optimized + learnv * avgError * b'; 
    end
end