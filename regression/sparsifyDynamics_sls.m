function Xi = sparsifyDynamics_sls(Theta, dXdt, lambda)

% compute Sparse regression: sequential least squares

% Xi = Theta\dXdt;  % initial guess: Least-squares
% 
% % lambda is our sparsification knob.
% for k=1:10
%     smallinds = (abs(Xi)<lambda);   % find small coefficients
%     Xi(smallinds)=0;                % and threshold
%     for ind = 1:n                   % n is state dimension
%         biginds = ~smallinds(:,ind);
%         % Regress dynamics onto remaining terms to find sparse Xi
%         Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind); 
%     end
% end

lassoLambda = [9.6128e-08 2.5671e-07 9.7338e-09 8.6769e-09 1.1747e-08];
rightDXDT = dXdt(:, 6:10);

% Xi = Theta\rightDXDT;
% for ind = 1:5
%     for k=1:10
%         smallinds = (abs(Xi)<lassoLambda(ind));     % find small coefficients
%         Xi(smallinds) = 0;                % and threshold
%         biginds = ~smallinds(:,ind);
%         Xi(biginds, ind) = Theta(:, biginds)\rightDXDT(:, ind);
% %         Xi(biginds,ind) = lsqminnorm(Theta(:,biginds), rightDXDT(:,ind)); 
%     end
% end

Xi = Theta\rightDXDT;
% Xi = lsqminnorm(Theta, rightDXDT);
for k=1:10
    smallinds = (abs(Xi)<lambda);     % find small coefficients
    Xi(smallinds) = 0;                % and threshold
    for ind = 1:5                     % n is state dimension
        biginds = ~smallinds(:,ind);
        Xi(biginds,ind) = Theta(:,biginds)\rightDXDT(:,ind);
%         Xi(biginds,ind) = lsqminnorm(Theta(:,biginds), rightDXDT(:,ind)); 
    end
end
