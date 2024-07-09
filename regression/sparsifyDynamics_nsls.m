function Xi = sparsifyDynamics_nsls(Theta, dXdt, lambda)

% compute Sparse regression: sequential least squares
% SINDy with library normalization

for norm_k = 1:size(Theta,2)
    normLib(norm_k) = norm(Theta(:,norm_k));
    Theta(:,norm_k) = Theta(:,norm_k)/normLib(norm_k);
end

rightDXDT = dXdt(:, 6:10);
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

for norm_k = 1:length(Xi)
    Xi(norm_k,:) = Xi(norm_k,:)/normLib(norm_k);
end
