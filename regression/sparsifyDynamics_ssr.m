function Xi = sparsifyDynamics_ssr(Theta, dXdt)

% Compute sparse regression of a discrepancy (ONE STATE): iterative least squares
% SINDy with Stepwise Sparse Regressor
% Sparse learning of stochastic dynamical equations. 

%% Peform sparse regression

rightDXDT = dXdt(:, 6:10);

% initial guess: Least-squares
Xi = Theta \ rightDXDT;  
[n, m] = size(Theta);

% get rid of zero terms in Xi and repsectively Theta and dXdt
if any(Xi == 1)
    Xi_hat = find(Xi == 0);
    Theta(:,Xi_hat) = [];    
    Xi = Theta\rightDXDT;  % initial guess: Least-squares
    
    M = m - length(Xi_hat);
else
    % if no coeffs with zeros OR once removed, find min(abs(Xi)) = 0 and re-gregress
    M = m;
end

for i = 1:M-1
    Xi(abs(Xi) == min(abs(Xi(:)))) = 0;
    % [row, col] = find(Xi == 0);
    % Theta(:,row) = [];
    Xi = Theta \ rightDXDT;  % initial guess: Least-squares
    
end
