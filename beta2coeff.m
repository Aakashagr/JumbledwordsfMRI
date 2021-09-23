% converting betas from non-linear fit to coefficients of ISI model
function coeff =  beta2coeff(beta)
global strlen
polyB = [beta(1) beta(2) beta(3)];
W1 = []; W2 = [];
for i = 1:strlen
    for j = 1:strlen - i+1; W1 = [W1; i^2 i 1]; end
    for j = 1:strlen - i;   W2 = [W2; i^2 i 1]; end
end
weightM = [W1;W2]';

dc1 = []; dc2 = [];
for i = 1:strlen - 1
    dc2 = [dc2; 0];
    for j = i:-1:1
        dc1 = [dc1; j];
        dc2 = [dc2; j];
    end
end
decayC = [0; flipud(dc2); flipud(dc1)];

coeff = [];
if beta(4) > 10^2; beta(4) = -beta(4); end % Avoiding nlinfit error due to non-convergence
for i = 1:nchoosek(strlen,2) + strlen
    coeff(end+1) = polyB*weightM(:,i)*exp(decayC(i)*beta(4));
end

for j = i+1: i + nchoosek(strlen,2)
    coeff(end+1) =  -polyB*weightM(:,j)*exp(decayC(j)*beta(4));
end

% coeff(end+1) = beta(5);
