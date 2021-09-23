function ypred =  batonfun(b,x)
global strlen
polyB = [b(1) b(2) b(3)];
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

ypred = zeros(size(x,1),1);
if b(4) > 10^2; b(4) = -b(4); end % Avoiding nlinfit error due to non-convergence
for i = 1:nchoosek(strlen,2) + strlen
    ypred = ypred + polyB*weightM(:,i)*exp(decayC(i)*b(4))*x(:,i);
end

for j = i+1: i + nchoosek(strlen,2)
    ypred = ypred - polyB*weightM(:,j)*exp(decayC(j)*b(4))*x(:,j);
end

ypred = ypred + b(5)*ones(size(x,1),1);
