if ~exist('L2_str'); load L2fmri_LDTw; end
allclearL2

dpsyw = L2_str.search.w;
dpsyn = L2_str.search.n;
Dw = squareform(dpsyw); [Yw] = mdscale(Dw,2);
Dn = squareform(dpsyn); [Yn] = mdscale(Dn,2);
[~,xYn] = procrustes(Yw,Yn,'scaling',0);

% Figure S15A
figure; plot(Yw(:,1),Yw(:,2),'.')
for i = 1:size(Yw,1)
    text(Yw(i,1),Yw(i,2),char(L2_str.stimid(i,:)+64));
end

% Figure S15B
figure; plot(xYn(:,1),xYn(:,2),'.')
for i = 1:size(xYn,1)
    text(xYn(i,1),xYn(i,2),char(L2_str.stimid(i+32,:)+64));
end


%% Figure S15C

load semanticfeat; dsemfeat = pdist(semanticfeat,'correlation')';
D = squareform(dsemfeat); [Y] = mdscale(D,2);
figure; plot(Y(:,1),Y(:,2),'.')
for i = 1:size(Y,1)
    text(Y(i,1),Y(i,2),char(L2_str.stimid(i,:)+64));
end

%% Estimating goodness of fit
nancorrcoef(pdist(Yw),dpsyw);
nancorrcoef(pdist(Yn),dpsyn);
nancorrcoef(pdist(Y),dsemfeat);