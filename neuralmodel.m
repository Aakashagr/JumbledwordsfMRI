
function dpred = neuralmodel(w,xinput)
allparts = xinput{1}; rates = xinput{2}; 
ncells = size(rates,2); 
npairs = size(allparts,1); nparts = size(allparts,2)/2; 
wmat = reshape(w(1:end-1),[],ncells); % different spatial weighting for each cell 
% wmat = repmat(w(1:end-1),1,nparts); % common spatial weighting across all cells
rates(63,:) = zeros(1,ncells); % adding an extra cell that has zero response (for unequal length strings)

% calculate distance between each pair of ngrams 
for pairid = 1:npairs
    p1 = allparts(pairid,1:nparts); 
    p2 = allparts(pairid,nparts+1:end); 
    r1 = sum(rates(p1,:).*wmat,1); % r1 should be ncells x 1 response vector for the first stim
    r2 = sum(rates(p2,:).*wmat,1); % r2 should be ncells x 1 response vector for the first stim
    dpred(pairid,1) = sum(abs(r1-r2)); % L1 norm of r1-r2
%     dpred(pairid,1) = sum((r1-r2).^2); % L2 norm of r1-rd2
end
dpred = dpred + w(end);% + sum(w.^2); % adding a constant term

return