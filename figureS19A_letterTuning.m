allclearL2
if ~exist('L2_str'); load L2fmri_LDTw; end
BETAs = L2_str.mergedevtbeta;
qw = 1:32; qn = 33:64; ql = 65:74;
[idx, ROIname] = getvoxind(L2_str);

%% load tuning curve from perceptual dissimilarity
load dis_letter; sl = 'ESAROLITND'-64; slpair = sl(nchoosek(1:10,2)); % Single letter
for i = 1:size(slpair,1); dpsy(i,1) = dis_uletter(find(ismember(nchoosek(1:26,2),[slpair(i,:); fliplr(slpair(i,:))],'rows')));end
D = squareform(dpsy);
opts = statset('MaxIter',3000);
psy_lt = mdscale(D,6,'Options',opts);
dpred = pdist(psy_lt); [c,p]=nancorrcoef(dpsy,dpred);

% Normalizing the data
psy_lt = zscore(psy_lt);
%% single letter tuning curves
for roi = 1:5
    actl = []; acts = [];
    for sub = 1:numel(L2_str.subjectid)
        betas = L2_str.mergedevtbeta{sub};
        if roi == 4; max_vox = 40; else, max_vox = Inf; end
        nvox = min(numel(idx{sub,roi}),max_vox);
        actl = [actl; betas(idx{sub,roi}(1:nvox),ql)];  % mean single letter responses
    end
    actl(isnan(sum(actl,2)),:) = [];% actl = zscore(actl);
    actl = zscore(actl');
    
    % Finding the best MDS dimension for each voxel
    mdsmatch = [];
    for mdim = 1:size(psy_lt,2)
        mdsmatch(mdim,:) = mean((actl - psy_lt(:,mdim)).^2);
    end
    
    % Selecting the best dimension
    [errval,id] = sort(mdsmatch);
    BestDim = id(1,:); actlSorted = [];
    for mdim = 1:size(psy_lt,2)
        qt = find(BestDim == mdim);
        errtemp = errval(qt);
        [~,id2] = sort(errtemp);
        qt = qt(id2);
        actlSorted = [actlSorted; actl(:,qt)'];
        voxids{mdim,roi} = qt;
        MactlSorted(:,mdim,roi) = nanmean(actl(:,qt),2);
    end    
    subplot(2,3,roi); imagesc(actlSorted); colormap(flipud(pink)); title(ROIname{roi});
    totalvox(roi,1) = size(id,2);
end

for i = 1:5, r(i) = nancorrcoef(psy_lt,MactlSorted(:,:,i)); end

%%
figure;
imagesc(zscore(psy_lt)'); colormap(flipud(pink))
% for i = 1:size(psy_lt,2)
%     subplot(2,3,i); plot(psy_lt(:,i));
% end
% suptitle('Comparing tuning curve')
% 
