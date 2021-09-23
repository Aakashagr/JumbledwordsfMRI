%% Fitting neural model string searches
if ~exist('L2_str'); load L2fmri_LDTw; end

%%
allclearL2
qw = 1:32; qn = 33:64; ql = 65:74;
stimid = L2_str.stimid;

%% load tuning curve from perceptual dissimilarity
load dis_letter; sl = 'ESAROLITND'-64; slpair = sl(nchoosek(1:10,2)); % Single letter
for i = 1:numel(sl); stimid(stimid == sl(i)) = 100+i; end; stimid = stimid - 100;

for i = 1:size(slpair,1); dpsy(i,1) = dis_uletter(find(ismember(nchoosek(1:26,2),[slpair(i,:); fliplr(slpair(i,:))],'rows')));end
D = squareform(dpsy);
opts = statset('MaxIter',3000); ncells = 6;
[rates,e] = mdscale(D,ncells,'Options',opts);

%% 
imgparts = stimid; 
srchpairs = [nchoosek([1:32],2); nchoosek([1:32],2)+32; [1:32;33:64]']; 
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
X = {allparts,rates}; 

dobs = [L2_str.search.w; L2_str.search.n; L2_str.search.wnw]; 
nparts = 5; 

opts = statset('nlinfit'); % opts.RobustWgtFun = 'bisquare'; opts.MaxIter = 1000;
% rng('default')
w0 = rand(ncells*nparts+1,1); 
[west,r,J] = nlinfit(X,dobs,@neuralmodel,w0,opts);
dpred = neuralmodel(west,X); 
wmat = reshape(west(1:end-1),[],ncells); % different spatial weighting for each cell 
%%
figure
% for neu = 1:size(wmat,2)
%     subplot(2,3,neu); plot(wmat(:,neu),'LineWidth',2); title(['Neuron -' num2str(neu)]);
%     ylabel('Coefficient'); xlabel('Letter position')
% end
imagesc(zscore(wmat)'); colormap(flipud(pink))

%% Voxel data
[idx, ROIname] = getvoxind(L2_str);
figure; 
for roi = 1:5
    Bw = []; cnt = 1;
    for sub = 1:numel(L2_str.subjectid)
        betas = L2_str.mergedevtbeta{sub};
        if roi == 4; max_vox = 40; else, max_vox = Inf; end
        
        nvox = min(numel(idx{sub,roi}),max_vox);
        Betas = betas(idx{sub,roi}(1:nvox),:);
        Betas(isnan(sum(Betas,2)),:) = [];
        
        acts  = Betas(:,[qw qn]);
        actl  = Betas(:,ql);
        
        for vox = 1:size(acts,1)
            yobs = (acts(vox,:)'); yl = actl(vox,:)';
            Xmat = [(yl(stimid(:,:))) ones(numel(yobs),1)];  w = regress(yobs, Xmat);
            Bw(:,cnt) = w(1:end-1); cnt = cnt+1;
        end
    end
    
    Bw = zscore(Bw); wmat = zscore(wmat);
    
    % Finding the best weight match for each voxel
    Wmatch = [];
    for mdim = 1:size(wmat,2)
        Wmatch(mdim,:) = mean((Bw - wmat(:,mdim)).^2);
    end
    
    % Selecting the best dimension
    [errval,id] = sort(Wmatch);
    BestDim = id(1,:); BwSorted = [];
    for mdim = 1:size(wmat,2)
        qt = find(BestDim == mdim);
        errtemp = errval(qt);
        [~,id2] = sort(errtemp);
        qt = qt(id2);
        BwSorted = [BwSorted; Bw(:,qt)'];
    end
    figure(2); subplot(2,3,roi); imagesc(BwSorted); colormap(flipud(pink)); title(ROIname{roi});
    totalvox(roi,1) = size(id,2);
end
