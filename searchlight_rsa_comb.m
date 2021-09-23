% Searchlight on dissimilarity using RSA toolbox
% allclear; load L2fmri_LDTw
allclearL2

qw = 1:32; qn = 33:64; ql = 65:74;
nsubs = numel(L2_str.subjectid);

%%%% Visual dissimilarity
load dis_letter; sl = 'ESAROLITND'-64; slpair = sl(nchoosek(1:10,2)); % Single letter
for i = 1:size(slpair,1); dpsyl(i,1) = dis_uletter(find(ismember(nchoosek(1:26,2),[slpair(i,:); fliplr(slpair(i,:))],'rows')));end

dpsyw = L2_str.search.w; dpsyn = L2_str.search.n; % Strings
dpsynw = L2_str.search.wnw; % Dissimilarity between words and nonwords
load semanticfeat; dsemfeat = pdist(semanticfeat,'correlation')'; % Semantic distance

%%%% Behaviour response
RT = L2_str.RT;
correctresp(1:32,1:2:nsubs) = 55; correctresp(1:32,2:2:nsubs) = 41;
correctresp(33:64,1:2:nsubs) = 41; correctresp(33:64,2:2:nsubs) = 55;
CR = repmat([correctresp; zeros(10,nsubs)],[1 1 8]);
iscorrect = L2_str.responsekey == CR; % Finding the correct response key
RT(iscorrect == 0) = nan; RT(RT == 0)= nan;  % excluding RTs with incorreact or no response
meanRTn = nanmean(nanmean(RT(qn,:,1:8),3),2);
meanRTw = nanmean(nanmean(RT(qw,:,1:8),3),2);

%%%% Initialsing weights
Rbeha = single(NaN(53,63,52));  Pbeha = Rbeha;  % dissimilarity between words with estimated dissimilarity
Rbehs = Rbeha;  Pbehs = Rbeha;                  % dissimilarity between word semantics with estimated dissimilarity
Rvmod  = single(NaN(53,63,52,nsubs));
ractRT = single(NaN(53,63,52,nsubs));  pactRT = ractRT;


wid = find(ismember(nchoosek(1:74,2),nchoosek(1:32,2),'rows'));
nid = find(ismember(nchoosek(1:74,2),nchoosek(1:32,2)+32,'rows'));
nwid = find(ismember(nchoosek(1:74,2),[1:32; 33:64]','rows'));
load median_droi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stimid = L2_str.stimid;
sl = 'ESAROLITND'-64; for i = 1:numel(sl); stimid(stimid == sl(i)) = 100+i; end; stimid = stimid - 100;

% Correlating mean dissimilarity with estimated behaviour dissimilarities
[x,y,z] = ind2sub([53 63 52],find(mask >0));

for sub = 1:nsubs
    beta = reshape(L2_str.mergedevtbeta{sub},[53 63 52 75]);  beta(:,:,:,75) = [];
    mask = reshape(L2_str.mergedevtbetash{sub},[53 63 52]);   mask(~isnan(mask)) = 1;
    [x,y,z] = ind2sub([53 63 52],find(mask >0));
    
    for idx = 1:numel(x)
        % Correlation between mean activity and response time
        [ractRT(x(idx),y(idx),z(idx),sub), pactRT(x(idx),y(idx),z(idx),sub)] = nancorrcoef(squeeze(beta(x(idx),y(idx),z(idx),[qw qn])), [meanRTw; meanRTn]);
        
        % Single voxel model
        robsm = squeeze(beta(x(idx),y(idx),z(idx),:));
        robs2  = robsm([qw]); r1 = robsm(ql); robs1 = r1(stimid([qw],:));
        X = [robs1 ones(length(robs2),1)]; b = regress(robs2,X); rpred2 = X*b;
        Rvmod1 = nancorrcoef(robs2,rpred2);

        robsm = squeeze(beta(x(idx),y(idx),z(idx),:));
        robs2  = robsm([qn]); r1 = robsm(ql); robs1 = r1(stimid([qn],:));
        X = [robs1 ones(length(robs2),1)]; b = regress(robs2,X); rpred2 = X*b;
        Rvmod2 = nancorrcoef(robs2,rpred2);

        Rvmod(x(idx),y(idx),z(idx),sub) = Rvmod1 - Rvmod2;
        
    end
    disp(sub);
end


for i = 1:17; mergedact(:,:,i) = L2_str.mergedevtbeta{i}; end; mergedact = reshape(nanmean(mergedact,3),[53 63 52 75]);
for idx = 1:numel(x)
    % correlating BEHAVIOUR dissimilarity with voxel
    [Rbeha(x(idx),y(idx),z(idx)), Pbeha(x(idx),y(idx),z(idx))]  = nancorrcoef([dpsyn;dpsyw;dpsynw],   dmedian(x(idx),y(idx),z(idx),[nid; wid;nwid]));
    [Rbehs(x(idx),y(idx),z(idx)), Pbehs(x(idx),y(idx),z(idx))]  = nancorrcoef(dsemfeat,dmedian(x(idx),y(idx),z(idx),wid));    
    disp(idx);
end

save normalised125_RSA_comb
return
%%
pvoxmod = NaN(53,63,52); Mpvoxmod = nanmean(Rvmod,4);
for x = 1:53
    for y = 1:63
        for z = 1:52
            if ~isnan(nanmean(squeeze(Rvmod(x,y,z,:))))                
                [pvoxmod(x,y,z)] = signrank(squeeze(Rvmod(x,y,z,:)));
            end
        end
    end
    disp(x)
end


for sub = 1:nsubs
    mask = reshape(L2_str.mergedevtbetash{sub},[53 63 52]);   mask(~isnan(mask)) = 1;
    MASK(:,:,:,sub) = mask;
end
mask = nanmean(MASK,4); mask(mask > 0) = 1;

Folder = 'RAW';
% Mpvoxmod(pvoxmod > .05 | mask == 0) = nan;

% saving different matrix as a mask
V = spm_vol('..\preprocessing\20190303_LDT_SUB01\glm\wdnlocmerged\spmT_0001.nii');
V = rmfield(V,'pinfo');


V.fname = ['searchlight\' Folder '\psy_W.nii'];         spm_write_vol(V,Rbeha);
V.fname = ['searchlight\' Folder '\psy_S.nii'];         spm_write_vol(V,Rbehs);
V.fname = ['searchlight\' Folder '\ractRT.nii'];        spm_write_vol(V,nanmedian(ractRT,4));
V.fname = ['searchlight\' Folder '\modeldiff.nii'];     spm_write_vol(V,nanmedian(Rvmod,4));

return
