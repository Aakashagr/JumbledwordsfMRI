if ~exist('L2_str'); load L2fmri_LDTw; end
allclearL2

%%
allclearL2
BETAs = L2_str.mergedevtbeta;
qw = 1:32; qn = 33:64; ql = 65:74;
RT = L2_str.RT;
nsubs = numel(L2_str.subjectid);
[idx, ROIname] = getvoxind(L2_str);

%% Figure S16A-F: Plotting mean activity values
meanRTn = nanmean(nanmean(RT(qn,:,1:8),3),2);  % Mean dissimilarities are better correlated with activity
meanRTw = nanmean(nanmean(RT(qw,:,1:8),3),2);

cnt = 0;
for roi = 1:5
    cnt = cnt + 1;
    for sub = 1:numel(L2_str.subjectid)
        betas = L2_str.mergedevtbeta{sub};
        if roi == 4; max_vox = 40; else, max_vox = Inf; end
        nvox = min(numel(idx{sub,roi}),max_vox);
        
        actw(roi,sub,:)  = nanmean(betas(idx{sub,roi}(1:nvox),qw));  % mean words responses
        actn(roi,sub,:)  = nanmean(betas(idx{sub,roi}(1:nvox),qn));  % mean nonword responses
        actl(roi,sub,:)  = nanmean(betas(idx{sub,roi}(1:nvox),ql));  % mean single letter responses
    end
    
    meanactw = nanmean(nanmean(actw(roi,:,:),3),2);
    meanactn = nanmean(nanmean(actn(roi,:,:),3),2);
    meanactl = nanmean(nanmean(actl(roi,:,:),3),2);
    
    semactw  = nansem(nanmean(actw(roi,:,:),3)');
    semactn  = nansem(nanmean(actn(roi,:,:),3)');
    semactl  = nansem(nanmean(actl(roi,:,:),3)');
    
    figure(1); subplot(1,5,cnt); barweb([meanactw meanactn meanactl], [semactw semactn semactl]);title(ROIname{roi}); legend('Words','Nonwords','Letters','Location','Best'); ylabel('Average Betas');
end

%% Regressing out RTs
cnt = 0;
meanRTn = nanmean(RT(qn,:,1:8),3);  % Mean dissimilarities are better correlated with activity
meanRTw = nanmean(RT(qw,:,1:8),3);


for roi = 1:5
    cnt = cnt + 1;
    for sub = 1:numel(L2_str.subjectid)
        betas = L2_str.mergedevtbeta{sub};
        if roi == 4; max_vox = 40; else, max_vox = Inf; end
        nvox = min(numel(idx{sub,roi}),max_vox);
        act = nanmean(betas(idx{sub,roi}(1:nvox),[qw; qn]));
        
        Xmat = [[meanRTw(:,sub);meanRTn(:,sub)], ones(64,1)];
        w = regress(act',Xmat);
        res(sub,:)  = act' - Xmat(:,1)*w(1);
    end
    
    P(roi) = signrank(nanmean(res(:,qn),2),nanmean(res(:,qw),2));
    
    meanactw(roi,1) = nanmean(nanmean(res(:,qw),2));
    meanactn(roi,1) = nanmean(nanmean(res(:,qn),2));
    
    semactw(roi,1) = nansem(nanmean(res(:,qw),2));
    semactn(roi,1) = nansem(nanmean(res(:,qn),2));
    
end

figure;  barweb([meanactw meanactn], [semactw semactn]);
xticklabels(ROIname); legend('Words','Nonwords','Location','Best'); ylabel('Average Betas');

%% Activity of transposed and substituted nonword
stimtype = vec(L2_str.stimtype');
qw1 = stimtype(1:16);  qn1 = stimtype(1:16)+32;  ql = 65:74;
qw2 = stimtype(17:32); qn2 = stimtype(17:32)+32; ql = 65:74;

for sub = 1:numel(L2_str.subjectid)
    betas = L2_str.mergedevtbeta{sub};
    max_vox = 40; 
    nvox = min(numel(idx{sub,roi}),max_vox);
    actn1(sub,1)  = nanmean(nanmean(betas(idx{sub,roi}(1:nvox),qn1360)),2);  % mean nonword responses
    actn2(sub,1)  = nanmean(nanmean(betas(idx{sub,roi}(1:nvox),qn2)),2);  % mean nonword responses
end

[mean(actn1) mean(actn2) nansem(actn1) nansem(actn2)]


%% Figure S16F: word/nonword classification

rng('default'); clear act
for roi = 1:5
    for sub = 1:numel(L2_str.subjectid)
        betas = L2_str.mergedevtbeta{sub};
        if roi == 4; max_vox = 40; else, max_vox = Inf; end
        nvox = min(numel(idx{sub,roi}),max_vox);
        
        act = betas(idx{sub,roi}(1:nvox),[qw1 qn1])'; act(:,isnan(mean(act))) = [];
        Mdl = fitcdiscr(act,[ones(16,1); ones(16,1)*2],'KFold',4); L1(roi,sub) = kfoldLoss(Mdl);
        
        act = betas(idx{sub,roi}(1:nvox),[qw2 qn2])'; act(:,isnan(mean(act))) = [];
        Mdl = fitcdiscr(act,[ones(16,1); ones(16,1)*2],'KFold',4); L2(roi,sub) = kfoldLoss(Mdl);
        
        act = betas(idx{sub,roi}(1:nvox),[qw1 qw2 qn1 qn2])'; act(:,isnan(mean(act))) = [];
        Mdl = fitcdiscr(act,[ones(32,1); ones(32,1)*2],'KFold',8); Lall(roi,sub) = kfoldLoss(Mdl);
        
    end
    disp(roi)
end
accuracy1 = 1-L1; accuracy2 = 1-L2; accuracyall = 1-Lall;
figure; barweb([mean(accuracy1,2) mean(accuracy2,2)],[nansem(accuracy1,2) nansem(accuracy2,2)]);
hold on; plot(1:5, ones(5,1)*.5,'k-'); ylabel('4 fold CV Accuracy'); title('Word vs nonword classification');
legend('Transposed pairs','Substituted pairs'); xticklabels({'V1-V3','V4','LOC','VWFA','TG'})
for i = 1:5
    [p(i,1)] = signrank(accuracy1(i,:)-.5);
    [p(i,2)] = signrank(accuracy2(i,:)-.5);
    [pall(i,1)] = signrank(accuracyall(i,:)-.5);
    [p(i,3)] = signrank(accuracy1(i,:),accuracy2(i,:));
end

%% Figure S16H - Linear model
stimid = L2_str.stimid;
sl = 'ESAROLITND'-64;
for i = 1:numel(sl)
    stimid(stimid == sl(i)) = 100+i;
end
stimid = stimid - 100;
froi = {'EVC','EVC','LOC','VWFA','WSCR'};


nrep = 10;
rw = NaN(5,numel(L2_str.subjectid),1000,nrep); p = rw; Bw = repmat(rw,[1 1 1 1 6]); rn= rw; Bn = Bw;rwe= rw;rne= rw;
csrw = rw; csrn = rw; psrw = rw; psrn = rw;

for roi = 1:5
    for sub = 1:numel(L2_str.subjectid)
        for rep = 1:nrep
            if roi == 4; max_vox = 40; else, max_vox = Inf; end
            nvox = min(numel(idx{sub,roi}),max_vox);
            
            runbetas = L2_str.fROI.(froi{roi}).runbetas{sub};
            id = L2_str.fROI.(froi{roi}).ids{sub};
            cidx = find(ismember(id,idx{sub,roi}(1:nvox))); runbetas = runbetas(cidx,:,:);
            
            qrm = isnan(nanmean(nanmean(runbetas,3),2));
            runbetas(qrm,:,:) = [];
            
            tr_split = randperm(8,4);
            acts_tr  = nanmean(runbetas(:,[qw qn],tr_split),3);
            acts_ts  = nanmean(runbetas(:,[qw qn],setdiff(1:8, tr_split)),3);
            actl  = nanmean(runbetas(:,ql,:),3);
            
            for vox = 1:size(acts_tr,1)
                yobs = (acts_tr(vox,:)'); yl = actl(vox,:)';
                Xmat = [(yl(stimid(:,:))) ones(numel(yobs),1)];  w = regress(yobs, Xmat);
                ypred = Xmat*w;
                
                [rw(roi,sub,vox,rep),p(roi,sub,vox,rep)] = nancorrcoef(ypred(qw), acts_ts(vox,qw)');
                [rn(roi,sub,vox,rep),p(roi,sub,vox,rep)] = nancorrcoef(ypred(qn), acts_ts(vox,qn)');
                [csrw(roi,sub,vox,rep),psrw(roi,sub,vox,rep)] = nancorrcoef(acts_ts(vox,qw),acts_tr(vox,qw));
                [csrn(roi,sub,vox,rep),psrn(roi,sub,vox,rep)] = nancorrcoef(acts_ts(vox,qn),acts_tr(vox,qn));
                
                Bw(roi,sub,vox,rep,1:numel(w)) = w;
            end
        end
    end
    disp(roi)
end

Crw = rw./csrw; Crn = rn./csrn;
figure; barweb([nanmedian(nanmedian(nanmedian(Crw,4),3),2) nanmedian(nanmedian(nanmedian(Crn,4),3),2)],[nansem(nanmedian(nanmedian(Crw,4),3),2) nansem(nanmedian(nanmedian(Crn,4),3),2)]); xticklabels(ROIname);
legend('Words','Nonwords'); ylabel('noise corrected median model fit')
