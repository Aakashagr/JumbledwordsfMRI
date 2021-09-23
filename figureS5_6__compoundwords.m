allclear
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
D = squareform(dis_uletter); ncells = 10;
opts.MaxIter = 400;
[rates6,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses

% 6 letter
load L2cb_6letter; L2cb_6 = L2_str; clear L2_str
acc = mean(L2cb_6.PC);
RT6 = rmRToutlier1(L2cb_6.RT,3, 3); dobs6 = 1./nanmean(nanmean(RT6,3),2); 
wwpairs = find(ismember(nchoosek(1:36,2), nchoosek([1 8 14 15 22 29 30],2),'rows'));
nnpairs = find(ismember(nchoosek(1:36,2), nchoosek(setdiff(1:36,[1 8 14 15 22 29 30]),2),'rows'));
csh = spearmanbrowncorrection(splithalfcorrd(nanmean(RT6,3)'),2);

% 3 letter
load L2cb_3letter; L2cb_3 = L2_str; clear L2_str
w = 1:66; nw = 67:132;
RT3 = rmRToutlier1(L2cb_3.RT,3, 3); dis3 = 1./nanmean(nanmean(RT3(w,:,:),3),2);
D = squareform(dis3); 
opts.MaxIter = 400;
[rates3,e] = mdscale(D,ncells,'Options',opts,'Start','random'); % rates = nstim x ncells matrix of single letter responses


% Data consistency
nancorrcoef(nanmean(nanmean(RT6(:,1:2:end,:),3),2),nanmean(nanmean(RT6(:,2:2:end,:),3),2));
nancorrcoef(nanmean(nanmean(RT3(w,1:2:end,:),3),2),nanmean(nanmean(RT3(w,2:2:end,:),3),2));
nancorrcoef(nanmean(nanmean(RT3(nw,1:2:end,:),3),2),nanmean(nanmean(RT3(nw,2:2:end,:),3),2));

%% Figure S5C
imgparts = L2cb_6.letterind; srchpairs = nchoosek(1:36,2);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
X6 = {allparts,rates6};  nparts = size(imgparts,2); 
opts = statset('nlinfit'); % opts.RobustWgtFun = 'bisquare'; opts.MaxIter = 1000;
% rng('default')
w0 = rand(ncells*nparts+1,1); 
west = nlinfit(X6,dobs6,@neuralmodel,w0,opts);
dpred = neuralmodel(west,X6); 

figure; corrplot(dpred,dobs6,'Upright bigram',1); hold all;
hold on;h = plot(dpred(wwpairs), dobs6(wwpairs),'ro'); 
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity');legend(h,'Word-Word pairs');
xlim([0.1 1.1]); ylim([0.1 1.1])

% OLD measure
for i = 1:size(srchpairs,1)
    dist(i)=editDistance(char(imgparts(srchpairs(i,1),:)+64),char(imgparts(srchpairs(i,2),:)+64),'SubstituteCost',2);
end

%% Figure S5DE
for rep = 1:1000
    q1 = randperm(size(L2cb_6.RT,2),size(L2cb_6.RT,2)/2); q2 = setdiff(1:size(L2cb_6.RT,2),q1);
    dis6_1 = 1./nanmean(nanmean(RT6(:,q1,:),3),2); dis6_2 = 1./nanmean(nanmean(RT6(:,q2,:),3),2);    
    
    imgparts = L2cb_6.letterind; srchpairs = nchoosek(1:36,2); 
    allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
    X6 = {allparts,rates6}; nparts = size(imgparts,2); w0 = rand(ncells*nparts+1,1);
    west6(rep,:) = nlinfit(X6,dis6_1,@neuralmodel,w0,opts); dpred = neuralmodel(west6(rep,:),X6); 
    [r(1,rep), p(1,rep)] = nancorrcoef(dpred, dis6_2);
    
    c(rep,1) = nancorrcoef(dpred(nnpairs), dis6_2(nnpairs));
    c(rep,2) = nancorrcoef(dpred(wwpairs), dis6_2(wwpairs));
    
    imgparts = L2cb_6.wordind; srchpairs = nchoosek(1:36,2);
    allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
    X3 = {allparts,rates3}; nparts = size(imgparts,2); w0 = rand(ncells*nparts+1,1);
    west = nlinfit(X3,dis6_1,@neuralmodel,w0,opts); dpred = neuralmodel(west,X3); 
    [r(2,rep), p(2,rep)] = nancorrcoef(dpred, dis6_2);
    [r(3,rep), p(3,rep)] = nancorrcoef(dist, dis6_2);
    [r(4,rep), p(4,rep)] = nancorrcoef(dis6_1, dis6_2); 
    disp(rep)
end

% Figure S5D
figure;
barweb(mean(r(1:3,:),2),std(r(1:3,:),[],2)); hold all; ylim([0 0.6]);
shadedErrorBar([.7; 1.3],[mean(r(4,:)); mean(r(4,:))], std(r(4,:)));  ylabel('Correlation Coefficient');

% Figure S5E
figure; barweb(mean(c),std(c)); xlabel('Model fits')


%% Figure S5F Comparing model performance on words and non-words
dis3_all = 1./nanmean(nanmean(RT3,3),2);

imgparts = L2cb_3.letterind(1:12,:); srchpairs = nchoosek(1:12,2);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
X3 = {allparts,rates6}; nparts = size(imgparts,2); w0 = ones(ncells*nparts+1,1);
west = nlinfit(X3,dis3_all(w),@neuralmodel,w0,opts); dpred_w = neuralmodel(west,X3);
    
imgparts = L2cb_3.letterind(13:24,:); srchpairs = nchoosek(1:12,2);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
X3 = {allparts,rates6}; nparts = size(imgparts,2); w0 = ones(ncells*nparts+1,1);
west = nlinfit(X3,dis3_all(nw),@neuralmodel,w0,opts); dpred_nw = neuralmodel(west,X3);
        

figure; corrplot(dpred_w,dis3_all(w),[],1); 
hold on; corrplot(dpred_nw,dis3_all(nw));    xlim([.4 1.1]); ylim([.4 1.1]);
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity')
%% Figure S6
figure
for neu = 1:ncells
    q = (neu-1)*6+1:neu*6;
    subplot(2,5,neu); shadedErrorBar(1:6,mean(west6(:,q)), std(west6(:,q)));xlim([1 6])
    title(['cell # ' num2str(neu)])
end

