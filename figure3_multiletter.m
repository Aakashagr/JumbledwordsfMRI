allclear
% Loading single letter data
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
D = squareform(dis_uletter); ncells = 10;
opts.MaxIter = 400;
[rates,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses

%% Loading bigram data
load L2_bigram.mat
images = L2_str.images;
image_pairs = L2_str.img_pairs;
acc = mean(L2_str.PC);
bigram_freq = L2_str.bigram_freq; % Bigram frequency counts

% segregating high and low frequency bigrams
% stimuli that forms words
word_id   = [ 2, 5, 3, 6,7,18,15,28,27]; % actual word
hfword_id = [45,37,42,29,46, 6,7,18,15,28,27]; % high frequency bigram word
pword_id  = [46,43,24,25,23,26,8,11,36,39,1,4,29,32]; % pseudo word
% nonword_id = setdiff(1:49,[word_id hfword_id pword_id]); % all other trigrams
nonword_id = setdiff(1:49,[word_id pword_id]); % all other trigrams

temp_word   = word_id(nchoosek(1:numel(word_id),2));
temp_hfword = hfword_id(nchoosek(1:numel(hfword_id),2));
temp_pword  = pword_id(nchoosek(1:numel(pword_id),2));
temp_others = nonword_id(nchoosek(1:numel(nonword_id),2));

word_pairs   = find(ismember(image_pairs,[temp_word; fliplr(temp_word)],'rows')==1);
hfword_pairs = find(ismember(image_pairs,[temp_hfword; fliplr(temp_hfword)],'rows')==1);
pword_pairs  = find(ismember(image_pairs,[temp_pword; fliplr(temp_pword)],'rows')==1);
nonword_pairs  = find(ismember(image_pairs,[temp_others; fliplr(temp_others)],'rows')==1);

word_others = find(ismember(image_pairs(:,1),word_id)|ismember(image_pairs(:,2),word_id));
nonword_others = find(ismember(image_pairs(:,1),nonword_id)|ismember(image_pairs(:,2),nonword_id));

% removing common pairs and word-word pairs & nonword-nonword pairs
temp1 = [intersect(nonword_others,word_others); intersect(word_pairs,word_others)];
temp2 = [intersect(nonword_others,word_others); intersect(nonword_pairs,nonword_others)];
word_others = setdiff(word_others,temp1);
nonword_others = setdiff(nonword_others,temp2);

%% Model fitting (Figure 3C)
tmp = [1 4 8 9 13 14 20];  imgparts = tmp(permn(1:7,2)); RT = rmRToutlier1(L2_str.RT,5,2); 

% Calculating split-half between odd and even subjects
nancorrcoef(nanmean(nanmean(RT(:,1:2:end,:),3),2),nanmean(nanmean(RT(:,2:2:end,:),3),2));
csh = spearmanbrowncorrection(splithalfcorrd(nanmean(RT,3)'),2);

srchpairs = nchoosek([1:49],2); 
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
X = {allparts,rates}; dobs = 1./nanmean(nanmean(RT,3),2); 
nparts = 2; 

opts = statset('nlinfit'); % opts.RobustWgtFun = 'bisquare'; opts.MaxIter = 1000;
% rng('default')
w0 = rand(ncells*nparts+1,1); 
[west,r,J] = nlinfit(X,dobs,@neuralmodel,w0,opts);
% ci = nlparci(west,r,'Jacobian',J);
dpred = neuralmodel(west,X); 

figure; corrplot(dpred,dobs,'Upright bigram',1); hold all;
h(1) = plot(dpred(word_pairs),dobs(word_pairs),'rd');
h(2) = plot(dpred(hfword_pairs),dobs(hfword_pairs),'bo');
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity');
legend(h,{'word-word pairs','High-freq bigram'});
xlim([0.2 1.8]); ylim([0.2 1.8])

%% Correlating error in model prediction with mean bigram frequency (Text)
avg_freq = mean(bigram_freq(nchoosek(1:49,2)),2);
[~,fid] = sort(avg_freq);
error = abs(dobs - dpred);
nancorrcoef(error,avg_freq);
error_lf = error(fid(1:20)); 
error_hf = error(fid(end:-1:end-19));
Merror_lf = mean(error_lf); Serror_lf = std(error_lf);
Merror_hf = mean(error_hf); Serror_hf = std(error_hf);
ranksum(error_lf,error_hf)
return

% Comparing difference in weights with random shuffle
w = reshape(west(1:end-1),[],10);
dobsR = dobs(randperm(length(dobs)));
[west_R] = nlinfit(X,dobsR,@neuralmodel,w0,opts);
wR = reshape(west_R(1:end-1),[],10);

w(abs(w) < .001) = .001; wR(abs(wR) < .01) = .001;
Wdiff = (w(1,:) - w(2,:));
Wdiff_R = (wR(1,:) - wR(2,:));
p = signrank(abs(Wdiff), abs(Wdiff_R))

%% %%%%%%%%       Upright vs Inverted bigrams (Experiment 3)      %%%%%%%%%%%%%%%
load L2_upinv_bigram.mat
RT = rmRToutlier1(L2_str.RT,4, 2);
% RT = L2_str.RT; % to generate full anova
img_pairs = L2_str.img_pairs;
up = 1:630; inv = 631:1260;
srch_pair = nchoosek(1:36,2);
tdparts = L2_str.partsTD;
qaabb = find(tdparts(:,1) == tdparts(:,2) & tdparts(:,3) == tdparts(:,4));
qabba = find(tdparts(:,1) == tdparts(:,4) & tdparts(:,3) == tdparts(:,2));
mRT = nanmean(RT,3);
%% Figure 3D
figure; barweb([mean(nanmean(mRT(up(qaabb),:))) mean(nanmean(mRT(inv(qaabb),:))); [mean(nanmean(mRT(up(qabba),:))) mean(nanmean(mRT(inv(qabba),:)))]]...
                    ,[nansem(nanmean(mRT(up(qaabb),:))') nansem(nanmean(mRT(inv(qaabb),:))');nansem(nanmean(mRT(up(qabba),:))') nansem(nanmean(mRT(inv(qabba),:))')])
ylabel('Observed Reaction Time, s');

% STATS
for sub = 1:8
    for qpairs = 1:numel(qaabb)
        for type = 1:2
            for rep = 1:2
                SUB(sub,qpairs,type,rep) = sub;
                QPAIRS(sub,qpairs,type,rep) = qpairs;
                TYPE(sub,qpairs,type,rep) = type;
                if type == 1
                    anova_RT(sub,qpairs,type,rep) = RT(up(qaabb(qpairs)),sub,rep);
%                     anova_RT(sub,qpairs,type,rep) = RT(up(qabba(qpairs)),sub,rep);
                else
                    anova_RT(sub,qpairs,type,rep) = RT(inv(qaabb(qpairs)),sub,rep);
%                     anova_RT(sub,qpairs,type,rep) = RT(inv(qabba(qpairs)),sub,rep);
                end
            end
        end
    end
end

p = anovan(vec(anova_RT),{vec(SUB), vec(QPAIRS), vec(TYPE)},'full');
%% Fig 3E
tmp = [1 12 14 18 19 20]; imgparts = tmp(permn(1:6,2)); nparts = size(imgparts,2); srchpairs = nchoosek([1:36],2);
RT = rmRToutlier1(L2_str.RT); RTu = RT(up,:,:); RTi = RT(inv,:,:);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 

% cross validated model fit
X = {allparts,rates};
opts = statset('nlinfit');  opts.UseParallel = 1; 

for rep = 1:30
    q1 = randperm(8,4); q2 = setdiff(1:8,q1); w0 = rand(ncells*nparts+1,1);
    d1 = 1./nanmean(nanmean(RTu(:,q1,:),3),2); 
    westu(rep,:) = nlinfit(X,d1,@neuralmodel,w0,opts);
    
    d1 = 1./nanmean(nanmean(RTi(:,q1,:),3),2); 
    westi(rep,:) = nlinfit(X,d1,@neuralmodel,w0,opts);
    disp(rep)
end

Wu = reshape(westu(:,1:end-1),[30 2 10]); Wu(abs(Wu) < .01) = .01;
Wi = reshape(westi(:,1:end-1),[30 2 10]); Wi(abs(Wi) < .01) = .01;
for i = 1:30
    MIu(i,:) = (Wu(i,1,:)-Wu(i,2,:))./ (Wu(i,1,:)+Wu(i,2,:));
    MIi(i,:) = (Wi(i,1,:)-Wi(i,2,:))./ (Wi(i,1,:)+Wi(i,2,:));
end

figure; statcomparemean(mean(abs(MIu)),mean(abs(MIi)))

%% %%%%%%%%%%    6 -letter strings (Experiment 4)      %%%%%%%%%%
load L2cb_6letter; L2cb_6 = L2_str;
acc = mean(L2cb_6.PC);
RT6 = rmRToutlier1(L2cb_6.RT,3, 3); dobs6 = 1./nanmean(nanmean(RT6,3),2); 
nancorrcoef(nanmean(nanmean(RT6(:,1:2:end,:),3),2),nanmean(nanmean(RT6(:,2:2:end,:),3),2));
wwpairs = find(ismember(nchoosek(1:36,2), nchoosek([1 8 14 15 22 29 30],2),'rows'));
nnpairs = find(ismember(nchoosek(1:36,2), nchoosek(setdiff(1:36,[1 8 14 15 22 29 30]),2),'rows'));
csh = spearmanbrowncorrection(splithalfcorrd(nanmean(RT6,3)'),2);

%% Figure 3F
imgparts = L2cb_6.letterind; srchpairs = nchoosek(1:36,2);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
X6 = {allparts,rates};  nparts = size(imgparts,2); 
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

%% Figure 3G Cross validated model fits
for rep = 1:1000
    q1 = randperm(size(L2cb_6.RT,2),size(L2cb_6.RT,2)/2); q2 = setdiff(1:size(L2cb_6.RT,2),q1);
    dis6_1 = 1./nanmean(nanmean(RT6(:,q1,:),3),2); dis6_2 = 1./nanmean(nanmean(RT6(:,q2,:),3),2);    
    
    imgparts = L2cb_6.letterind; srchpairs = nchoosek(1:36,2); 
    allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
    X6 = {allparts,rates}; nparts = size(imgparts,2); w0 = rand(ncells*nparts+1,1);
    west6(rep,:) = nlinfit(X6,dis6_1,@neuralmodel,w0,opts); dpred = neuralmodel(west6(rep,:),X6); 
    [r(1,rep), p(1,rep)] = nancorrcoef(dpred, dis6_2);
    
    c(rep,1) = nancorrcoef(dpred(nnpairs), dis6_2(nnpairs));
    c(rep,2) = nancorrcoef(dpred(wwpairs), dis6_2(wwpairs));
    
    [r(2,rep), p(2,rep)] = nancorrcoef(dist, dis6_2);
    [r(3,rep), p(3,rep)] = nancorrcoef(dis6_1, dis6_2); 
    disp(rep)
end

figure;
barweb(mean(r(1:2,:),2),std(r(1:2,:),[],2)); hold all; ylim([0 0.6]);
shadedErrorBar([.7; 1.3],[mean(r(3,:)); mean(r(3,:))], std(r(3,:)));  ylabel('Correlation Coefficient');

% Figure 3H
figure; barweb(mean(c),std(c)); xlabel('Model fits')
