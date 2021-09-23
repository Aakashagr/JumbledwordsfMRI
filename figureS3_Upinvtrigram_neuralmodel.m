allclear

load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
D = squareform(dis_uletter); ncells = 10;
opts.MaxIter = 400;
[rates,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses
%%
load L2_trigram.mat
RT = rmRToutlier1(L2_str.RT,4, 2);
images = L2_str.images;
image_pairs = L2_str.img_pairs(1:500,:);
nimg = 6; img_len = 3;
obj_parts = permn(1:nimg,img_len);

up = 1:500; inv = 501:1000;
mRT = nanmean(RT,3);
nancorrcoef(nanmean(mRT(up,1:2:end),2),nanmean(mRT(up,2:2:end),2));
nancorrcoef(nanmean(mRT(inv,1:2:end),2),nanmean(mRT(inv,2:2:end),2));

% Average across subjects
RT_up_mean = nanmean(mRT(up,:),2);
RT_inv_mean = nanmean(mRT(inv,:),2);

% Identifying word-word search pairs
words = [1 3 5; 1 3 6; 1 4 5; 2 1 2; 2 1 6; 3 1 2; 3 1 6;...
         4 1 2; 4 1 3; 4 1 5; 4 1 6; 5 1 2; 5 1 3; 5 1 4; 5 1 6];

% 3 letter transposition stimuli
trans1 = perms([4 3 2]);
trans2 = perms([5 1 4]);

% identifying pair locations in object parts list
words_id = find(ismember(obj_parts,words,'rows') == 1);
trans1_id = find(ismember(obj_parts,trans1,'rows') == 1);
trans2_id = find(ismember(obj_parts,trans2,'rows') == 1);

% word search pairs to be included in the stimulus set
words_pairs = words_id(nchoosek(1:numel(words_id),2));
trans1_pairs = trans1_id(nchoosek(1:numel(trans1_id),2));
trans2_pairs = trans2_id(nchoosek(1:numel(trans2_id),2));

% identifying whether seach pairs are already included in the list
words_pairs_id = find(ismember(image_pairs, words_pairs,'rows') == 1); 
trans1_pairs_id = find(ismember(image_pairs, trans1_pairs,'rows') == 1); 
trans2_pairs_id = find(ismember(image_pairs, trans2_pairs,'rows') == 1); 
trans = [trans1_pairs_id; trans2_pairs_id];

nw_pairs_id = setdiff(1:500, words_pairs_id);

for i = 1:500
    cimages{i} = [images{image_pairs(i,1)} images{image_pairs(i,2)}];
end

[~,ci] = splithalfcorrd(nanmean(RT(up,:,:),3)'); 
SH = spearmanbrowncorrection(ci,2);
%% Fig S3B
figure; corrplot(1./RT_up_mean,1./RT_inv_mean,[],1); hold on
h(1) = plot(1./RT_up_mean(words_pairs_id),1./RT_inv_mean(words_pairs_id),'or');
h(2) = plot(1./RT_up_mean(trans),1./RT_inv_mean(trans),'db');
ylim([0.2 1.6]); xlim([0.2 1.6])
%% Fig S3C; model fit
lid = [1 7 14 18 20 25]; 
allparts = [lid(obj_parts(image_pairs(:,1),:)) lid(obj_parts(image_pairs(:,2),:))]; 
X = {allparts,rates}; dobs = 1./nanmean(nanmean(RT(up,:,:),3),2); 
nparts = size(obj_parts,2); 

opts = statset('nlinfit'); % opts.RobustWgtFun = 'bisquare'; opts.MaxIter = 1000;
rng('default')
w0 = rand(ncells*nparts+1,1); 
[west,r,J] = nlinfit(X,dobs,@neuralmodel,w0,opts);
y_pred = neuralmodel(west,X); 

figure; corrplot(y_pred,dobs,'Upright stimuli',1); hold on
h(1) = plot(y_pred(words_pairs_id),dobs(words_pairs_id),'or');
h(2) = plot(y_pred(trans),dobs(trans),'db');
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity');
legend(h,{'Word-Word pairs','Transposition pairs'})
xlim([.2 1.5]);ylim([.2 1.5])

% Error
error = abs(y_pred- dobs);
word_error = mean(error(words_pairs_id));
trans_error = mean(error(trans));
nw_pairs_id = setdiff(1:500, [words_pairs_id; trans]);
others_error = mean(error(nw_pairs_id));


%% Fig S3D; Cross validated model fits
RTu = RT(up,:,:); RTi = RT(inv,:,:);
lid = [1 7 14 18 20 25]; 
allparts = [lid(obj_parts(image_pairs(:,1),:)) lid(obj_parts(image_pairs(:,2),:))]; 
nparts = size(obj_parts,2); 

% cross validated model fit
X = {allparts,rates};
opts = statset('nlinfit');  opts.UseParallel = 1; opts.MaxIter = 1000; 

for rep = 1:30
    % Upright trigrams
    q1 = randperm(9,4); q2 = setdiff(1:9,q1);
    d1 = 1./nanmean(nanmean(RTu(:,q1,:),3),2); d2 = 1./nanmean(nanmean(RTu(:,q2,:),3),2); 
    w0 = rand(ncells*nparts+1,1);
    westu(rep,:) = nlinfit(X,d1,@neuralmodel,w0,opts);
    dpred = neuralmodel(westu(rep,:),X);
    [ru(rep,1), pu(rep,1)] = nancorrcoef(dpred,d2);
    [ru(rep,2), pu(rep,2)] = nancorrcoef(d1,d2);
    
    % Inverted trigrams
    d1 = 1./nanmean(nanmean(RTi(:,q1,:),3),2); d2 = 1./nanmean(nanmean(RTi(:,q2,:),3),2);
    westi(rep,:) = nlinfit(X,d1,@neuralmodel,w0,opts);
    dpred = neuralmodel(westi(rep,:),X);
    [ri(rep,1), pi(rep,1)] = nancorrcoef(dpred,d2);
    [ri(rep,2), pi(rep,2)] = nancorrcoef(d1,d2);
    disp(rep)
end

figure;
barweb([mean(ru(:,1)) mean(ri(:,1))],[std(ru(:,1)) std(ri(:,1))]); hold all;
shadedErrorBar([1-.3 1] ,[mean(ru(:,2)) mean(ru(:,2))], std(ru(:,2)));
shadedErrorBar([1 1+.3] ,[mean(ri(:,2)) mean(ri(:,2))], std(ri(:,2)));
ylim([0 .9])

find(ri(:,1)./ri(:,2) > ru(:,1)./ru(:,2))

%% Fig S3EFG;  Comparing weights
Wu = reshape(westu(:,1:end-1),[30 3 ncells]); Wu(abs(Wu) < .001) = .001;
Wi = reshape(westi(:,1:end-1),[30 3 ncells]); Wi(abs(Wi) < .001) = .001;
MIu= []; MIi= []; cnt = 1;
for i = 1:30 
    MIu(cnt,:,1) = (Wu(i,1,:)-Wu(i,2,:))./ (Wu(i,1,:)+Wu(i,2,:));
    MIu(cnt,:,2) = (Wu(i,1,:)-Wu(i,3,:))./ (Wu(i,1,:)+Wu(i,3,:));
    MIu(cnt,:,3) = (Wu(i,3,:)-Wu(i,2,:))./ (Wu(i,3,:)+Wu(i,2,:));
    
    MIi(cnt,:,1) = (Wi(i,1,:)-Wi(i,2,:))./ (Wi(i,1,:)+Wi(i,2,:));
    MIi(cnt,:,2) = (Wi(i,1,:)-Wi(i,3,:))./ (Wi(i,1,:)+Wi(i,3,:));
    MIi(cnt,:,3) = (Wi(i,3,:)-Wi(i,2,:))./ (Wi(i,3,:)+Wi(i,2,:));
    cnt = cnt+1;
end
MIu(abs(MIu) > 1) = nan;
MIi(abs(MIi) > 1) = nan;
for i = 1:3
    figure; statcomparemean(nanmean(abs(MIu(:,:,i))),nanmean(abs(MIi(:,:,i))))
end
