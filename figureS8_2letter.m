allclear
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
lid = [1 4 8 9 13 14 20]; % Letter id
lpairs = find(ismember(nchoosek(1:26,2),lid(nchoosek(1:7,2)),'rows'));

load L2_bigram.mat
RT = rmRToutlier1(L2_str.RT,5, 2);
images = L2_str.images;
image_pairs = L2_str.img_pairs;

mRT = mean(RT,3);
RT_mean = nanmean(mRT,2);

% segregating high and low frequency bigrams & stimuli that form words
word_id   = [ 2, 5, 3, 6,7,18,15,28,27]; % actual word
hfword_id = [45,37,42,29,46, 6,7,18,15,28,27]; % high frequency bigram word
pword_id  = [46,43,24,25,23,26,8,11,36,39,1,4,29,32]; % pseudo word
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

%% Fig S8B
[y_pred,b,bint,X] = Batonmodel(1./RT_mean, 7, 2,1,1);
figure; corrplot(y_pred,1./RT_mean,'Upright bigram',1); hold all;
h(1) = plot(y_pred(word_pairs),1./RT_mean(word_pairs),'rd');
h(2) = plot(y_pred(hfword_pairs),1./RT_mean(hfword_pairs),'bo');
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity');
legend(h,{'word-word pairs','High-freq bigram'});
xlim([0.2 1.6]); ylim([0.2 1.6])

% Calculating error
avg_freq = mean(L2_str.bigram_freq(nchoosek(1:49,2)),2);
[~,fid] = sort(avg_freq);
error = abs(1./RT_mean - y_pred);
error_lf = error(fid(1:20)); 
error_hf = error(fid(end:-1:end-19));
Merror_lf = mean(error_lf); Serror_lf = std(error_lf);
Merror_hf = mean(error_hf); Serror_hf = std(error_hf);
ranksum(error_lf,error_hf)
%% Fig S8C; Correlating corresponding terms with across and within terms
figure;
corrplot(b(1:21), b(43:63)); xlim([0 .4]); ylim([-.2 .3]); hold all
corrplot(b(1:21), b(22:42),[],[],'+');xlim([0 .4]); ylim([-.2 .3]); unityslope
corrplot(b(1:21), b(64:84));xlim([0 .4]); ylim([-.2 .3]);

%% Fig S8D
nancorrcoef(b(1:21), dis_uletter(lpairs));
figure; simplot(b(1:21), L2_str_letter.images(lid))

%% Fig S8F; comparing both part sum and ISI models
for S = 1:30
    q1 = randperm(8,4); q2 = setdiff(1:8,q1);
    dis1 = 1./nanmean(mRT(:,q1),2); dis2 = 1./nanmean(mRT(:,q2),2);
    dpred1 = Batonmodel(dis1, 7, 2,1,1);
    dpred2 = Batonmodel(dis1, 7, 2,3,1,[],[],dis_uletter(lpairs));
    R(S,1) = nancorrcoef(dpred1,dis2);
    R(S,2) = nancorrcoef(dpred2,dis2);
    SH(S,1)= nancorrcoef(dis1,dis2);
end
  
figure; barweb(mean(R),std(R)); ylabel('Correlation Coefficient'); hold on
shadedErrorBar([.7; 1.3],[mean(SH); mean(SH)], std(SH))

