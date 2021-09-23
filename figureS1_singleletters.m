%% Analysing single letter dissimilarities
allclear
load L2_letter

L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.
RT = rmRToutlier1(L2_str.RT,1.7, 5);
images = L2_str.images;
for i = 1:numel(images)
    images{i}(:,[1:200 350:end]) = [];
end

%% MDS visualization
dpsy = 1./nanmean(nanmean(RT,3),2);
D = squareform(dpsy); 
[Y,e] = mdscale(D,2);
dpred = pdist(Y); 
[cmatch,pmatch]=nancorrcoef(dpsy,dpred); 
% figure; simplot(dpsy,images,.07);  % Upper case

%%
dletters = 1./nanmean(nanmean(RT,3),2);
[~, shi] = splithalfcorrd(nanmean(RT,3)');
csh = spearmanbrowncorrection(shi,2);

ndimarray = [1:15]; 
D = squareform(dletters); 
opts.MaxIter = 400;

for ndim = ndimarray
    [Y,e(ndim)] = mdscale(D,ndim,'Options',opts); dpred = pdist(Y); 
    [cmatch(ndim,1),pmatch(ndim,1)]=nancorrcoef(dletters,dpred); 
end

figure; plot(ndimarray,cmatch); hold on; 
shadedErrorBar(ndimarray,mean(csh)*ones(size(ndimarray)),std(csh)*ones(size(ndimarray))); 
xlabel('Number of MDS dimensions'); ylabel('Correlation with Observed distances'); 
min(find(cmatch>=mean(csh)))
ylim([0.5 1]); xlim([1 15])

%% match with pixel distance
clear IMG
for i = 1:numel(images)
    IMG(i,:) = double(vec(images{i}));
    v1_img(i,:) = vec(v1model(images{i}));
end

figure; corrplot(pdist(IMG,'spearman'), dletters);
figure; corrplot(pdist(v1_img,'spearman'), dletters);


%%
allclear
load dis_letter
[num,txt,raw] = xlsread('uppercase_subjectiverating.xlsx');

Udis = num(:,5); % Storing single letter dissimilarities from the excel
% Getting indices of the stimuli pairs
stimpairs = txt(2:end,[3 4]); 
stimids = cellfun(@(x)  x-64, stimpairs);

% Selecting pairs used in visual search data
idx = nchoosek(1:26,2);
for i = 1:size(idx,1)
    q = find(ismember(stimids,[idx(i,:);fliplr(idx(i,:))],'rows'));
    SR(i,1) = mean(Udis(q));
end
corrplot(7-SR,dis_uletter,'Uppercase');
xlabel('Subjective ratings'); ylabel('Search dissimilarities')

%%
allclear
load dis_letter
num = xlsread('letter_confusability.xlsx');
pairs = nchoosek(1:26,2);

for i = 1:size(pairs,1)
    RT(i,:) = [num(pairs(i,1),pairs(i,2)), num(pairs(i,2),pairs(i,1))];
end
figure; corrplot(1./mean(RT,2),dis_uletter)

