%% Analysing single letter dissimilarities
allclear
load L2_letter

L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.
RT = rmRToutlier1(L2_str.RT,1.7, 5);
uppercase = find(ismember(L2_str.img_pairs, nchoosek(1:26,2),'rows'));
lowercase = find(ismember(L2_str.img_pairs, nchoosek(27:52,2),'rows'));
digits = find(ismember(L2_str.img_pairs, nchoosek(53:62,2),'rows'));

images = L2_str.images;
for i = 1:numel(images)
    images{i}(:,[1:200 350:end]) = [];
end

% Calculating split-half between odd and even subjects
nancorrcoef(nanmean(nanmean(RT(:,1:2:end,:),3),2),nanmean(nanmean(RT(:,2:2:end,:),3),2));
acc = mean(L2_str.PC);
%% Fig 2B - MDS visualization for uppercase
dpsy = 1./nanmean(nanmean(RT,3),2);
D = squareform(dpsy(uppercase)); 
[Y,e] = mdscale(D,2);
simplot(dpsy(uppercase),images(1:26),.07);  % Upper case
hold on
plot(Y(:,1),ones(size(Y,1),1)*-.7,'b+')
plot(ones(size(Y,1),1)*-.7,Y(:,2),'r+')
figure; plot(Y(:,1)); xlim([1 26]); hold on
plot(Y(:,2)); xlim([1 26])
%% Fig 2D
RT_letter = rmRToutlier1(L2_str.RT,1.7, 5);
uppercase = find(ismember(L2_str.img_pairs, nchoosek(1:26,2),'rows'));
RT = RT_letter(uppercase,:,:); 
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


