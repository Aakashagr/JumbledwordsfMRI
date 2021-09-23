allclear

load L2_upinv_bigram.mat
RT = rmRToutlier(L2_str.RT,4, 2);
images = L2_str.images;
image_pairs = L2_str.img_pairs;
up = 1:630; inv = 631:1260;

mRT = reshape(RT,[size(RT,1) size(RT,2)*size(RT,3)]);

% Average across subjects
RT_up_mean = nanmean(mRT(up,:),2);
RT_inv_mean = nanmean(mRT(inv,:),2);
% figure; statcomparemean(RT_up_mean, RT_inv_mean); ylabel('Mean Reaction Time (s)');
% figure; corrplot(1./RT_up_mean, 1./RT_inv_mean,[],1); xlabel('Dissimilarity 1/RT, Upright stimuli'); ylabel('Dissimilarity 1/RT, Inverted stimuli')
[Upred,Ub,~] = Batonmodel(1./RT_up_mean, 6, 2,1,1);
[Ipred,Ib,~] = Batonmodel(1./RT_inv_mean, 6, 2,1,1);

% Panel C
figure; barweb([mean(reshape(Ub(1:end-1), [15 4]))', mean(reshape(Ib(1:end-1), [15 4]))'],...
               [nansem(reshape(Ub(1:end-1), [15 4]))',nansem(reshape(Ib(1:end-1), [15 4]))']);  ylim([-.3, .4])
xlabel('Parameters');  ylabel('Mean Estimated coefficient')
set(gca,'Xticklabel',{'C1','C2','Across','Within'},'XTickLabelRotation',45)


% panel B
[~,c] = splithalfcorrd(nanmean(RT(up,:,:),3)'); SH(:,1) = spearmanbrowncorrection(c,2);
[~,c] = splithalfcorrd(nanmean(RT(inv,:,:),3)'); SH(:,2) = spearmanbrowncorrection(c,2);
Rup = nancorrcoef(Upred, 1./RT_up_mean); Rin = nancorrcoef(Ipred, 1./RT_inv_mean);
figure; bar([Rup Rin]); hold on
for i = 1:2
    shadedErrorBar([i-.3 i+.3] ,[mean(SH(:,i)) mean(SH(:,i))], std(SH(:,i)));
end
ylabel('Correlation coefficient')

%% Text 
R = corr(reshape(Ub(1:end-1), [15 4]),reshape(Ib(1:end-1), [15 4])); % look at the diagonal values