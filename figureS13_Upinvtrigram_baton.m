allclear
load L2_trigram.mat
RT = rmRToutlier1(L2_str.RT,4, 2);

images = L2_str.images;
image_pairs = L2_str.img_pairs(1:500,:);
nimg = 6; img_len = 3;
obj_parts = permn(1:nimg,img_len);

up = 1:500; inv = 501:1000;
mRT = mean(RT,3);

% Average across subjects
RT_up_mean = nanmean(mRT(up,:),2);
RT_inv_mean = nanmean(mRT(inv,:),2);
% figure; statcomparemean(RT_up_mean, RT_inv_mean);
% figure; corrplot(1./RT_up_mean, 1./RT_inv_mean,[],1); xlabel('Dissimilarity 1/RT, Upright stimuli'); ylabel('Dissimilarity 1/RT, Inverted stimuli')

[Upred,Ub,Ubint,X] = Batonmodel(1./RT_up_mean, nimg, 3,1,1,obj_parts,image_pairs);
[Ipred,Ib,Ibint] = Batonmodel(1./RT_inv_mean, nimg, 3,1,1,obj_parts,image_pairs);

% Figure S13B
% The error bars depends on whether subjects are sampled with replacement or not
[~,c] = splithalfcorrd(nanmean(RT(up,:,:),3)'); SH(:,1) = spearmanbrowncorrection(c,2); 
[~,c] = splithalfcorrd(nanmean(RT(inv,:,:),3)'); SH(:,2) = spearmanbrowncorrection(c,2);
Rup = nancorrcoef(Upred, 1./RT_up_mean); Rin = nancorrcoef(Ipred, 1./RT_inv_mean);
figure; bar([Rup Rin]); hold on
for i = 1:2
    shadedErrorBar([i-.3 i+.3] ,[mean(SH(:,i)) mean(SH(:,i))], std(SH(:,i)));
end
ylabel('Correlation coefficient')

% Figure S13C
figure; barweb([mean(reshape(Ub(1:end-1), [15 9]))', mean(reshape(Ib(1:end-1), [15 9]))'],...
               [nansem(reshape(Ub(1:end-1), [15 9]))',nansem(reshape(Ib(1:end-1), [15 9]))']);  ylim([-.2, .3])
xlabel('Parameters');  ylabel('Mean Estimated coefficient')
set(gca,'Xticklabel',{'C1','C2','C3','A-near1', 'A-near2', 'A-far','W-near1', 'W-near2', 'W-far'},'XTickLabelRotation',45)


%%
binv= reshape(Ib(1:end-1),[],9);
bup= reshape(Ub(1:end-1),[],9);
for i =1:9
    [r(i) p(i)] = nancorrcoef(binv(:,i),bup(:,i));
end


