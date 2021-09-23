%% fitting a non-linear model and comparing coefficients with baton model.
allclear
global strlen ncells rates
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);

lowercase = find(ismember(L2_str_letter.img_pairs, nchoosek(27:52,2),'rows'));
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
mixedcase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:52,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
dis_lletter = 1./nanmean(nanmean(RT_letter(lowercase,:,:),3),2);
dis_mletter = 1./nanmean(nanmean(RT_letter(mixedcase,:,:),3),2);

% get the single cell shape tuning required
D = squareform(dis_mletter); ncells = 6;
[rates,e] = mdscale(D,ncells); % rates = nstim x ncells matrix of single letter responses

% 2 letter searches
load L2_bigram
tmp = [1 4 8 9 13 14 20];  imgparts = tmp(permn(1:7,2)); RT = rmRToutlier1(L2_str.RT,5,2); 
strlen = 2; [sh2(1,:), rm2(1,:), b2(1,:)] = partmodel_isi(imgparts, RT);

load L2_Rbigram % rank deficient
tmp = [1 5 9 16 19 20 25]; imgparts = tmp(permn(1:7,2)); imgparts = [imgparts(:,1) imgparts(:,2)]; RT = rmRToutlier1(L2_str.RT,5,2);
strlen = 2; [sh2(2,:), rm2(2,:), b2(2,:)] = partmodel_isi(imgparts, RT);

load L2_upinv_bigram
tmp = [1 12 14 18 19 20]; imgparts = tmp(permn(1:6,2)); RT = rmRToutlier1(L2_str.RT,4,2);  RT = RT(1:630,:,:);
strlen = 2; [sh2(3,:), rm2(3,:), b2(3,:)] = partmodel_isi(imgparts, RT);

% 3 letter searches
load L2_trigram.mat
image_pairs = L2_str.img_pairs(1:500,:); tmp = [1 7 14 18 20 25]; imgparts = tmp(permn(1:6,3)); RT = rmRToutlier1(L2_str.RT,4, 2);  RT = RT(1:500,:,:);
strlen = 3; [sh3(1,:), rm3(1,:), b3(1,:)] = partmodel_isi(imgparts, RT, image_pairs);

load L2cb_3letter.mat
RT = rmRToutlier1(L2_str.RT,2,3); imgparts = L2_str.letterind(1:12,:); RT = RT(1:66,:,:);
strlen = 3; [sh3(2,:), rm3(2,:), b3(2,:)] = partmodel_isi(imgparts, RT);
% non-words
RT = rmRToutlier1(L2_str.RT,2,3); imgparts = L2_str.letterind(13:24,:); RT = RT(67:132,:,:);
strlen = 3; [sh3(3,:), rm3(3,:), b3(3,:)] = partmodel_isi(imgparts, RT);

% 4 letters
load L2_4letter
RT = rmRToutlier1(L2_str.RT,4,4); imgparts = L2_str.letterind;  image_pairs = L2_str.img_pairs(1:300,:);
strlen = 4; [sh4(1,:), rm4(1,:), b4(1,:)] = partmodel_isi(imgparts, RT, image_pairs);

% string
load L2_6letter
RT = rmRToutlier1(L2_str.RT,5,4); imgparts = L2_str.letterind;  image_pairs = L2_str.img_pairs(1:300,:);
strlen = 6; [sh6(1,:), rm6(1,:), b6(1,:)] = partmodel_isi(imgparts, RT, image_pairs);

load L2cb_6letter
RT = rmRToutlier1(L2_str.RT,3,3); imgparts = L2_str.letterind;  
strlen = 6; [sh6(2,:), rm6(2,:), b6(2,:)] = partmodel_isi(imgparts, RT);

load L2_string; 
RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{1}; image_pairs = L2_str.img_pairs(1:300,:); RT = RT(1:300,:,:);
strlen = 6; [sh6(3,:), rm6(3,:), b6(3,:)] = partmodel_isi(imgparts, RT, image_pairs);
RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{13}+26; image_pairs = L2_str.img_pairs(701:800,:)-1400; RT = RT(701:800,:,:);
strlen = 6; [sh6(4,:), rm6(4,:), b6(4,:)] = partmodel_isi(imgparts, RT, image_pairs,dis_mletter);
RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{14}; image_pairs = L2_str.img_pairs(801:900,:)-1600; RT = RT(801:900,:,:);
strlen = 6; [sh6(5,:), rm6(5,:), b6(5,:)] = partmodel_isi(imgparts, RT, image_pairs,dis_mletter);

RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{2}; image_pairs = L2_str.img_pairs(301:400,:)-600; RT = RT(301:400,:,:);
strlen = 5; [sh5(1,:), rm5(1,:), b5(1,:)] = partmodel_isi(imgparts, RT, image_pairs);
RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{3}; image_pairs = L2_str.img_pairs(401:450,:)-800; RT = RT(401:450,:,:);
strlen = 4; [sh4(2,:), rm4(2,:), b4(2,:)] = partmodel_isi(imgparts, RT, image_pairs);
RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{4}; image_pairs = L2_str.img_pairs(451:500,:)-900; RT = RT(451:500,:,:);
strlen = 3; [sh3(4,:), rm3(4,:), b3(4,:)] = partmodel_isi(imgparts, RT, image_pairs);


% Unequal length data
RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{5}; L5 = L2_str.letterind{6}; l5(:,[1 2 4 5 6]) = L5;   l5(:,3) = 63; 
imgparts = [L6;l5]; image_pairs = [1:50 ;51:100]'; RT = RT(501:550,:,:);
strlen = 6; [shue(1,:), rmue(1,:), bue(1,:)] = partmodel_isi(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{7}; L4 = L2_str.letterind{8}; l4(:,[1 2 5 6]) = L4; l4(:,3:4) = 63;  
imgparts = [L6;l4]; image_pairs = [1:50 ;51:100]'; RT = RT(551:600,:,:);
strlen = 6; [shue(2,:), rmue(2,:), bue(2,:)] = partmodel_isi(imgparts, RT, image_pairs);

[RT,frac] = rmRToutlier1(L2_str.RT,3, 5); L5 = L2_str.letterind{9}; L3 = L2_str.letterind{10}; l3(:,[1 2 5]) = L3; l3(:,3:4) = 63; 
imgparts = [L5;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(601:650,:,:);
strlen = 5; [shue(3,:), rmue(3,:), bue(3,:)] = partmodel_isi(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L4 = L2_str.letterind{11}; L3 = L2_str.letterind{12}; clear l3; l3(:,[1 2 4]) = L3; l3(:,3) = 63;
imgparts = [L4;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(651:700,:,:);
strlen = 4; [shue(4,:), rmue(4,:), bue(4,:)] = partmodel_isi(imgparts, RT, image_pairs);


%%
order = [1 5 6 12 2 7 9 10 13 14 15 16 17 18 19 8 11 3 4];
% order = [1:7 9 8 10 13:15 12 11 16:19];
names = {'2L','2L(R)','2Lu','3Lu','3Lw','3Lnw','3Lr','4L','4Lr','5Lr','6L','6Lc','6Lr','6Ll','6Lm','UE65','UE64','UE53','UE43'};

% mean_sh = mean([sh2; sh3; sh4; sh5; sh6; shue],2);
% std_sh = std([sh2; sh3; sh4; sh5; sh6; shue],[],2);
% Mrm    = mean([rm2; rm3; rm4; rm5; rm6; rmue],2);  Srm = std([rm2; rm3; rm4; rm5; rm6; rmue],[],2);
mean_sh = mean([sh2; sh3; sh4; sh5; sh6; shue],2);
std_sh = std([sh2; sh3; sh4; sh5; sh6; shue],[],2);
brm    = mean([b2; b3; b4; b5; b6; bue],2);
bSrm = std([b2; b3; b4; b5; b6; bue],[],2);
Mrm    = mean([rm2; rm3; rm4; rm5; rm6; rmue],2);
Srm = std([rm2; rm3; rm4; rm5; rm6; rmue],[],2);


%%% Figure S7
% figure; bar(brm(order)); hold on
% errorbar(brm(order), bSrm(order));

%%% Figure S9
figure; barweb([brm(order) Mrm(order)],[bSrm(order) Srm(order)]); hold on

for i = 1:length(order)
    shadedErrorBar([i-.3 i+.3] ,[mean_sh(order(i)) mean_sh(order(i))], std_sh(order(i)),[]);
end

set(gca,'Xtick',1:length(order),'XtickLabel',names(order),'XTickLabelRotation', 45)
ylabel('Correlation coefficient'); legend('neural model','ISI model')
title('Model fit')

% NOTE: the barplots are rearranged in the paper 