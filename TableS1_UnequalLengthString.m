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

%%
load L2_string; 
% Unequal length data
% 6 vs 5
RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{5}; L5 = L2_str.letterind{6}; l5(:,[1 2 3 4 5]) = L5;   l5(:,6) = 63; 
imgparts = [L6;l5]; image_pairs = [1:50 ;51:100]'; RT = RT(501:550,:,:);
strlen = 6; [sh65(1,:), b65(1,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{5}; L5 = L2_str.letterind{6}; l5(:,[2 3 4 5 6]) = L5;   l5(:,1) = 63; 
imgparts = [L6;l5]; image_pairs = [1:50 ;51:100]'; RT = RT(501:550,:,:);
strlen = 6; [sh65(2,:), b65(2,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{5}; L5 = L2_str.letterind{6}; l5(:,[1 2 4 5 6]) = L5;   l5(:,3) = 63; 
imgparts = [L6;l5]; image_pairs = [1:50 ;51:100]'; RT = RT(501:550,:,:);
strlen = 6; [sh65(3,:), b65(3,:)] = partmodel_neu(imgparts, RT, image_pairs);
%%
% 6 vs 4
RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{7}; L4 = L2_str.letterind{8}; l4(:,[1 2 3 4]) = L4; l4(:,5:6) = 63;  
imgparts = [L6;l4]; image_pairs = [1:50 ;51:100]'; RT = RT(551:600,:,:);
strlen = 6; [sh64(1,:), b64(1,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{7}; L4 = L2_str.letterind{8}; l4(:,[3 4 5 6]) = L4; l4(:,1:2) = 63;  
imgparts = [L6;l4]; image_pairs = [1:50 ;51:100]'; RT = RT(551:600,:,:);
strlen = 6; [sh64(2,:), b64(2,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{7}; L4 = L2_str.letterind{8}; l4(:,[2 3 4 5]) = L4; l4(:,[1 6]) = 63;  
imgparts = [L6;l4]; image_pairs = [1:50 ;51:100]'; RT = RT(551:600,:,:);
strlen = 6; [sh64(3,:), b64(3,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L6 = L2_str.letterind{7}; L4 = L2_str.letterind{8}; l4(:,[1 2 5 6]) = L4; l4(:,3:4) = 63;  
imgparts = [L6;l4]; image_pairs = [1:50 ;51:100]'; RT = RT(551:600,:,:);
strlen = 6; [sh64(4,:), b64(4,:)] = partmodel_neu(imgparts, RT, image_pairs);

%% 5 vs 3 
[RT] = rmRToutlier1(L2_str.RT,3, 5); L5 = L2_str.letterind{9}; L3 = L2_str.letterind{10}; l3(:,[1 2 3]) = L3; l3(:,4:5) = 63; 
imgparts = [L5;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(601:650,:,:);
strlen = 5; [sh53(1,:),b53(1,:)] = partmodel_neu(imgparts, RT, image_pairs);

[RT] = rmRToutlier1(L2_str.RT,3, 5); L5 = L2_str.letterind{9}; L3 = L2_str.letterind{10}; l3(:,[3 4 5]) = L3; l3(:,1:2) = 63; 
imgparts = [L5;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(601:650,:,:);
strlen = 5; [sh53(2,:),b53(2,:)] = partmodel_neu(imgparts, RT, image_pairs);

[RT] = rmRToutlier1(L2_str.RT,3, 5); L5 = L2_str.letterind{9}; L3 = L2_str.letterind{10}; l3(:,[2 3 4]) = L3; l3(:,[1 5]) = 63; 
imgparts = [L5;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(601:650,:,:);
strlen = 5; [sh53(3,:),b53(3,:)] = partmodel_neu(imgparts, RT, image_pairs);

[RT] = rmRToutlier1(L2_str.RT,3, 5); L5 = L2_str.letterind{9}; L3 = L2_str.letterind{10}; l3(:,[1 2 5]) = L3; l3(:,3:4) = 63; 
imgparts = [L5;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(601:650,:,:);
strlen = 5; [sh53(4,:),b53(4,:)] = partmodel_neu(imgparts, RT, image_pairs);
%% 4 vs 3
RT = rmRToutlier1(L2_str.RT,3, 5); L4 = L2_str.letterind{11}; L3 = L2_str.letterind{12}; clear l3; l3(:,[1 2 3]) = L3; l3(:,4) = 63;
imgparts = [L4;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(651:700,:,:);
strlen = 4; [sh43(1,:), b43(1,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L4 = L2_str.letterind{11}; L3 = L2_str.letterind{12}; clear l3; l3(:,[ 2 3 4]) = L3; l3(:,1) = 63;
imgparts = [L4;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(651:700,:,:);
strlen = 4; [sh43(2,:), b43(2,:)] = partmodel_neu(imgparts, RT, image_pairs);

RT = rmRToutlier1(L2_str.RT,3, 5); L4 = L2_str.letterind{11}; L3 = L2_str.letterind{12}; clear l3; l3(:,[1 2 4]) = L3; l3(:,3) = 63;
imgparts = [L4;l3]; image_pairs = [1:50 ;51:100]'; RT = RT(651:700,:,:);
strlen = 4; [sh43(3,:), b43(3,:)] = partmodel_neu(imgparts, RT, image_pairs);
