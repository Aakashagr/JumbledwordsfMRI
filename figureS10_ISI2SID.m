%% Figure S10AB
allclear; global strlen
load L2_string; load dis_letter
RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{1}; image_pairs = L2_str.img_pairs(1:300,:); RT = RT(1:300,:,:);
strlen = 6; dis = 1./nanmean(nanmean(RT,3),2);
[~,b,~,~] = Batonmodel(dis, 26, strlen,3,1,imgparts,image_pairs,dis_uletter);
figure(1); plot(b(1:strlen)); hold on; plot(1:strlen, b(1:strlen),'*')
figure(2); q = b(varorder(strlen)); plot(q(1:strlen)); hold on; plot(1:strlen, q(1:strlen),'*')

RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{2}; image_pairs = L2_str.img_pairs(301:400,:)-600; RT = RT(301:400,:,:);
strlen = 5; dis = 1./nanmean(nanmean(RT,3),2);
[~,b,~,~] = Batonmodel(dis, 26, strlen,3,1,imgparts,image_pairs,dis_uletter);
figure(1);plot(b(1:strlen)); hold on; plot(1:strlen, b(1:strlen),'d')
figure(2); q = b(varorder(strlen)); plot(q(1:strlen)); hold on; plot(1:strlen, q(1:strlen),'d')

RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{3}; image_pairs = L2_str.img_pairs(401:450,:)-800; RT = RT(401:450,:,:);
strlen = 4; dis = 1./nanmean(nanmean(RT,3),2);
[~,b,~,~] = Batonmodel(dis, 26, strlen,3,1,imgparts,image_pairs,dis_uletter);
figure(1);plot(b(1:strlen)); hold on; plot(1:strlen, b(1:strlen),'s')
figure(2); q = b(varorder(strlen)); plot(q(1:strlen)); hold on; plot(1:strlen, q(1:strlen),'s')

RT = rmRToutlier1(L2_str.RT,3, 5); imgparts = L2_str.letterind{4}; image_pairs = L2_str.img_pairs(451:500,:)-900; RT = RT(451:500,:,:);
strlen = 3; dis = 1./nanmean(nanmean(RT,3),2);
[~,b,~,~] = Batonmodel(dis, 26, strlen,3,1,imgparts,image_pairs,dis_uletter);
figure(1); plot(b(1:strlen)); hold on; plot(1:strlen, b(1:strlen),'o')
figure(2); q = b(varorder(strlen)); plot(q(1:strlen)); hold on; plot(1:strlen, q(1:strlen),'o')

xlabel('String Length'); ylabel('Model coefficient')
%%
clear; global strlen
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);

uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);

load L2cb_6letter
[RT] = rmRToutlier1(L2_str.RT,3,3); imgparts = L2_str.letterind;
wwpairs = find(ismember(nchoosek(1:36,2), nchoosek([1 8 14 15 22 29 30],2),'rows'));

mRT = nanmean(RT,3); opts = statset('nlinfit'); opts.RobustWgtFun = 'bisquare'; opts.MaxIter = 1000;
strlen = 6; dis6 = 1./nanmean(mRT,2);
[dpred6,b6,b6int,X6] = Batonmodel(dis6, 26, 6,3,1,imgparts,[],dis_uletter);
Wisi = b6(varorder(6));

%% Fig S10CD Non-linear fit
opts = statset('nlinfit');  opts.RobustWgtFun = 'bisquare';  opts.MaxIter = 1000;
b0 = [.8 -.1 -.15 -2 2];
% b0 = rand(1,5);
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(X6(:,varorder(6)),dis6,@batonfun,b0,opts);
[Ypred,delta] = nlpredci(@batonfun,X6(:,varorder(6)),beta,R,'Covar',CovB);

figure; corrplot(Ypred, dis6,[],1);
hold on;h = plot(Ypred(wwpairs), dis6(wwpairs),'ro'); xlim([.1 1.1]);ylim([.1 1.1]);
xlabel('Predicted dissimilarity'); ylabel('Observed dissimilarity');legend(h,'Word-Word pairs');

clear X
for i = 1:strlen
    X{i} = [beta(1).*i^2 + beta(2).*i + beta(3)]*exp(beta(4)*[0 1:strlen-i]);
    X{i+strlen} = -X{i}(2:end);
end
Wsid = cell2mat(X);
figure; corrplot(Wsid,Wisi); xlabel('SID model coeff'); ylabel('ISI model coeff');

%% Fig S10E Comparing ISI and SID models
opts = statset('nlinfit');  opts.RobustWgtFun = 'bisquare';  opts.MaxIter = 1000; b0 = [.8 -.1 -.15 -2 2];
for S = 1:30
    q1 = randperm(8,4); q2 = setdiff(1:8,q1);
    dis1 = 1./nanmean(mRT(:,q1),2); dis2 = 1./nanmean(mRT(:,q2),2);
    [dpred1,~,~,X] = Batonmodel(dis1, 26, 6,3,1,L2_str.letterind,[],dis_uletter);    
    [beta,res,J,CovB,MSE,ErrorModelInfo] = nlinfit(X(:,varorder(6)),dis1,@batonfun,b0,opts);
    dpred2 = nlpredci(@batonfun,X(:,varorder(6)),beta,res,'Covar',CovB); 
    
    r(S,1) = nancorrcoef(dpred1,dis2);
    r(S,2) = nancorrcoef(dpred2,dis2);
    SH(S,1)= nancorrcoef(dis1,dis2);
    
end
figure; barweb(mean(r),std(r)); ylabel('Correlation Coefficient'); hold on
shadedErrorBar([.7; 1.3],[mean(SH); mean(SH)], std(SH))

