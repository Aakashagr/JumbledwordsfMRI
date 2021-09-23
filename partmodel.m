function [sh, Rb, cp, C] = partmodel(imgparts, RT,image_pairs, Ldis)
global strlen
b0 = [.1 -.1 .2, -1 2]; opts = statset('nlinfit'); opts.RobustWgtFun = 'bisquare'; opts.MaxIter = 1000;
if ~exist('image_pairs'); image_pairs = []; end
load dis_letter; 
if ~exist('Ldis'); Ldis = dis_uletter; end
if numel(Ldis) > 340; nletters = 52;  else, nletters = 26; end
for S = 1:30
    q1 = randperm(size(RT,2), round(size(RT,2)/2)); q2 = setdiff(1:size(RT,2), q1);
    dis1 = 1./nanmean(nanmean(RT(:,q1,:),3),2); dis2 = 1./nanmean(nanmean(RT(:,q2,:),3),2);
    [pred_dis,b,~,X2] = Batonmodel(dis1, nletters, strlen,3,1,imgparts,image_pairs,Ldis);
    sh(1,S) = nancorrcoef(dis1,dis2);
    Rb(1,S) = nancorrcoef(dis2,pred_dis); % ISI model
    
    regmat = X2(:,varorder(strlen)); clear X
    f = fit(vec(1:6),vec([0.1255 0.0507 0.0397 0.0482 0.0408 0.0593]),'poly2');
    B = polyval([f.p1 f.p2 f.p3],linspace(1,6,strlen));
    for j = 1:strlen; X{j} = B(j)*exp(-2*[0 1:strlen-j]); X{j+strlen} = -X{j}(2:end); end
    W = cell2mat(X); pred_dis = regmat*vec(W);
    cp(1,S) = nancorrcoef(dis2,pred_dis); % SID on one experiment
    
    [beta,R,~,CovB,~,~] = nlinfit(regmat,dis1,@batonfun,b0,opts);
    [pred_dis,~] = nlpredci(@batonfun,regmat,beta,R,'Covar',CovB);
    C(1,S) = nancorrcoef(dis2,pred_dis); % SID model
end
