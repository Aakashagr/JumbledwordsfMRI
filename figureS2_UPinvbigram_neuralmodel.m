allclear
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
D = squareform(dis_uletter); ncells = 10;
opts.MaxIter = 400;
[rates,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses
%% Fig S2C
load L2_upinv_bigram.mat
RT = rmRToutlier1(L2_str.RT,4, 2);
% RT = L2_str.RT; % to generate full anova
images = L2_str.images;
img_pairs = L2_str.img_pairs;
up = 1:630; inv = 631:1260;
srch_pair = nchoosek(1:36,2);
tdparts = L2_str.partsTD;
qaabb = find(tdparts(:,1) == tdparts(:,2) & tdparts(:,3) == tdparts(:,4));
qabba = find(tdparts(:,1) == tdparts(:,4) & tdparts(:,3) == tdparts(:,2));
mRT = nanmean(RT,3);

% Average across subjects
RT_up_mean = nanmean(mRT(up,:),2);
RT_inv_mean = nanmean(mRT(inv,:),2);
figure; barweb([mean(nanmean(mRT(up(qaabb),:))) mean(nanmean(mRT(inv(qaabb),:))); [mean(nanmean(mRT(up(qabba),:))) mean(nanmean(mRT(inv(qabba),:)))]]...
                    ,[nansem(nanmean(mRT(up(qaabb),:))') nansem(nanmean(mRT(inv(qaabb),:))');nansem(nanmean(mRT(up(qabba),:))') nansem(nanmean(mRT(inv(qabba),:))')])
ylabel('Mean Reaction Time (s)');

% STATS
for sub = 1:8
    for qpairs = 1:numel(qaabb)
        for type = 1:2
            for rep = 1:2
                SUB(sub,qpairs,type,rep) = sub;
                QPAIRS(sub,qpairs,type,rep) = qpairs;
                TYPE(sub,qpairs,type,rep) = type;
                if type == 1
                    anova_RT(sub,qpairs,type,rep) = RT(up(qaabb(qpairs)),sub,rep);
%                     anova_RT(sub,qpairs,type,rep) = RT(up(qabba(qpairs)),sub,rep);
                else
                    anova_RT(sub,qpairs,type,rep) = RT(inv(qaabb(qpairs)),sub,rep);
%                     anova_RT(sub,qpairs,type,rep) = RT(inv(qabba(qpairs)),sub,rep);
                end
            end
        end
    end
end

p = anovan(vec(anova_RT),{vec(SUB), vec(QPAIRS), vec(TYPE)},'full');
%% Fig S2D
figure; corrplot(1./RT_up_mean, 1./RT_inv_mean,[],1); hold on
plot(1./RT_up_mean(qaabb), 1./RT_inv_mean(qaabb),'or'); plot(1./RT_up_mean(qabba), 1./RT_inv_mean(qabba),'ob')
xlabel('Dissimilarity 1/RT, Upright stimuli'); ylabel('Dissimilarity 1/RT, Inverted stimuli'); ylim([.2 2]);xlim([.2 2])

%% Fig S2E
tmp = [1 12 14 18 19 20]; imgparts = tmp(permn(1:6,2)); nparts = size(imgparts,2); srchpairs = nchoosek([1:36],2);
RT = rmRToutlier1(L2_str.RT); RTu = RT(1:630,:,:); RTi = RT(631:end,:,:);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 

% cross validated model fit
X = {allparts,rates};
opts = statset('nlinfit');  opts.UseParallel = 1; %opts.MaxIter = 1000;

for rep = 1:30
    q1 = randperm(8,4); q2 = setdiff(1:8,q1);
    d1 = 1./nanmean(nanmean(RTu(:,q1,:),3),2); d2 = 1./nanmean(nanmean(RTu(:,q2,:),3),2); 
    w0 = rand(ncells*nparts+1,1);
    westu(rep,:) = nlinfit(X,d1,@neuralmodel,w0,opts);
    dpred = neuralmodel(westu(rep,:),X);
    [ru(rep,1), pu(rep,1)] = nancorrcoef(dpred,d2);
    [ru(rep,2), pu(rep,2)] = nancorrcoef(d1,d2);
    
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

invratio = ri(:,1)./ri(:,2);
upratio = ru(:,1)./ru(:,2);
%% Fig S2GH; Comparing weights
Wu = reshape(westu(:,1:end-1),[30 2 10]); Wu(abs(Wu) < .01) = .01;
Wi = reshape(westi(:,1:end-1),[30 2 10]); Wi(abs(Wi) < .01) = .01;
for i = 1:30
    MIu(i,:) = (Wu(i,1,:)-Wu(i,2,:))./ (Wu(i,1,:)+Wu(i,2,:));
    MIi(i,:) = (Wi(i,1,:)-Wi(i,2,:))./ (Wi(i,1,:)+Wi(i,2,:));
end

barweb([mean(abs(MIu)); mean(abs(MIi))]',[std(abs(MIu)); std(abs(MIi))]');
xlabel('Neurons')
figure; statcomparemean(mean(abs(MIu)),mean(abs(MIi)))

%% Fig S2F 

tmp = [1 12 14 18 19 20]; imgparts = tmp(permn(1:6,2)); nparts = size(imgparts,2); srchpairs = nchoosek([1:36],2);
RT = rmRToutlier1(L2_str.RT); RTu = RT(1:630,:,:); RTi = RT(631:end,:,:);
allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 

% cross validated model fit
X = {allparts,rates};
opts = statset('nlinfit');  opts.UseParallel = 1; %opts.MaxIter = 1000;

q1 = randperm(8,4); q2 = setdiff(1:8,q1);
d1 = 1./nanmean(nanmean(RTu(:,q1,:),3),2); d2 = 1./nanmean(nanmean(RTu(:,q2,:),3),2); 
w0 = rand(ncells*nparts+1,1);
westu = nlinfit(X,d1,@neuralmodel,w0,opts);
RTpredu = 1./neuralmodel(westu,X);

d1 = 1./nanmean(nanmean(RTi(:,q1,:),3),2); d2 = 1./nanmean(nanmean(RTi(:,q2,:),3),2);
westi = nlinfit(X,d1,@neuralmodel,w0,opts);
RTpredi= 1./neuralmodel(westi,X);

figure; barweb([mean(RTpredu(qaabb)) mean(RTpredi(qaabb)); mean(RTpredu(qabba)) mean(RTpredi(qabba))]...
                    ,[nansem(RTpredu(qaabb)) nansem(RTpredi(qaabb)); nansem(RTpredu(qabba)) nansem(RTpredi(qabba))])
ylabel('Mean Reaction Time (s)');

