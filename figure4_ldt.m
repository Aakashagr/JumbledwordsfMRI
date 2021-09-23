allclear
% fitting a non-linear model and comparing coefficients with baton model.
global strlen ncells rates
load L2_letter
L2_str_letter = L2_str; clear L2_str
[RT_letter,frac] = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);

% get the single cell shape tuning required
D = squareform(dis_uletter); ncells = 10;
opts.MaxIter = 400;
[rates,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses

opts = statset('nlinfit'); opts.UseParallel = true; opts.MaxIter = 1000;
rng('default')
%%
load L2_ldt;
RT = L2_str.RT; PC = L2_str.PC; RT(PC < .3) = NaN; RT(RT < .3) = NaN;
RT(isoutlier(RT')') = nan;
% % Setting inaccurate RTs as nans
for i = 1:900; xx(i) = numel(find(isnan(RT(i,:)))); end
RT(xx > 5,:)  = NaN;

widx = L2_str.word_id';
wRT{1,1} = RT(1:100,:);    nwRT{1,1} = RT(101:200,:);
wRT{2,1} = RT(201:350,:);  nwRT{2,1} = RT(351:500,:);
wRT{3,1} = RT(501:700,:);  nwRT{3,1} = RT(701:900,:);

qwords = [1:100 201:350 501:700];
qnwords = [101:200 351:500 701:900];

%% extracting lexical features
for type = 1:3
    % predicting non-word RTs using lexical features
    [num, str] = xlsread('feat_nw_ldt.xlsx'); strlen = type+3; clear varmat
    nw = L2_str.word_id{strlen - 3}(50*(strlen -2)+1:end,:);
    for i = 1:length(nw); varmat(i,:) = num(ismember(str(2:end,2), char(nw(i,:)+64)),4:7); end
    varmat(isnan(varmat)) = 0; nwordfeat{type,1} = varmat;
    
    % Word features
    [num, str] = xlsread('feat_w_ldt.xlsx'); clear varmat
    w = L2_str.word_id{strlen - 3}(1:50*(strlen -2),:);
    for i = 1:length(w); varmat(i,:) = num(ismember(str(2:end,1), char(w(i,:)+64)),2:12); end
    varmat(isnan(varmat)) = 0; wordfeat{type,1} = varmat;
    
    % Letter frequency
    let_freq = [8.04 1.48 3.34 3.82 12.49 2.4 1.87 5.05 7.57 0.16 0.54 4.07 2.51 7.23 7.64 2.14 0.12 6.28 6.51 9.28 2.73 1.05 1.68 0.23 1.66 0.09];
    for i = 1:length(w); avg_lf{type,1}(i,1) = mean(let_freq(w(i,:))); end
end
xword = cell2mat(wordfeat); xnword = cell2mat(nwordfeat);

%% Main text 
acc_words = mean(L2_str.PC([1:100 201:350 501:700],:));
acc_nwords = mean(L2_str.PC([101:200 351:500 701:900],:));
mean_acc = [mean(acc_words) mean(acc_nwords)];
std_acc  = [std(acc_words) std(acc_nwords)];
nancorrcoef(nanmean(RT(qwords,1:2:end),2),nanmean(RT(qwords,2:2:end),2));
nancorrcoef(nanmean(RT(qnwords,1:2:end),2),nanmean(RT(qnwords,2:2:end),2));

%% Fig 4A and 4C
NWRT = RT(qnwords,:); WRT = RT(qwords,:);
% Storing word ids
NWidx = NaN(450,6); NWidx(1:100,1:4) = L2_str.word_id{1}(101:200,:);
NWidx(101:250,1:5) = L2_str.word_id{2}(151:300,:); NWidx(251:450,1:6) = L2_str.word_id{3}(201:400,:);
Widx = NaN(450,6); Widx(1:100,1:4) = L2_str.word_id{1}(1:100,:);
Widx(101:250,1:5) = L2_str.word_id{2}(1:150,:); Widx(251:450,1:6) = L2_str.word_id{3}(1:200,:);

% Sorting data in decreasing order of RT
[nwb,nwI] = sort(nanmean(NWRT,2),'descend');
[wb,wI]   = sort(nanmean(WRT,2),'descend');

cnt = 1;
for i = 1:2:10
    NWindex(cnt) = (i-1)*45+8; Windex(cnt) = (i-1)*45+8;
    nwT{cnt} = char(NWidx(nwI(NWindex(cnt)),:)+64);
    wT{cnt}  = char(Widx(wI(Windex(cnt)),:)+64); cnt = cnt +1;
end

% Changing a few index for better illustration
NWindex(1) = 9; nwT{1} = char(NWidx(nwI(9),:)+64);
NWindex(6) = 440;  nwT{6} = char(NWidx(nwI(440),:)+64);
Windex(6) = 440;  wT{6} = char(Widx(wI(440),:)+64);

% Figure 4A
figure; shadedErrorBar([],nanmean(WRT(wI,:),2), nansem(WRT(nwI,:),2),'r'); hold on;
plot(Windex,0.5*ones(6,1),'.r'); ylim([.4 1]); axis square; ylabel('Response time, s');

% Figure 4C
figure; shadedErrorBar([],nanmean(NWRT(nwI,:),2), nansem(NWRT(nwI,:),2)); hold on; 
plot(NWindex,0.5*ones(6,1),'.k'); ylim([.4 1]); axis square; ylabel('Response time, s');


%% Figure 4B
wordfeat1 = cell2mat(wordfeat); avg_lf1 = cell2mat(avg_lf);
wRT1 = cell2mat(wRT);

for rep = 1:1000
    q1 = randperm(size(wRT1,2), round(size(wRT1,2)/2)); q2 = setdiff(1:size(wRT1,2), q1);
    rt1 = nanmean(wRT1(:,q1),2); rt2 = nanmean(wRT1(:,q2),2);
    
    MAT = [wordfeat1(:,2) ones(length(rt1),1)]; W = regress(rt1, MAT); Rc(1,rep) = nancorrcoef(MAT*W, rt2);
    MAT = [wordfeat1(:,3) ones(length(rt1),1)]; W = regress(rt1, MAT); Rc(2,rep) = nancorrcoef(MAT*W, rt2);
    MAT = [log(wordfeat1(:,8)) ones(length(rt1),1)]; W = regress(rt1, MAT); Rc(3,rep) = nancorrcoef(MAT*W, rt2);
    MAT = [log(avg_lf1) ones(length(rt1),1)];        W = regress(rt1, MAT); Rc(4,rep) = nancorrcoef(MAT*W, rt2);
    MAT = [wordfeat1(:,[2 3]) log(wordfeat1(:,8)) log(avg_lf1) ones(length(rt1),1)]; W = regress(rt1, MAT); Rc(5,rep) = nancorrcoef(MAT*W, rt2);
    SH(rep) = nancorrcoef(rt1,rt2);
end
figure;
barweb(mean(Rc,2),std(Rc,[],2)); hold all; ylim([-0.1 0.6]);
shadedErrorBar([.7; 1.3],[mean(SH); mean(SH)], std(SH));  ylabel('Correlation Coefficient');
% [r p ] = partialcorri(yobs,MAT)

%% Figure 4E
for type = 1:3
    clear mat allparts Ypred
    strlen = numel(widx{type,1}(1,:));
    for L = 1:size(wRT{type,1},1)
        allparts(L,:) = [widx{type,1}(L,:) widx{type,1}(L+(50*(1+type)),:)];
    end
    rt = 1./nanmean(nwRT{type},2);
    w0 = rand(ncells*strlen+1,1);
    X = {allparts,rates};
    west = nlinfit(X,rt,@neuralmodel,w0,opts);
    Ypred = neuralmodel(west,X);
    NWrt_pred{type} = Ypred;
    NWrt_obs{type} = rt;
end

x1 = cell2mat(NWrt_pred'); x2 = cell2mat(NWrt_obs'); 
figure; corrplot(x1,x2,[],1); xlim([.8 2]); ylim([.8 2]); hold on; plot(x1(nwI(NWindex)),x2(nwI(NWindex)),'ro')
xlabel('Predicted 1/RT'); ylabel('Observed 1/RT')

%% Model fit using visual search weights
for type = 1:3
    clear mat allparts Ypred
    strlen = numel(widx{type,1}(1,:));
    for L = 1:size(wRT{type,1},1)
        allparts(L,:) = [widx{type,1}(L,:) widx{type,1}(L+(50*(1+type)),:)];
    end
    rt = 1./nanmean(nwRT{type},2);
    pred_dis = gendiss_neural(widx{type,1},allparts);
    NWrt_pred{type} = pred_dis;
    NWrt_obs{type} = rt;
end

x1 = cell2mat(NWrt_pred'); x2 = cell2mat(NWrt_obs'); 
figure; corrplot(x1,x2); 

%% Figure 4F
% Predicting on full dataset
for type = 1:3
    clear mat allparts Ypred
    strlen = numel(widx{type,1}(1,:));
    for L = 1:size(wRT{type,1},1)
        allparts(L,:) = [widx{type,1}(L,:) widx{type,1}(L+(50*(1+type)),:)];
    end
    rt = 1./nanmean(nwRT{type},2);
    w0 = rand(ncells*strlen+1,1);
    X = {allparts,rates};
    west = nlinfit(X,rt,@neuralmodel,w0,opts);
    Ypred = neuralmodel(west,X);
    nwordfeat{type}(nwordfeat{type}(:,3) == 0,3) = nan;
    MAT = [Ypred wordfeat{type}(:,2) log(wordfeat{type}(:,8)) nwordfeat{type}(:,1) log(nwordfeat{type}(:,3)) log(avg_lf{type}) ones(length(rt),1)]; W = regress(rt, MAT);
    NWrt_pred{type} = MAT*W; NWrt_pred{type} = 1./NWrt_pred{type};
    
    % Predicting word RTs
    rt = 1./nanmean(wRT{type},2);
    MAT = [wordfeat{type}(:,[2 3]) log(wordfeat{type}(:,8)) log(avg_lf{type}) ones(length(rt),1)]; W = regress(rt, MAT);
    Wrt_pred{type} = MAT*W;   Wrt_pred{type} = 1./Wrt_pred{type};
end

types{1,1} = 1:15;      types{1,2} = 1:15;      types{1,3} = 1:20;
types{2,1} = 16:30;     types{2,2} = 16:30;     types{2,3} = 21:40;
types{3,1} = [];        types{3,2} = 31:50;     types{3,3} = 41:70;
types{4,1} = [];        types{4,2} = 51:70;     types{4,3} = 71:100;
types{5,1} = 31:55;     types{5,2} = 71:105;    types{5,3} = 101:140;
types{6,1} = 56:70;     types{6,2} = 106:120;   types{6,3} = 141:160;
types{7,1} = 71:85;     types{7,2} = 121:135;   types{7,3} = 161:180;
types{8,1} = 86:100;    types{8,2} = 136:150;   types{8,3} = 181:200;

rawRT = cell(4,1); predRT = cell(4,1); ArawRT =cell(4,1);

cnt = 1;
for stimT = [2 1 6 7] % These 4 are the conditions of interest
    for L = 1:3
        ArawRT{cnt} = [ArawRT{cnt}; (nwRT{L}(types{stimT,L},:) - wRT{L}(types{stimT,L},:))./wRT{L}(types{stimT,L},:)]; % USED IN STATS

        % Data averaged across subjects
        rawRT{cnt} = [rawRT{cnt}; (nanmean(nwRT{L}(types{stimT,L},:),2) - nanmean(wRT{L}(types{stimT,L},:),2))./nanmean(wRT{L}(types{stimT,L},:),2)];
        predRT{cnt} = [predRT{cnt};  (NWrt_pred{L}(types{stimT,L},:)- Wrt_pred{L}(types{stimT,L},:))./ Wrt_pred{L}(types{stimT,L},:)];
    end
    cnt = cnt + 1;
end

avg_rt = cellfun(@nanmean, rawRT); sem_rt = cellfun(@nansem, rawRT);
avg_prt = cellfun(@nanmean, predRT); sem_prt = cellfun(@nansem, predRT);

figure; barweb([avg_rt avg_prt]',[sem_rt sem_prt]'); set(gca,'Xtick',1:3,'XTickLabel',{'Observed','Predicted'})
legend('Middle trans','Edge trans','Edge Subs','Middle subs'); ylabel('% change in RT'); ylim([0 0.2])

% STATS
cnt = 1; clear anova_rt subid condid
for cond = [1 2] % Change the condition number to perform appropriate comparisons. 
    for sub = 1:16
        for trial = 1:50
            anova_rt(cnt,1) =  ArawRT{cond}(trial,sub);
            subid(cnt,1) = sub;
            condid(cnt,1) = cond;
            cnt = cnt + 1;
        end
    end
end
p = anovan(anova_rt,{subid, condid},'full');
% p = anovan(anova_rt,{condid})

%% Figure 4G
cnt=1;
for type = 1:3
    for i = 1:size(widx{type},1)/2
        dist(cnt,1) = editDistance(char(widx{type}(i,:)+64),char(widx{type}(i+size(widx{type},1)/2,:)+64),'SubstituteCost',2); % OLD measure
        dpsy(cnt,1) = 1./nanmean(nwRT{type}(i,:),2);
        cnt = cnt+1;
    end
end

figure; corrplot(dist,dpsy); xlim([1 9]); hold on; plot(dist(nwI(NWindex)),dpsy(nwI(NWindex)),'ro')
xlabel('Orthographic Levenshtein Distance'); ylabel('Observed 1/RT' )

%% figure 4H
for rep = 1:1000
    q1 = randperm(size(nwRT{1},2), round(size(nwRT{1},2)/2)); q2 = setdiff(1:size(nwRT{1},2), q1);
    for type = 1:3
        % Modelling RTs as a linear combination of dissimilarity and lexical features
        clear mat allparts
        strlen = numel(widx{type,1}(1,:));
        for L = 1:size(wRT{type,1},1)
            allparts(L,:) = [widx{type,1}(L,:) widx{type,1}(L+(50*(1+type)),:)];
        end
        % dissimilarity
        w0 = rand(ncells*strlen+1,1);
        X = {allparts,rates};
        west = nlinfit(X,1./nanmean(nwRT{type,1}(:,q1),2),@neuralmodel,w0,opts);
        Ypred{type,1} = neuralmodel(west,X);
        rt1{type,1} = 1./nanmean(nwRT{type,1}(:,q1),2); rt2{type,1} = 1./nanmean(nwRT{type,1}(:,q2),2);
    end
    
    [Rc(1,rep), Pc(1,rep)] = nancorrcoef(cell2mat(Ypred), cell2mat(rt2));
    
    % lexical factors
    xnword(xnword(:,3) == 0,3) = nan;
    MAT = [xword(:,2) log(xword(:,8)) xnword(:,1) log(xnword(:,3)) log(cell2mat(avg_lf)) ones(size(xword,1),1)];  
    [W, Wint] = regress(cell2mat(rt1), MAT);
    [Rc(2,rep), Pc(2,rep)] = nancorrcoef(MAT*W, cell2mat(rt2));
    
    % Combined
    MAT = [cell2mat(Ypred) xword(:,2) log(xword(:,8)) xnword(:,1) log(xnword(:,3)) log(cell2mat(avg_lf)) ones(size(xword,1),1)];     W = regress(cell2mat(rt1), MAT);
    [Rc(3,rep), Pc(3,rep)] = nancorrcoef(MAT*W, cell2mat(rt2));
    
    % OLD distance - This bar is shifted to second position while generating the figure
    [Rc(4,rep), Pc(4,rep)] = nancorrcoef(dist, cell2mat(rt2)); 
    SH(rep) = nancorrcoef(cell2mat(rt1),cell2mat(rt2));
    
    disp(rep)
end

figure;
barweb(mean(Rc,2),std(Rc,[],2)); hold all; ylim([0 0.75]);
shadedErrorBar([.7; 1.3],[mean(SH); mean(SH)], std(SH));  ylabel('Correlation Coefficient');


%% non-word 1/RT prediction using SID model (text)
clear Rc Pc Yobs Ypred  
for rep = 1:30
    q1 = randperm(size(nwRT{1},2), round(size(nwRT{1},2)/2)); q2 = setdiff(1:size(nwRT{1},2), q1);
    for type = 1:3
        clear mat
        for L = 1:size(wRT{type},1)
            strlen = numel(widx{type}(1,:));
            [~,~,~,mat(L,:)] = Batonmodel([], 26, strlen,3,1, [widx{type}(L,:); widx{type}(L+(50*(1+type)),:)],[],dis_uletter);
        end
        opts = statset('nlinfit');
        opts.RobustWgtFun = 'bisquare';
        opts.MaxIter = 10000;
        b0 = [.8 -.1 .5 -2 1];
        %     b0 = [rand(1,5)];
        regmat = mat(:,varorder(strlen));
        [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(regmat,1./nanmean(nwRT{type}(:,q1),2),@batonfun,b0,opts);
        [Ypred{type,1},delta] = nlpredci(@batonfun,regmat,beta,R,'Covar',CovB);
        Yobs{type,1} = 1./nanmean(nwRT{type}(:,q2),2);
    end
    [Rc(rep), Pc(rep)] = nancorrcoef(cell2mat(Ypred), cell2mat(Yobs));
end