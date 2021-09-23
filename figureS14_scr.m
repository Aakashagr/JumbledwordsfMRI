allclear
% fitting a non-linear model and comparing coefficients with baton model.
global strlen ncells rates
load L2_letter
L2_str_letter = L2_str; clear L2_str
RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);

% get the single cell shape tuning required
D = squareform(dis_uletter); ncells = 10;
opts.MaxIter = 400;
[rates,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses

opts = statset('nlinfit');

%%
load L2_scrwrd; 

% Data processing
L2_str.RT(L2_str.RT < .3) = NaN; % accidental button press
L2_str.typed_word = cellfun(@lower, L2_str.typed_word,'UniformOutput',0);

% Wrongly decoded words
L2_str.typed_word(cellfun(@isempty, L2_str.typed_word)) = {'1'}; % filling up empty cells with nonsensical number
for wrd = 1:300
    L2_str.RT(wrd,~(ismember(L2_str.typed_word(wrd,:),[char(L2_str.word_id{wrd} + 96)]))) = NaN;
    nNan(wrd,1) = sum(isnan(L2_str.RT(wrd,:)));
end

% Extracting parameters
RT = L2_str.RT;
words = L2_str.word_id;
scrwords = L2_str.scrword_id;

% Excluding noisy trials
nsub = 10;
RT(nNan > nsub,:) = []; words(nNan > nsub) = []; scrwords(nNan > nsub) = [];  L2_str.images(nNan > nsub) = [];
Wlen = cellfun(@length,words);

% consistency across subjects
nancorrcoef(nanmean(RT(:,1:2:end),2),nanmean(RT(:,2:2:end),2));
rdata = spearmanbrowncorrection(splithalfcorrd(RT'),2);

% Accuracy
meanacc = 1-mean(sum(isnan(RT))/size(RT,1));
%% Fig S14A. Identifying easy and hard trials
[b,I] = sort(nanmean(RT,2),'descend');
cnt = 1;
for i = 1:2:10
    index(cnt) = (i-1)*24 +6;
    T{cnt} = char(scrwords{I(index(cnt))}+64);
    t(cnt) = b(index(cnt)); cnt = cnt+1;
end
index(6) = 237;  T{6} = char(scrwords{I(index(6))}+64);
figure; shadedErrorBar([],nanmean(RT(I,:),2), nansem(RT(I,:),2)); hold on
plot( index,0.5*ones(6,1), 'r.');
ylabel('Unscrambling time (s)'); xlabel('Scrambled words');
%% Fig S14C
for strlen = 4:6
    clear diss mat allparts
    idx = find(Wlen == strlen);
    for i = 1:numel(idx)
        allparts(i,:) = [words{idx(i)} scrwords{idx(i)}]; 
        diss(i,:) = RT(idx(i),:);
    end
    w0 = ones(ncells*strlen+1,1);
    X = {allparts,rates};    
    west = nlinfit(X,nanmean(diss,2),@neuralmodel,w0,opts);
    Ypred{strlen-3} = neuralmodel(west,X);        
    rt1{strlen-3} = nanmean(diss,2);  rt2{strlen-3} = nanmean(diss,2);
end
obs_rt2 = [rt2{1}; rt2{2}; rt2{3}];
pred_dis = [Ypred{1}; Ypred{2}; Ypred{3}];

figure; corrplot(pred_dis, obs_rt2,[],1); hold on; plot(pred_dis(I(index)), obs_rt2(I(index)),'ro'); ylim([.5 4.5]);xlim([.5 4.5])
% seescatterplot(pred_dis, obs_rt2,L2_str.images)
%% Fig S14D
% OLD measure
for i = 1:238, dist(i) = editDistance(char(words{i}+64),char(scrwords{i}+64),'SubstituteCost',2); end
figure; corrplot(dist,nanmean(RT,2)); hold on; plot(dist(I(index)),nanmean(RT(I(index),:),2),'ro'); xlim([1 9])


%% Fig S14D
for rep = 1:1000
    q1 = randperm(size(RT,2), round(size(RT,2)/2)); q2 = setdiff(1:size(RT,2), q1);
    
    for strlen = 4:6
        clear diss mat allparts
        idx = find(Wlen == strlen);
        for i = 1:numel(idx)
            allparts(i,:) = [words{idx(i)} scrwords{idx(i)}]; 
            diss(i,:) = RT(idx(i),:);
        end
        w0 = ones(ncells*strlen+1,1);
        X = {allparts,rates};    
        west = nlinfit(X,nanmean(diss(:,q1),2),@neuralmodel,w0,opts);
        Ypred{strlen-3} = neuralmodel(west,X);        
        rt1{strlen-3} = nanmean(diss(:,q1),2);  rt2{strlen-3} = nanmean(diss(:,q2),2);
    end
    obs_rt2 = [rt2{1}; rt2{2}; rt2{3}];
    obs_rt1 = [rt1{1}; rt1{2}; rt1{3}];
    pred_dis = [Ypred{1}; Ypred{2}; Ypred{3}];

    % Adding lexical features
    % Mean letter frequency
    let_freq = [8.04 1.48 3.34 3.82 12.49 2.4 1.87 5.05 7.57 0.16 0.54 4.07 2.51 7.23 7.64 2.14 0.12 6.28 6.51 9.28 2.73 1.05 1.68 0.23 1.66 0.09];
    for i = 1:numel(words); avg_lf(i,1) = mean(let_freq(words{i})); end
    
    % Non-Word features
    [num, str] = xlsread('feat_nw_scr.xlsx'); clear varmat
    for i = 1:length(scrwords); varmat(i,:) = num(ismember(str(2:end,2), char(scrwords{i}+64)),4:7); end
    varmat(isnan(varmat)) = 0; nwordfeat = varmat;
    nwordfeat(nwordfeat(:,3) == 0,3) = nan;

    % Word features
    [num, str] = xlsread('feat_w_ldt.xlsx'); clear varmat
    for i = 1:length(words); varmat(i,:) = num(ismember(str(2:end,1), char(words{i}+64)),2:12); end
    varmat(isnan(varmat)) = 0; wordfeat = varmat;  old20 = varmat(:,6);
    
    Rc(1,rep) = nancorrcoef(pred_dis, obs_rt2);

    % lexical factors
    MAT = [wordfeat(:,2) log(wordfeat(:,8)) nwordfeat(:,1) log(nwordfeat(:,3)) log(avg_lf) ones(length(obs_rt2),1)];  W = regress(obs_rt1, MAT); lex_pred = MAT*W;
    Rc(2,rep) = nancorrcoef(lex_pred, obs_rt2);
    
    % Combined
    MAT = [pred_dis wordfeat(:,2) log(wordfeat(:,8)) nwordfeat(:,1) log(nwordfeat(:,3)) log(avg_lf) ones(length(obs_rt2),1)]; W = regress(obs_rt1, MAT);
    Rc(3,rep) = nancorrcoef(MAT*W, obs_rt2);
    
    % OL distance
    Rc(4,rep) = nancorrcoef(dist, obs_rt2);

    SH(rep) = nancorrcoef(obs_rt2,obs_rt1);
    disp(rep)
end
figure;
barweb(mean(Rc,2),std(Rc,[],2)); hold all; ylim([0 0.6]);
shadedErrorBar([.7; 1.3],[mean(SH); mean(SH)], std(SH));  ylabel('Correlation Coefficient');

%% Text data 
for strlen = 4:6
    clear diss mat allparts
    idx = find(Wlen == strlen);
    for i = 1:numel(idx)
        allparts(i,:) = [words{idx(i)} scrwords{idx(i)}]; 
        diss(i,:) = RT(idx(i),:);
    end
    w0 = ones(ncells*strlen+1,1);
    X = {allparts,rates};    
    west = nlinfit(X,nanmean(diss,2),@neuralmodel,w0,opts);
    Ypred{strlen-3} = neuralmodel(west,X);        
    rt1{strlen-3} = nanmean(diss,2);  rt2{strlen-3} = nanmean(diss,2);
end
obs_rt = [rt2{1}; rt2{2}; rt2{3}];
pred_dis = [Ypred{1}; Ypred{2}; Ypred{3}];
MAT = [wordfeat(:,2) log(wordfeat(:,8)) nwordfeat(:,1) log(nwordfeat(:,3)) log(avg_lf) ones(length(obs_rt),1)];  W = regress(obs_rt1, MAT); lex_pred = MAT*W;


[rl, pl ] = partialcorri(obs_rt,MAT); % for lexical model
[rm, pm ] = partialcorri(obs_rt,[lex_pred pred_dis]);

