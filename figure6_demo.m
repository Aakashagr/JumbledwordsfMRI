%% figure 6A
%create MDS for a given set of words

allclear
word_id = ['FORGET'-64; 'FGROET'-64;'FGEORT'-64; 'FOGRET'-64; 'OFRGET'-64];
word_id = [word_id; 'FORCET'-64; 'FORXET'-64; 'FGORET'-64; 'FROGET'-64;[ 6 18 63 63 7 20];[20 7 63 63 18 6]];

load w6_mds
srchpairs = nchoosek(1:size(word_id,1),2);
allparts = [word_id(srchpairs(:,1),:) word_id(srchpairs(:,2),:)]; 
X6 = {allparts,rates6}; 
pred_dis = neuralmodel(west,X6); 


Y = mdscale(squareform(pred_dis),2);
for i = 1:size(Y,1)
    plot(Y(i,1),Y(i,2),'.');     
    text(Y(i,1),Y(i,2),char(word_id(i,:)+64));hold on
end
% figure; corrplot(pdist(Y),pred_dis)

%%   Figure 6B
% Creating sentences with increasing or decreasing order of reading difficulty
% HUMAN MIND DOES NOT READ EVERY LETTER BY ITSELF
% HUAMN MNID DEOS NOT RAED ERVEY LTETER BY ISTLEF

% allclear 
clearvars -except pred_dis

% % Sentence 1- with increasing difficulty 
words = {'HUMAN' 'MIND' 'DOES' 'READ' 'EVERY' 'LETTER' 'ITSELF'};
scr_stim = {'HUAMN' 'MIDN' 'DOSE' 'RDEA' 'YVERE' 'TRETLE' 'FSLTEI'};
for i = 1:numel(words)    
    pred_dis(i) = gendiss([words{i}-64; scr_stim{i}-64]);
end

% Sentence 2- with HLHLH
% words = {'HUMAN' 'MIND' 'DOES' 'READ' 'EVERY' 'LETTER' 'ITSELF'};
% scr_stim = {'HUAMN' 'DNMI' 'DEOS' 'DAER' 'EVREY' 'ETTELR' 'ITSLEF'};
% for i = 1:numel(words)   
%     pred_dis(i) = gendiss([words{i}-64; scr_stim{i}-64]);
% end

% Sentence 3- with decreasing difficulty
% words = {'HUMAN' 'MIND' 'DOES' 'READ' 'EVERY' 'LETTER' 'ITSELF'};
% scr_stim = {'ANMHU' 'DINM' 'SOED' 'RDEA' 'ERVEY' 'LTETER' 'ITSLEF'};
% for i = 1:numel(words)   
%     pred_dis(i) = gendiss([words{i}-64; scr_stim{i}-64]);
% end

%% Figure 6C

% stimuli = [18 5 1 4 9 14 7; 18 4 14 9 5 1 7; 18 5 1 15 54 14 7];  % READING
% stimuli = [18 5 3 5 14 20; 18 5 5 14 3 20; 18 56 3 56 14 20];  % RECENT
stimuli = [3 21 12 20 21 18 1 12; 3 12 18 20 21 1 21 12; 3 21 12 60 21  18 57 12];  % CULTURAL
% stimuli = [9 14 22 5 14 20 9 15 14; 9 15 14 5 20 14 9 22 14; 9 14 22 56 14 60 54 53 14];  % INVENTION

pred_dis = gendiss(stimuli)

