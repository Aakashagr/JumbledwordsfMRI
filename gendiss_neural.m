function pred_dis = gendiss_neural(stimuli,allparts)


% allclear
% load L2_letter
% L2_str_letter = L2_str; clear L2_str
% RT_letter = rmRToutlier1(L2_str_letter.RT,1.7, 5);
% uppercase = find(ismember(L2_str_letter.img_pairs, nchoosek(1:26,2),'rows'));
% dis_uletter = 1./nanmean(nanmean(RT_letter(uppercase,:,:),3),2);
% D = squareform(dis_uletter); ncells = 10;
% opts.MaxIter = 400;
% [rates6,e] = mdscale(D,ncells,'Options',opts); % rates = nstim x ncells matrix of single letter responses
% 
% % 6 letter
% load L2cb_6letter; L2cb_6 = L2_str; clear L2_str
% acc = mean(L2cb_6.PC);
% RT6 = rmRToutlier1(L2cb_6.RT,3, 3); dobs6 = 1./nanmean(nanmean(RT6,3),2); 
% %%
% 
% dis6_1 = 1./nanmean(nanmean(RT6,3),2);  
% imgparts = L2cb_6.letterind; srchpairs = nchoosek(1:36,2); 
% allparts = [imgparts(srchpairs(:,1),:) imgparts(srchpairs(:,2),:)]; 
% X6 = {allparts,rates6}; nparts = size(imgparts,2);
% rng('default'); w0 = rand(ncells*nparts+1,1);
% west = nlinfit(X6,dis6_1,@neuralmodel,w0,opts); 

%%
% load W_gendiss_neural
load w6_mds
nstr = size(stimuli,2); clear B
for neu = 1:10
    q = (neu-1)*6+1:neu*6;
    f = fit(vec(1:6),west(q),'poly2');
    B(:,neu) = polyval([f.p1 f.p2 f.p3],linspace(1,6,nstr));
end

Wsamp = vec(B); Wsamp(end+1) = west(end);

srchpairs = nchoosek(1:size(stimuli,1),2); 
if ~exist('allparts'), allparts = [stimuli(srchpairs(:,1),:) stimuli(srchpairs(:,2),:)];  end
X = {allparts,rates6};
pred_dis = neuralmodel(Wsamp,X);