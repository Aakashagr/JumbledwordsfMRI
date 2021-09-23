function pred_dis = gendiss(stimuli,srchpairs)
% pred_dis = gendiss(stimuli)
% Calculates pair-wise dissimlarities for any pair of stimulus
% stimuli should be a vector comprising of numbers between 1-62
% 1:26 for A-Z, 27:52 for a-z, and 53:62 for 0-9

load L2_letter
RT_letter = rmRToutlier(L2_str.RT,1.7, 5);
uppercase = find(ismember(L2_str.img_pairs, nchoosek(1:26,2),'rows'));
dis_letter = 1./nanmean(nanmean(RT_letter,3),2);
if ~exist('srchpairs'), srchpairs = []; end
nstr = size(stimuli,2);
clear X
f = fit(vec(1:6),vec([0.1352 0.0681 0.0493 0.0372 0.0538 0.0670]),'poly2');
B = polyval([f.p1 f.p2 f.p3],linspace(1,6,nstr));
for i = 1:nstr
    X{i} = B(i)*exp(-1*[0 1:nstr-i]);
    X{i+nstr} = -X{i}(2:end);
end
W = cell2mat(X);
[~,~,~,mat] = Batonmodel([], 62, nstr,3,1, stimuli ,srchpairs,dis_letter);
pred_dis = mat(:,varorder(nstr))*vec(W);

