function [sh, Rn] = partmodel_neu(imgparts, RT,image_pairs)
global strlen rates ncells
if ~exist('image_pairs'); image_pairs = nchoosek(1:size(imgparts,1),2); end
for S = 1:100
    q1 = randperm(size(RT,2), round(size(RT,2)/2)); q2 = setdiff(1:size(RT,2), q1);
    dis1 = 1./nanmean(nanmean(RT(:,q1,:),3),2); dis2 = 1./nanmean(nanmean(RT(:,q2,:),3),2);
    sh(1,S) = nancorrcoef(dis1,dis2);   
    
    opts = statset('nlinfit'); opts.UseParallel = true;
    w0 = rand(ncells*strlen+1,1);
    allparts = [imgparts(image_pairs(:,1),:) imgparts(image_pairs(:,2),:)]; 
    X = {allparts,rates};    
    west = nlinfit(X,dis1,@neuralmodel,w0,opts);
    dpred = neuralmodel(west,X);
    
    Rn(1,S) = nancorrcoef(dpred,dis2);
    disp(S)
end
