%% Generating RDMs

if ~exist('L2_str'); load L2fmri_LDTw; end

%% Figure 1
allclearL2
qw = 1:32; qn = 33:64; ql = 65:74;
RT = L2_str.RT;
nsubs = numel(L2_str.subjectid);
[idx, ROIname] = getvoxind(L2_str);

% Loading behaviour and fMRI dissimilarities
load dfmri_LDC; droi(isoutlier(droi,2)) = nan; % Removing outliers
dpsyw = L2_str.search.w;
dpsyn = L2_str.search.n;
dpsynw = L2_str.search.wnw;

% Loading semantic dissimilarity
load semanticfeat; dsemfeat = pdist(semanticfeat,'correlation')';


%% RDMs for behaviour and semantic data
wid = find(ismember(nchoosek(1:64,2),nchoosek(1:32,2),'rows'));
nid = find(ismember(nchoosek(1:64,2),nchoosek(1:32,2)+32,'rows'));
nwid = find(ismember(nchoosek(1:64,2),[1:32; 33:64]','rows'));
psyRDM = nan(nchoosek(64,2),1);
psyRDM([wid; nid; nwid]) = [dpsyw;dpsyn;dpsynw];
psyRDM = squareform(psyRDM); psyRDM(psyRDM ==0) = nan;
figure; heatmap((psyRDM)); colormap(jet); title('Perceptual space')

semRDM = nan(nchoosek(64,2),1);
semRDM(wid) = dsemfeat;
semRDM = squareform(semRDM); semRDM(semRDM ==0) = nan;
figure; heatmap((semRDM)); colormap(jet); title('Semantic space')

%% RDMs for fMRI data
sid = find(ismember(nchoosek(1:74,2),nchoosek(1:64,2),'rows'));
for roi = 1:5
    roiRDM = squareform(nanmedian(droi(sid,:,roi),2));
    roiRDM(roiRDM ==0) = nan;
    figure; heatmap((roiRDM)); colormap(jet); title(ROIname{roi}); 
end