% Searchlight on dissimilarity using RSA toolbox
% allclear; load L2fmri_LDTw
allclearL2

conditions = repmat([1:74 zeros(1,7)], [1,8]); conditions = conditions';    
partition = [ones(81,1);ones(81,1)*2;ones(81,1)*3;ones(81,1)*4;ones(81,1)*5;ones(81,1)*6;ones(81,1)*7;ones(81,1)*8];
nsubs = numel(L2_str.subjectid);

for sub = 1:nsubs
    beta = reshape(L2_str.mergedevtbeta{sub},[53 63 52 75]);  beta(:,:,:,75) = [];
    mask = reshape(L2_str.mergedevtbetash{sub},[53 63 52]);   mask(~isnan(mask)) = 1;
    [x,y,z] = ind2sub([53 63 52],find(mask >0));
    
    droi = single(nan(53,63,52,2701));
    
    % Trial level data for RSA toolbox
    fspm = dir(['..\preprocessing\*', L2_str.subjectid{sub}]);
    load(['..\preprocessing\', fspm.name,'\glm\wdnevt\SPM.mat']);
    f = dir(['..\preprocessing\',fspm.name,'\nii\dwars*_evt_*.nii']);

    cnt = 1; actRDM = single(nan(53, 63, 52,1504));
    for run = 1:8
        rawact = spm_read_vols(spm_vol(['..\preprocessing\', fspm.name,'\nii\' f(run).name]));
        for sid = 1:size(rawact,4)
            actRDM(:,:,:,cnt)   = rawact(:,:,:,sid);
            cnt = cnt + 1;
        end
    end
    
    for idx = 1:numel(x)
        x_range = max(x(idx)-2,1):min(x(idx)+2,53);
        y_range = max(y(idx)-2,1):min(y(idx)+2,63);
        z_range = max(z(idx)-2,1):min(z(idx)+2,52);
        
        % Identifying regions inside the brain
        act = reshape(beta(x_range,y_range,z_range,:),[numel(x_range)*numel(y_range)*numel(z_range),74]);
        actr = reshape(actRDM(x_range,y_range,z_range,:),[numel(x_range)*numel(y_range)*numel(z_range),size(actRDM,4)]);
        actr(isnan(nanmean(act,2)) ,:) = [];  % removing voxels which are not included in the mask.
             
        if size(actr,1) > 30
            betaRDM = rsa.spm.noiseNormalizeBeta(double(actr'),SPM);    
            droi(x(idx),y(idx),z(idx),:) = rsa.distanceLDC(betaRDM,partition,conditions,SPM.xX.xKXs.X);
        end
        disp(['X :', num2str(idx)])
    end
    
    save(['droi_sub',num2str(sub,'%02d')],'droi')
    disp(sub)
end
