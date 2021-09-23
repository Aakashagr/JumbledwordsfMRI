% allclear; load L2fmri_LDTw

%%
allclearL2
BETAs = L2_str.mergedevtbeta;
qw = 1:32; qn = 33:64; ql = 65:74;
nsubs = numel(L2_str.subjectid);
[idx, ROIname] = getvoxind(L2_str);
%%
if 0
for sub = 1:nsubs
    fspm = dir(['..\preprocessing\*', L2_str.subjectid{sub}]);
    load(['..\preprocessing\', fspm.name,'\glm\wdnevt\SPM.mat']);
    f = dir(['..\preprocessing\',fspm.name,'\nii\dwars*_evt_*.nii']);
    
    cnt = 1; act_EVC = [];act_V4 = [];act_LOC = [];act_VWFA = [];act_TG = [];
    for run = 1:8
        rawact = spm_read_vols(spm_vol(['..\preprocessing\', fspm.name,'\nii\' f(run).name]));
        for sid = 1:size(rawact,4)
            temp = rawact(:,:,:,sid);
            act_EVC(cnt,:)   = temp(idx{sub,1});
            act_V4(cnt,:)    = temp(idx{sub,2});
            act_LOC(cnt,:)   = temp(idx{sub,3});
            act_VWFA(cnt,:)  = temp(idx{sub,4});
            act_TG(cnt,:)    = temp(idx{sub,5});
            cnt = cnt + 1;
        end
    end
    %     conditions = 1:74; conditions = conditions';
    conditions = repmat([1:74 zeros(1,7)], [1,8]); conditions = conditions';    
    partition = [ones(81,1);ones(81,1)*2;ones(81,1)*3;ones(81,1)*4;ones(81,1)*5;ones(81,1)*6;ones(81,1)*7;ones(81,1)*8];
    
    beta = rsa.spm.noiseNormalizeBeta(act_EVC,SPM);     droi(:,sub,1) = rsa.distanceLDC(beta,partition,conditions,SPM.xX.xKXs.X);
    beta = rsa.spm.noiseNormalizeBeta(act_V4,SPM);      droi(:,sub,2) = rsa.distanceLDC(beta,partition,conditions,SPM.xX.xKXs.X);
    beta = rsa.spm.noiseNormalizeBeta(act_LOC,SPM);     droi(:,sub,3) = rsa.distanceLDC(beta,partition,conditions,SPM.xX.xKXs.X);
    beta = rsa.spm.noiseNormalizeBeta(act_VWFA,SPM);    droi(:,sub,4) = rsa.distanceLDC(beta,partition,conditions,SPM.xX.xKXs.X);
    beta = rsa.spm.noiseNormalizeBeta(act_TG,SPM);      droi(:,sub,5) = rsa.distanceLDC(beta,partition,conditions,SPM.xX.xKXs.X);
    disp(sub)
end
end
% save dfmri_LDC droi
