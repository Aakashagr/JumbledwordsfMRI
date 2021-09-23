% function [RT frac] = rmRToutlier1(rt)
% This function removes outliers in Reaction time using a MATLAB function is outlier. Then those RTs are
% set to nan. 
% Input:
%       rt = Input raw reaction time
% Output:
%       RT = Denoised Reaction Time 
%       frac = fraction of data points set to NaN;

function [RT, frac] = rmRToutlier1(rt,r1,r2)

% rt(rt < .3) = NaN;  % removing accidental key press.
for i = 1:size(rt,1)
    tmp(i,:,:) = isoutlier(rt(i,:,:));
end

RT = rt;
RT(tmp == 1) = NaN;
frac = numel(find(isnan(RT(:))))/numel(RT);