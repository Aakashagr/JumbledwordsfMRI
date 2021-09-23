% function RT = rmRToutlier(rt,rt_thresh, ncounts)
% This function removes outliers in Reaction time, If the frequency of
% Reaction time at given thershold is less than ncounts. Then those RTs are
% set to nan. 
% Input:
%       rt = Input raw reaction time
%       rt_thresh = Cut-off threshold of reaction time
%       ncounts = minimum number of occurrence of RT above the mentioned threshold
% Output:
%       RT = Denoised Reaction Time 

function RT = rmRToutlier(rt,rt_thresh, ncounts)

rt(rt < .3) = NaN;  % removing accidental key press.
mat = zeros(size(rt)); mat(rt>rt_thresh) = 1; % identifying harder search pairs
temp = sum(sum(mat,2),3);
% figure; stem(temp) % identifying mistakes.
id = find(temp > 0 & temp < ncounts);

% removing higher RT outlier
for i = 1:length(id)
    [c, v] =ind2sub([size(rt,2),2],find(mat(id(i),:,:)==1));
    rt(id(i),c,v) = NaN;
end
RT = rt;