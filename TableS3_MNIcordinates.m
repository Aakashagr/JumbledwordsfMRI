if ~exist('L2_str'); load L2fmri_LDTw; end
allclearL2
[ids, ROIname] = getvoxind(L2_str); % Extracting voxel IDs of various Region of Interest (ROI)-interesection of functional and anatomical regions
%% Extracting Peak ROI location

for roi = 1:5
    for sub = 1:17
        [x,y,z] = ind2sub([53 63 52], ids{sub,roi}(1));
        idx(roi,sub,:) = cor2mni([x,y,z]);
    end
end

std(idx(:,:,3),[],2)