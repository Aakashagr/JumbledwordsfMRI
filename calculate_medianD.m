allclear
% The dissimilarity mat files are generated using searchlight_rsa_gendis.m
f = dir('*.mat');

for i = 1:numel(f)
    load(f(i).name);
    Xall{i} = droi;
end
%%
dmedian = single(nan(53,63,52,2701));
for x = 1:53
    for y = 1:63
        for z = 1:52
            clear dtemp
            for n = 1:17
                dtemp(n,:) = Xall{n}(x,y,z,:);
            end
            dmedian(x,y,z,:) = nanmedian(dtemp);
        end
    end
    disp(x)
end


save median_droi dmedian 