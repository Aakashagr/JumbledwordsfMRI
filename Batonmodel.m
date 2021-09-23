% Batonmodel  --> Fit Baton model for stimuli of any length
% [dpred,b,bint,X] = Batonmodel(dobs, npart, stimlen,modeltype,sepflag,imgparts,srcpair,pwdiss)
% Required inputs
%    npart         = Number of unique image parts.
%    stimlen       = Number of parts used in each stimulus
%    Model type:   1 - Full model, parameters for absolute position
%                  2 - Reduced model, parameters for relative position
%                  3 - Reduced model-1, uses letter pairs dissimilarities and to estimate scaled model (type 1) coefficients
%                  4 - Reduced model-2, uses letter pairs dissimilarities and to estimate scaled model (type 2) coefficient
% Optional inputs:
%    imgparts      = part ids of each stimulus used in the experiment. (Default: all possible part combinations)
%    srcpair       = search pair ids used in the experiment. (Default: all possible searches i.e. nchoosek(1:imgparts,2))
%    dobs          = npairs x 1 vector of observed dissimilarities
%    pwdiss        = pair-wise part dissimilarity, Used only in Model type 3 and 4
% Outputs:
%    dpred         = predicted dissimilarities from the baton model
%    b             = estimated part distances
%    X             = npairs x nparameters regression matrix


function [dpred,b,bint,X] = Batonmodel(dobs, npart, stimlen,modeltype,sepflag,imgparts,srcpair,pwdiss)

if ~exist('imgparts') || isempty(imgparts), imgparts = permn(1:npart,stimlen); end  % default: all possible arrangement with replacement
if ~exist('srcpair')  || isempty(srcpair),  srcpair = nchoosek(1:size(imgparts,1),2); end % default: all possible search pairs

partpairs = nchoosek(1:stimlen,2); % all possible part pairs
[~,qid] = sort(partpairs(:,2) - partpairs(:,1)); % sorting part pairs to group N
partpairs = partpairs(qid,:);

% Generating index of all position pairs
pospairs = nchoosek(1:2*stimlen,2); % all position pairs between target and distractor
if modeltype == 1 ||  modeltype == 3  % Full model
    parm_len = stimlen^2; % total number of parameters
    for P = 1:size(pospairs,1)
        qt = pospairs(P,:);
        qdiff = abs(qt(1) - qt(2));
        
        % within parts
        if (qt(1)<stimlen+1 && qt(2) < stimlen+1)
            partpos(P) = stimlen + nchoosek(stimlen,2) + find(ismember(partpairs,sort(qt),'rows'));
        elseif (qt(1) > stimlen && qt(2) > stimlen)
            partpos(P) = stimlen + nchoosek(stimlen,2) + find(ismember(partpairs,sort(qt)-stimlen,'rows'));
            
            % corresponding parts
        elseif qdiff == stimlen
            partpos(P) = min(qt);
            
            % across parts
        elseif (qt(1)<stimlen+1 && qt(2) > stimlen) &&  qdiff ~= stimlen
            qt(2) = qt(2) - stimlen;
            partpos(P) = stimlen + find(ismember(partpairs,sort(qt),'rows'));
        else
            qt(1) = qt(1) - stimlen;
            partpos(P) = stimlen + find(ismember(partpairs,sort(qt),'rows'));
        end
    end
    
elseif  modeltype == 2 ||  modeltype == 4
    % reduced model
    parm_len = stimlen + (stimlen - 1)*2;
    for P = 1:size(pospairs,1)
        qt = pospairs(P,:);
        qdiff = abs(qt(1) - qt(2));
        
        % within parts
        if (qt(1)<stimlen+1 && qt(2) < stimlen+1) || (qt(1) > stimlen && qt(2) > stimlen)
            partpos(P) = stimlen + stimlen - 1 + qdiff;  % ncorr + n-1 acorss + qwithin
            
            % corresponding parts
        elseif qdiff == stimlen, partpos(P) = min(qt);
            
            % Across parts
        else,  partpos(P) = stimlen + abs(stimlen - qdiff);  end  % ncorr + q acorss
    end
    
end

% Building up the design matrix
part_pairs = nchoosek(1:npart,2);
if modeltype < 3
    X = zeros(size(srcpair,1),size(part_pairs,1)*parm_len+1);
    X(:,end) = 1;
    xx = nchoosek(npart,2);
    
    for  i = 1:size(srcpair,1)
        q1 = [imgparts(srcpair(i,1),:), imgparts(srcpair(i,2),:)];
        for  pid = 1:numel(partpos)
            [~,id] = ismember(sort([q1(pospairs(pid,1)), q1(pospairs(pid,2))]), part_pairs,'rows');
            if id>0; X(i,xx*(partpos(pid)-1)+ id)  = X(i,xx*(partpos(pid)-1)+ id) + 1; end
        end
    end
else
    % using part position from full model
    X = zeros(size(srcpair,1),parm_len+1);
    X(:,end) = 1;
    for  i = 1:size(srcpair,1)
        q1 = [imgparts(srcpair(i,1),:), imgparts(srcpair(i,2),:)];
        for  pid = 1:numel(partpos)
            [~,id] = ismember(sort([q1(pospairs(pid,1)), q1(pospairs(pid,2))]), part_pairs,'rows');
            if id>0; X(i,partpos(pid))  = X(i,partpos(pid)) + pwdiss(id); end
        end
    end
end

if(~sepflag)
    if modeltype < 3
        temp = sum(reshape(X(:,1:size(part_pairs,1)*stimlen),[size(X,1) size(part_pairs,1) stimlen]),3);
        temp = [temp, X(:, size(part_pairs,1)*stimlen +1:end)];
        X = temp;
    else
        temp = sum(X(:,1:stimlen),2);
        temp = [temp X(:,stimlen+1 :end)];
        X = temp;
    end
end

% predicting dissimilarities
if isempty(dobs)
    dpred = []; b = []; bint = [];
else
    [b, bint] = regress(dobs, X); dpred = X*b;
end
