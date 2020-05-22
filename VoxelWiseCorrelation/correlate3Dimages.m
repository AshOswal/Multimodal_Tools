function correlate3Dimages(Xi,Yi,design_matrix,prefix,mask,perm,avout)

% X is cell array of nifti images which will be the independent variable
% Y is cell array of nifti images which will be the dependent variable
% design_matrix for entry to spm_ancova
% nifti mask for cluster based permutation testing
% whether to meanscale - this is usually advised in this situation
% whether or not to perform permutation test for null distribution of
% cluster extent threshold

XX = [];
YY = [];

if perm
    Nperm = 1000;
else
    Nperm = 1;
end

if ~isempty(mask)
    mY = 0;
    for m = 1:numel(mask)
        vm = spm_vol(mask{m});
        [Ym,XYZm] = spm_read_vols(vm);
        mY = mY+Ym;
    end
    mind = find(mY>0);
end

shuffvoxels = 1;
for n = 1:numel(Xi)
    v0 = spm_vol(Xi{n});
    [X,XYZ] = spm_read_vols(v0);
    
    v1 = spm_vol(Yi{n});
    [Y,XYZ] = spm_read_vols(v1);
    
    if ~isempty(mask)
       X = X(mind);
       Y = Y(mind);      
    else
       mind = 1:numel(X); 
    end
    XX = [XX,X(:)];
    YY = [YY,Y(:)];
end
rmfield(v0,'pinfo')
if exist('avout','var')
    avX = mean(XX,2);
    avY = mean(YY,2);
    avX = avX > spm_percentile(avX,95);
    avY = avY > spm_percentile(avY,95);
    
    imX       = zeros(prod(v0.dim),1);
    imX(mind) = avX;
    
    imY       = zeros(prod(v0.dim),1);
    imY(mind) = avY;
    
    v0.fname = [prefix 'thresh_avX.nii'];
    spm_write_vol(v0,reshape(imX,v0.dim));
    v0.fname = [prefix 'thresh_avY.nii'];
    spm_write_vol(v0,reshape(imY,v0.dim));
end

p              = zeros(prod(v0.dim),Nperm);
r_squared      = zeros(prod(v0.dim),Nperm);
F_statistic    = zeros(prod(v0.dim),Nperm);
permclustcount = [];
clustcount     = [];
%%

for n = 1:Nperm
    disp(['permutation' num2str(n)])
    for i = 1:size(XX,1)
        
        if n == 1
            X = XX(i,:)';
            Y = YY(i,:)';
        else
            ind = randperm(size(XX(i,:),2));
            X = XX(i,ind)';
            Y = YY(i,randperm(numel(ind)))';
            % shuffle over voxels as well
            if shuffvoxels
               ind2= randperm(size(XX,1));
               X = XX(ind2(i),ind)';
            end            
        end
        
        % dont really need this
        if all(isnan(Y))||all(isnan(X))||sum(Y)==0||sum(X)==0
            continue
        end
      
        DM     = design_matrix;

        % GLM
        con = [1;zeros(size(DM,2),1)];
        [T,df,beta,~,c]=spm_ancova([X,DM],[],Y,con);
        
        F = T^2;
        prob=1-spm_Fcdf(F,df(1),df(2));
        % r sequared coefficient
        %TSS  =  sum((Y - mean(Y)).^2);
        %RSS  =  sum((Y - ([X,DM]*beta)).^2);
        
        
        if  prob<0.01 && beta(1)>0
            if c.c(1) ~= 1 
               keyboard;
            end
            p(mind(i),n) = 1;
            %r_squared(mind(i),n) = 1 - (RSS/TSS);
            F_statistic(mind(i),n) = F;
        end
        
    end
    % cluster count
    [L,num] = spm_bwlabel(reshape(p(:,n),v0.dim));
    for N = 1:num
        if isequal(n,1)
            L1   = L;
            clustcount = [clustcount, sum(L(:) == N)];
        else
            if sum(L(:) == N) > 1
            permclustcount = [permclustcount,sum(L(:) == N) ];
            end
        end
    end
end

if ~isempty(permclustcount)
    cind = clustcount>spm_percentile(permclustcount,99);
else
    cind = 1:numel(clustcount);
end

% extent_threshold = 100;
extent_threshold = 50;
cind = cind & (clustcount>extent_threshold);

L1        = ismember(L1,find(cind)).*ones(size(L1));
L1(L1==0) = NaN;
% output images

v0.fname = [prefix 'probability_map.nii'];
v0.dt    = [64 0];
p1       = p(:,1);
p1(p1==0)= NaN;
spm_write_vol(v0,reshape(p1,v0.dim).*L1);

%{
v0.fname = [prefix 'rsq_map.nii'];
spm_write_vol(v0,reshape(r_squared(:,1),v0.dim).*L1);
%}

v0.fname = [prefix 'f_statistic.nii'];
v0.dt    = [64 0];
F_statistic(:,1) = F_statistic(:,1).*p1;
spm_write_vol(v0,reshape(F_statistic(:,1),v0.dim).*L1);



