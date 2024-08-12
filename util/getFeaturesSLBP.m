function [ features_ref, features_dist ] = getFeaturesSLBP( featureFcn, featureOpts, img_path_ref, img_path_dist )
%function getFeatures( featureFcn, featureOpts, img_path_ref, img_path_dist )
n_ref = numel(img_path_ref);
feas_ref = cell( n_ref, 1);
smap_path_ref = cell( n_ref, 1);
for i=1:n_ref
A=img_path_ref{i};
A1=A(1:7);
A2=A(8:end);
a=[A1,'src',A2];
smap_path_ref{i}=a;
end
% feas_ref_lbp = cell(n_ref,1);
% feas_ref_cnt = cell(n_ref,1);
if isParallel()
    parfor i = 1:n_ref
        img=imread(img_path_ref{i});
        simg=imread(smap_path_ref{i});
        fea = feval( featureFcn,img,simg,featureOpts );
        feas_ref{i}=reshape( fea, [ 1,numel(fea) ] );

        % [fea_lbp,fea_cnt]=feval( featureFcn,img,featureOpts );
        % feas_ref_lbp{i}=reshape( fea_lbp, [ 1,numel(fea_lbp) ] );
        % feas_ref_cnt{i}=reshape( fea_cnt, [ 1,numel(fea_cnt) ] );
    end % parfor img_path_ref
else
    for i = 1:n_ref
        img=imread(img_path_ref{i});
        simg=imread(smap_path_ref{i});
        fea = feval( featureFcn,img,simg,featureOpts );
        feas_ref{i}=reshape( fea, [ 1,numel(fea) ] );

        % [fea_lbp,fea_cnt]=feval( featureFcn,img,featureOpts );
        % feas_ref_lbp{i}=reshape( fea_lbp, [ 1,numel(fea_lbp) ] );
        % feas_ref_cnt{i}=reshape( fea_cnt, [ 1,numel(fea_cnt) ] );
    end % for img_path_ref
end

n_dist = numel(img_path_dist);
feas_dist = cell( n_dist, 1);
smap_path_dist = cell( n_dist, 1);
for i=1:n_dist
A=img_path_dist{i};
A1=A(1:7);
A2=A(8:end);
a=[A1,'src',A2];
smap_path_dist{i}=a;
end
% feas_dist_lbp = cell( n_dist,1);
% feas_dist_cnt = cell( n_dist,1);

if isParallel()
    parfor i = 1:n_dist
        img=imread(img_path_dist{i});
        simg=imread(smap_path_dist{i});
        fea = feval( featureFcn,img,simg,featureOpts );
        feas_dist{i}=reshape( fea, [ 1,numel(fea) ] );

        % [fea_lbp,fea_cnt]=feval( featureFcn,img,featureOpts );
        % feas_dist_lbp{i}=reshape( fea_lbp, [ 1,numel(fea_lbp) ] );
        % feas_dist_cnt{i}=reshape( fea_cnt, [ 1,numel(fea_cnt) ] );
    end % parfor img_path_dist
else
    for i = 1:n_dist
        img=imread(img_path_dist{i});
        simg=imread(smap_path_dist{i});
        fea = feval( featureFcn,img,simg,featureOpts );
        feas_dist{i}=reshape( fea, [ 1,numel(fea) ] );

        % [fea_lbp,fea_cnt]=feval( featureFcn,img,featureOpts );
        % feas_dist_lbp{i}=reshape( fea_lbp, [ 1,numel(fea_lbp) ] );
        % feas_dist_cnt{i}=reshape( fea_cnt, [ 1,numel(fea_cnt) ] );
    end % for img_path_dist
end

features_ref = cell2mat(feas_ref);
features_dist = cell2mat(feas_dist);

% features_ref = {cell2mat(feas_ref_lbp) cell2mat(feas_ref_cnt)};
% features_dist = {cell2mat(feas_dist_lbp) cell2mat(feas_dist_cnt)};

end
