function [ features_ref, features_dist ] = getFeatures( featureFcn, featureOpts, img_path_ref, img_path_dist )
n_ref = numel(img_path_ref);
feas_ref = cell( n_ref, 1);
if isParallel()
    parfor i = 1:n_ref
        img=imread(img_path_ref{i});
        fea = feval( featureFcn,img,featureOpts );
        feas_ref{i}=reshape( fea, [ 1,numel(fea) ] );
    end % parfor img_path_ref
else
    for i = 1:n_ref
        img=imread(img_path_ref{i});
        fea = feval( featureFcn,img,featureOpts );
        feas_ref{i}=reshape( fea, [ 1,numel(fea) ] );
    end % for img_path_ref
end

n_dist = numel(img_path_dist);
feas_dist = cell( n_dist, 1);

if isParallel()
    parfor i = 1:n_dist
        img=imread(img_path_dist{i});
        fea = feval( featureFcn,img,featureOpts );
        feas_dist{i}=reshape( fea, [ 1,numel(fea) ] );
   end % parfor img_path_dist
else
    for i = 1:n_dist
        img=imread(img_path_dist{i});
        fea = feval( featureFcn,img,featureOpts );
        feas_dist{i}=reshape( fea, [ 1,numel(fea) ] );
    end % for img_path_dist
end
features_ref = cell2mat(feas_ref);
features_dist = cell2mat(feas_dist);
end
