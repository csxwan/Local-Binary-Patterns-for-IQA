function runDatasetSVR_eldp
%% Setup environment
close all;
clear all;

setExpEnv;

%% Feature type and parameters
featureFcnPath = genpath( 'feature/ELDP' );
% featureType = 'LBP_DBC_SCHARR';
% featureFcn = @(img, opts) lbp_dbc_scharr( opts.dim_mfs, opts.radius, opts.n_sampling, img );
featureType = 'ELDP';
featureFcn = @(img, opts)eldp(opts.radius, opts.n_sampling, img );

saveDir = [ 'result', filesep, featureType ];
featuresDir = [ saveDir, filesep, 'features' ];

featureOpts2StringFcn = @(opts) [num2str(opts.radius), ...
    '_', num2str(opts.n_sampling) ];

% featureOptsAll is a cell that collects all feature options
featureOptsAll = {};
for radius = 1
    for n_sampling = 8
%         if n_sampling > (2*radius-1)*(2*radius-1)-1
%             continue;
%         end
        featureOpts.radius=radius;
        featureOpts.n_sampling=n_sampling;
        featureOpts.toString = featureOpts2StringFcn( featureOpts );
        featureOptsAll = [ featureOptsAll, {featureOpts} ];
    end
end




featureInfo.featureType = featureType;
featureInfo.featureFcn = featureFcn;
featureInfo.featureOptsAll = featureOptsAll;
featureInfo.featuresDir = featuresDir;
featureInfo.getFeatureFileNameFcn = @getFeatureFileName;


%% Datasets
% datasetNames = { 'A57', 'WIQ', 'MICT', 'IVC', 'LIVE', 'CSIQ', 'TID2008','TID2013'};
datasetNames = { 'LIVE'};

% Check all the dataset
for idx_dataset = 1:numel( datasetNames )
    datasetName = datasetNames{idx_dataset};
    checkGetDataset( datasetName );
end
tic;
%% Loop for feature extraction on all datasets
addpath( featureFcnPath );
runDataset_extractFeatures( datasetNames, featureInfo );
rmpath( featureFcnPath );

%% Loop for computing SVRscores on all datasets
runDataset_SVR(datasetNames, featureInfo);


toc;
end


