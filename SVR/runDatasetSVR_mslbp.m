function runDatasetSVR_mslbp
%% Setup environment
close all;
clear all;

setExpEnv;

%% Feature type and parameters
featureFcnPath = genpath( 'feature/MSLBP' );
% featureType = 'LBP_DBC_SCHARR';
% featureFcn = @(img, opts) lbp_dbc_scharr( opts.dim_mfs, opts.radius, opts.n_sampling, img );
featureType = 'MSLBP';
mapping4=getmapping(4,'riu2');
mapping8=getmapping(8,'riu2');
mapping16=getmapping(16,'riu2');
mapping24=getmapping(24,'riu2');
featureFcn = @(img, opts)mslbp(opts.radius, opts.n_sampling, img ,mapping4 ,mapping8,mapping16,mapping24);

saveDir = [ 'result', filesep, featureType ];
featuresDir = [ saveDir, filesep, 'features' ];

featureOpts2StringFcn = @(opts) [num2str(opts.radius), ...
    '_', num2str(opts.n_sampling) ];

% featureOptsAll is a cell that collects all feature options
featureOptsAll = {};
for radius = 3
    for n_sampling = 24
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


