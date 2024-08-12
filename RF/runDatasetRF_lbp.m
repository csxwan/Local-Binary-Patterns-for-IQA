function runDatasetRF_lbp
%% Setup environment
close all;
clear all;

setExpEnv;

%% Feature type and parameters
featureFcnPath = genpath( 'feature/LBP2' );
featureType = 'LBP';
featureFcn = @(img, opts)lbp(opts.radius, opts.n_sampling, img );

saveDir = [ 'result', filesep, featureType ];
featuresDir = [ saveDir, filesep, 'features' ];

featureOpts2StringFcn = @(opts) [num2str(opts.radius), ...
    '_', num2str(opts.n_sampling) ];

% featureOptsAll is a cell that collects all feature options
featureOptsAll = {};
for radius = 1:3
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
% datasetNames = {  'LIVE', 'CSIQ', 'TID2013'};
datasetNames = { 'awgn','blur','contrast','fnoise','jpeg','jpeg2000'};
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

%% Loop for computing RFscores on all datasets
runDataset_RF(datasetNames, featureInfo);
toc;
end