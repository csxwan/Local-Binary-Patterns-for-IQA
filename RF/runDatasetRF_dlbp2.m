function runDatasetRF_dlbp2
%% Setup environment
close all;
clear all;

setExpEnv;

%% Feature type and parameters
featureFcnPath = genpath( 'feature/DLBP2' );
featureType = 'DLBP2';
% Image size filtered by circular Gabor filter
size1=[10 10];
featureFcn = @(img, opts)dlbp2(opts.radius, opts.n_sampling, img,size1 );

saveDir = [ 'result', filesep, featureType ];
featuresDir = [ saveDir, filesep, 'features' ];

featureOpts2StringFcn = @(opts) [num2str(opts.radius), ...
    '_', num2str(opts.n_sampling) ];

% featureOptsAll is a cell that collects all feature options
featureOptsAll = {};
for radius = 2
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
runDataset_RFdlbp2(datasetNames, featureInfo,size1);
toc;
end


