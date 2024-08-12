function runDatasetRF_ocpp
%% Setup environment
close all;
clear all;

setExpEnv;

%% Feature type and parameters
featureFcnPath = genpath( 'feature/OCPP' );
featureType = 'OCPP';
featureFcn = @(img, opts)ocpp(opts.radius,opts.neighborsxy,opts.neighborsxz,opts.neighborsyz, img );

saveDir = [ 'result', filesep, featureType ];
featuresDir = [ saveDir, filesep, 'features' ];

featureOpts2StringFcn = @(opts) [num2str(opts.radius), ...
    '_', num2str(opts.neighborsxy),  '_', num2str(opts.neighborsxz),  '_', num2str(opts.neighborsyz) ];

% featureOptsAll is a cell that collects all feature options
featureOptsAll = {};
for radius = 3
    for neighborsxy = 16
        for neighborsxz = 8
            for neighborsyz = 10
%         if n_sampling > (2*radius-1)*(2*radius-1)-1
%             continue;
%         end
        featureOpts.radius=radius;
        featureOpts.neighborsxy =neighborsxy ;
        featureOpts.neighborsxz =neighborsxz ;
        featureOpts.neighborsyz =neighborsyz ;
        featureOpts.toString = featureOpts2StringFcn( featureOpts );
        featureOptsAll = [ featureOptsAll, {featureOpts} ];
            end
        end
    end
end




featureInfo.featureType = featureType;
featureInfo.featureFcn = featureFcn;
featureInfo.featureOptsAll = featureOptsAll;
featureInfo.featuresDir = featuresDir;
featureInfo.getFeatureFileNameFcn = @getFeatureFileName;


%% Datasets
% datasetNames = { 'LIVE', 'CSIQ','TID2013'};
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


