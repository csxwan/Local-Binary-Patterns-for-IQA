function rundataset_hamminglbp
%% Setup environment
close all;
clear all;

setExpEnv;

%% Feature type and parameters
featureFcnPath = genpath( 'feature/HammingLBP' );
% featureType = 'LBP_DBC_SCHARR';
% featureFcn = @(img, opts) lbp_dbc_scharr( opts.dim_mfs, opts.radius, opts.n_sampling, img );
featureType = 'HammingLBP_INT';
featureFcn = @(img, opts) hamminglbp(opts.radius, opts.n_sampling, img );

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
% Recomend [ 8 5 8 ]

% Add one more options
% featureOpts.dim_mfs=8;
% featureOpts.radius=5;
% featureOpts.n_sampling=32;
% featureOpts.toString = featureOpts2StringFcn( featureOpts );
% featureOptsAll = [ featureOptsAll, {featureOpts} ];


featureInfo.featureType = featureType;
featureInfo.featureFcn = featureFcn;
featureInfo.featureOptsAll = featureOptsAll;
featureInfo.featuresDir = featuresDir;
featureInfo.getFeatureFileNameFcn = @getFeatureFileName;

%% Score type and parameters
scoreType = 'norm';
scoreFcn = @norm_util;

scoresDir = [ saveDir, filesep, 'scores_', scoreType ];

% scoreOpts_L0.p = 0;
% scoreOpts_L0.toString = 'L0';
scoreOpts_L1.p = 1;
scoreOpts_L1.toString = 'L1';
% scoreOpts_L2.p = 2;
% scoreOpts_L2.toString = 'L2';
% scoreOpts_Linf.p = Inf;
% scoreOpts_Linf.toString = 'Linf';
% scoreOptsAll = { scoreOpts_L0, scoreOpts_L1, scoreOpts_L2, scoreOpts_Linf };
scoreOptsAll = { scoreOpts_L1 };

scoreInfo.scoreType = scoreType;
scoreInfo.scoreFcn = scoreFcn;
scoreInfo.scoreOptsAll = scoreOptsAll;
scoreInfo.scoresDir = scoresDir;
scoreInfo.getScoreFileNameFcn = @getScoreFileName;

%% Mapped score
mapScoresDir = [ scoreInfo.scoresDir '_map' ];
mapScoreInfo.scoresDir = mapScoresDir;
mapScoreInfo.mapFcn = @(b,x) ( b(1).*(0.5-1./( 1+exp( b(2).*(x-b(3)) ) ))+b(4).*x+b(5) );
mapScoreInfo.initBetaFcn = @() rand(5,1);

%% Datasets
% datasetNames = { 'A57', 'WIQ', 'MICT', 'IVC', 'LIVE', 'CSIQ', 'TID2008','TID2013'};
datasetNames = {'LIVE'};

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

%% Loop for computing scores on all datasets
runDataset_computeScores( datasetNames, featureInfo, scoreInfo );

%% Loop for evaluation on all datasets
runDataset_simpleEvaluation( datasetNames, featureInfo, scoreInfo );

runDataset_computeMapScores( datasetNames, featureInfo, scoreInfo, mapScoreInfo );
scoreInfo2 = scoreInfo;
scoreInfo2.scoresDir = mapScoreInfo.scoresDir;
runDataset_simpleEvaluation( datasetNames, featureInfo, scoreInfo2 );
toc;
end
