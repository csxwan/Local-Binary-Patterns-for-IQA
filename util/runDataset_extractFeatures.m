function runDataset_extractFeatures( datasetNames, featureInfo )
featureType = featureInfo.featureType;
featureFcn = featureInfo.featureFcn;
featureOptsAll = featureInfo.featureOptsAll;
featuresDir = featureInfo.featuresDir;
getFeatureFileNameFcn = featureInfo.getFeatureFileNameFcn;

if ~exist( featuresDir, 'dir' )
    mkdir( featuresDir );
end

disp( 'Extracting features ...' );
for idx_dataset = 1:numel( datasetNames )
    datasetName = datasetNames{idx_dataset};
    disp( [ '************************' datasetName '************************' ] );
    dataset = getDataset( datasetName );
    img_path_ref = dataset.img_path_ref;
    img_path_dist = dataset.img_path_dist;
    for idx_opts = 1:numel( featureOptsAll )
        featureOpts = featureOptsAll{idx_opts};
        featureFname = getFeatureFileNameFcn( datasetName, featureType, featureOpts );
        saveFeaturesPath=[ featuresDir,filesep,featureFname];
        disp( [ 'Dealing with ', featureFname, '...' ] );
        if exist( saveFeaturesPath, 'file' )
            disp( [ 'Skip running feature extraction! ' saveFeaturesPath ' exists!' ]  );
        else
            %getFeatures( featureFcn, featureOpts, img_path_ref, img_path_dist );
            [ features_ref, features_dist ] = getFeatures( featureFcn, featureOpts, img_path_ref, img_path_dist );
            save( saveFeaturesPath, 'features_ref', 'features_dist' );
            clear features_ref features_dist;
            disp( [ saveFeaturesPath,' is OK!' ] );
        end
    end    
end
disp( 'Finished feature extraction.' );
end