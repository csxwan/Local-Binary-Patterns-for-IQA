function  runcrossDataset_NRRF( datasetNames,datasetNames1, featureInfo )
featureType = featureInfo.featureType;
featureOptsAll = featureInfo.featureOptsAll;
featuresDir = featureInfo.featuresDir;
getFeatureFileNameFcn = featureInfo.getFeatureFileNameFcn;

disp( 'Computing crossdataset NRRFscores ...' );
index_of_xlsx = 1;
for idx_dataset = 1:numel( datasetNames )
    datasetName = datasetNames{idx_dataset};
    datasetName1 = datasetNames1{idx_dataset};
    disp( [ '************************' datasetName '************************' ] );
    dataset = getDataset( datasetName );
    scores_ref = dataset.scores;
    dataset1 = getDataset( datasetName1 );
    scores_ref1 = dataset1.scores;
    for idx_opts = 1:numel( featureOptsAll )
        featureOpts = featureOptsAll{idx_opts};
        featureFname = getFeatureFileNameFcn( datasetName, featureType, featureOpts );
        featureFname1 = getFeatureFileNameFcn( datasetName1, featureType, featureOpts );
        saveFeaturesPath = [ featuresDir,filesep,featureFname];
        saveFeaturesPath1 = [ featuresDir,filesep,featureFname1];
        nTree = 50;
        kfold =5;
        feat = load( saveFeaturesPath );
        feat1 = load( saveFeaturesPath1 );
        features_dist = feat.features_dist;
        features_dist1 = feat1.features_dist;

        train_x=features_dist;
        test_x=features_dist1;
        train_label=scores_ref;
        test_label=scores_ref1;

        %% III. data normalization
        [Train_matrix,PS] = mapminmax(train_x');
        Train_matrix = Train_matrix';
        Test_matrix = mapminmax('apply',test_x',PS);
        Test_matrix = Test_matrix';
        Factor = TreeBagger(nTree, Train_matrix, train_label,'Method','regression');
        [predict_label_2,~] = predict(Factor, Test_matrix);
        criteria = getCriteria( test_label, predict_label_2 );
        

        mcriteria.plcc= mean([criteria.plcc]);
        mcriteria.srocc = mean([criteria.srocc]);
        mcriteria.krocc = mean([criteria.krocc]);
        mcriteria.rmse = mean([criteria.rmse]);
        mcriteria.mae =mean([ criteria.mae]);
        disp( mcriteria );
        disp( [ saveFeaturesPath,' is OK!' ] );

                %write the criteria into the excel file
                cri_info{1} = featureType;
                cri_info{3} = featureOpts.toString;
                cri_info{4} = datasetName;
                cri_info{5} = mcriteria.plcc;
                cri_info{6} = mcriteria.srocc;
                cri_info{7} = mcriteria.krocc;
                cri_info{8} = mcriteria.rmse;
                cri_info{9} = mcriteria.mae;
                               
                xlswrite('E:\Users\a\temp\SRP\result.xlsx',cri_info, 'Sheet3', num2str(index_of_xlsx));
                index_of_xlsx = index_of_xlsx + 1;

    end
end
disp( 'Finished RFscores computation.' );
end

