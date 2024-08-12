function  runDataset_NRRF( datasetNames, featureInfo )
featureType = featureInfo.featureType;
featureOptsAll = featureInfo.featureOptsAll;
featuresDir = featureInfo.featuresDir;
getFeatureFileNameFcn = featureInfo.getFeatureFileNameFcn;

disp( 'Computing RFscores ...' );
index_of_xlsx = 1;
for idx_dataset = 1:numel( datasetNames )
    datasetName = datasetNames{idx_dataset};
    disp( [ '************************' datasetName '************************' ] );
    dataset = getDataset( datasetName );
    scores_ref = dataset.scores;
    for idx_opts = 1:numel( featureOptsAll )
        featureOpts = featureOptsAll{idx_opts};
        featureFname = getFeatureFileNameFcn( datasetName, featureType, featureOpts );
        saveFeaturesPath = [ featuresDir,filesep,featureFname];
        nTree = 50;
        kfold =5;
        feat = load( saveFeaturesPath );
        features_dist = feat.features_dist;
        [mdata,~]=size(features_dist );
        data=1:mdata;
        %The data samples are randomly divided into k parts.
        indices = crossvalind('Kfold',mdata,kfold);
        for k = 1 : kfold
        %Get the index logical value of the kth test data.
        test = (indices == k);
        % Inverse to obtain the index logical value of the kth training
        % data.
        train = ~test;
        % 1 test, k-1 training
        test_index = data(1,test);
        train_index = data(1,train);
        % Code using data
        num_test=length(test_index);
        num_train=length(train_index);
        temp=1;
        for idex=1:num_test 
                idn=test_index(idex); 
                test_x(temp,:) = features_dist(idn,:); 
                test_label(temp,1) =scores_ref(idn);
                temp=temp+1;
        end

        temp=1;
       for idex=1:num_train
                idn=train_index(idex); 
                train_x(temp,:) = features_dist(idn,:); 
                train_label(temp,1) =scores_ref(idn);
                temp=temp+1;
       end
        %% III. data normalization
        [Train_matrix,PS] = mapminmax(train_x');
        Train_matrix = Train_matrix';
        Test_matrix = mapminmax('apply',test_x',PS);
        Test_matrix = Test_matrix';
        Factor = TreeBagger(nTree, Train_matrix, train_label,'Method','regression');
        [predict_label_2,~] = predict(Factor, Test_matrix);
        criteria(k) = getCriteria( test_label, predict_label_2 );
        end

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

