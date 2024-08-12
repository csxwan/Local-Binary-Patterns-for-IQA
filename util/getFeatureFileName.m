function fname = getFeatureFileName( datasetName, featureType, featureOpts )
fname = [ datasetName,'_',featureType, '_', featureOpts.toString,'.mat'];
end