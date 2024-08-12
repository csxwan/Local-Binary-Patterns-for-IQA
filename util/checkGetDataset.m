function checkGetDataset( datasetName)

dataset = getDataset( datasetName );
img_path_ref = dataset.img_path_ref;
img_path_dist = dataset.img_path_dist;
to_ref = dataset.to_ref;
scores_ref = dataset.scores;
scores_ref_type = dataset.scores_type;

%% test result of getDataset
disp( datasetName );
disp( dataset );

if numel( img_path_ref ) ~= numel( unique(img_path_ref) )
    error( 'path of reference image is not unique!' );
end

if numel( img_path_dist ) ~= numel( unique(img_path_dist) )
    error( 'path of distortion image is not unique!' );
end

for i = 1:numel(img_path_ref)
    img_path = img_path_ref{i};
    if ~exist( img_path, 'file' )
        error( [ img_path, 'does not exists!' ] );
    end
end

for i = 1:numel(img_path_dist)
    img_path = img_path_dist{i};
    if ~exist( img_path, 'file' )
        error( [ img_path, 'does not exists!' ] );
    end
end

for i = 1:numel(to_ref)
    img_path = img_path_ref{to_ref(i)};
    if ~exist( img_path, 'file' )
        error( [ img_path, 'does not exists!' ] );
    end
end

disp( [ 'getDataset(''', datasetName, ''') is OK!' ] );

end
