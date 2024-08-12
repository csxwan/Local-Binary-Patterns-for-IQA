function criteria = getCriteria( referenceScore, mapScore )    
score_ref = referenceScore(:);
score = mapScore(:);

% PLCC: Pearson's Linear Correlation Coefficient
plcc = corr(score_ref,score,'type','Pearson');

% SROCC: Spearman's Rank Correlation Coefficient
srocc = corr(score_ref,score,'type','Spearman');

% KROCC: Kendall rank-order correlation coefficient
krocc=corr(score_ref,score,'type','Kendall');

% RMSE: root mean square error
rmse = sqrt( mean( (score_ref-score).^2 ) );

% MAE: mean absolute error
mae = mean( abs(score_ref-score) );

criteria.plcc = abs(plcc);
criteria.srocc = abs(srocc);
criteria.krocc = abs(krocc);
criteria.rmse = rmse;
criteria.mae = mae;

end