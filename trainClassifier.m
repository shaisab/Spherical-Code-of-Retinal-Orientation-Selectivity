function [trainedClassifier, validationAccuracy] = trainClassifier(datasetTable)
% Extract predictors and response
predictorNames = {'Var1', 'Var2', 'Var3', 'Var4', 'Var5', 'Var6', 'Var7', 'Var8', 'Var9', 'Var10', 'Var11', 'Var12', 'Var13', 'Var14', 'Var15'};
predictors = datasetTable(:,predictorNames);
predictors = table2array(varfun(@double, predictors));
response = datasetTable.Var16;
% Train a classifier
trainedClassifier = fitensemble(predictors, response, 'Subspace', 200, 'KNN', 'Type', 'Classification', 'NPredToSample', 8, 'PredictorNames', {'Var1' 'Var2' 'Var3' 'Var4' 'Var5' 'Var6' 'Var7' 'Var8' 'Var9' 'Var10' 'Var11' 'Var12' 'Var13' 'Var14' 'Var15'}, 'ResponseName', 'Var16', 'ClassNames', [1 2]);

% Perform cross-validation
partitionedModel = crossval(trainedClassifier, 'KFold', 5);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

save('trainedClassifier.mat','trainedClassifier','partitionedModel','validationAccuracy')

%% Uncomment this section to compute validation predictions and scores:
% % Compute validation predictions and scores
% [validationPredictions, validationScores] = kfoldPredict(partitionedModel);