%% Input:
% indata: Feature matrix with rows = trials, columns = features
% responseData: Vector containing labels
% Classnames: Labels to decode, e.g. [0;1] for binary classification
%% Output:
% trainedClassifier: returns the trained classifier
% results: vector containing classifier performance (accuracy, AUC, ...),
% see below
% predictorWeights: weights assigned by the classifier for each feature

function [trainedClassifier, results, predictorWeights] = ft_subspacediscriminant(indata, responseData, Classnames)

for i = 1:size(indata,2)
    predictorNames{1,i} = ['feature_',num2str(i)];
end
isCategoricalPredictor = false(1,size(predictorNames,2));

inputTable = array2table(indata, 'VariableNames', predictorNames);
predictors = inputTable(:, predictorNames);
response = responseData(:);
% This code specifies all the classifier options and trains the classifier.
subspaceDimension = max(1, min((length(isCategoricalPredictor)/2), width(predictors) - 1));
NumLearningCycles = 30;

t = templateDiscriminant('DiscrimType','linear');

% Set up holdout validation
cvp = cvpartition(response, 'Holdout', 0.1);
trainingPredictors = predictors(cvp.training, :);
trainingResponse = response(cvp.training, :);

classificationEnsemble = fitcensemble(...
    trainingPredictors, ...
    trainingResponse, ...
    'Method', 'Subspace', ...
    'NumLearningCycles', NumLearningCycles, ...
    'Learners', t, ...
    'NPredToSample', subspaceDimension, ...
    'ClassNames', Classnames);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationEnsemble = classificationEnsemble;

% Create the result struct with predict function
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
validationPredictFcn = @(x) ensemblePredictFcn(x);

% Compute validation predictions
validationPredictors = predictors(cvp.test, :);
validationResponse = response(cvp.test, :);
[validationPredictions, validationScores] = validationPredictFcn(validationPredictors);

%% Determine classifier performance
% Compute validation accuracy
correctPredictions = (validationPredictions == validationResponse);
isMissing = isnan(validationResponse);
correctPredictions = correctPredictions(~isMissing);
validationAccuracy = sum(correctPredictions)/length(correctPredictions);

C = confusionmat(validationPredictions,validationResponse);
% figure
% confusionchart(C)
acc = (C(2,2) + C(1,1)) / (C(1,1) + C(1,2) + C(2,1) + C(2,2));
class1precision = C(1,1) / (C(1,1) + C(2,1));
class1recall = C(1,1) / (C(1,1) + C(1,2));
class2precision = C(2,2) / (C(2,2) + C(1,2));
class2recall = C(2,2) / (C(2,2) + C(2,1));

F1 = (2*class1precision*class1recall / (class1precision + class1recall));
F2 = (2*class2precision*class2recall / (class2precision + class2recall));

sensitivity = class1recall;
specificity = class2recall;
% [tpr,fpr,thresholds] = roc(validationPredictions, validationResponse);
% [x,y,t,AUC] = perfcurve(validationPredictions, validationResponse, 1);
% figure
% plot(x,y)

%% Extract feautre weights
all_featureweights = nan(size(classificationEnsemble.Trained,1), length(predictorNames));
for j = 1:size(classificationEnsemble.Trained,1)
    for i = 1:subspaceDimension
        x = classificationEnsemble.Trained{j,1}.PredictorNames{1,i};
        if length(x) == 9
            x = str2num(x(end));
        elseif length(x) == 10
            x = str2num(x(end - 1:end));
        elseif length(x) == 11
            x = str2num(x(end - 2:end));
        elseif length(x) == 12
            x = str2num(x(end - 3:end));
        elseif length(x) == 13
            x = str2num(x(end - 4:end));
        end
        all_featureweights(j,x) = classificationEnsemble.Trained{j,1}.DeltaPredictor(1,i);
    end
end
predictorWeights = nanmean(all_featureweights,1);

%% Decoding results vector
results = [acc, AUC, F1, F2, sensitivity, specificity, class1precision, class1recall, class2precision, class2recall];