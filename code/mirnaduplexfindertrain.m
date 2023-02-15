function miRnaDuplexFinderModel = mirnaduplexfindertrain ... 
   (miRnaDuplexFinderTrainConfig)
%MIRNADUPLEXFINDERTRAIN Run miRNA:miRNA*-duplex finder training configurations
import org.mensxmachina.mirna.*;

% get number of trainings
numTrainings = size(miRnaDuplexFinderTrainConfig, 1);

% initialize model dataset
miRnaDuplexFinderModel = dataset( ...
    {cell(numTrainings, 1), 'model'} ...
);

parfor i = 1:numTrainings % for each training

    fprintf('\nPerforming miRNA:miRNA*-duplex finder training %d of %d', ...
        i, numTrainings);
    
    miRnaDuplexFinderModel(i, :) = ...
        mirnaduplexfindertrain_single(miRnaDuplexFinderModel(i, :), ...
        miRnaDuplexFinderTrainConfig(i, :));

end

end



function miRnaDuplexFinderModel = mirnaduplexfindertrain_single...
   (miRnaDuplexFinderModel, miRnaDuplexFinderTrainConfig)

% create cache filename

cacheFilename = ['cache' filesep 'mirnaduplexfindermodel_trainsample_' ...
    lower(miRnaDuplexFinderTrainConfig.sampleName{1}) ...
    '_finder_' lower(miRnaDuplexFinderTrainConfig.finderName{1})]; 

miRnaDuplexFinderTrainParamStr = ...
    miRnaDuplexFinderTrainConfig.trainParam2StrFun{1}...
    (miRnaDuplexFinderTrainConfig.trainParam{1});

if ~isempty(miRnaDuplexFinderTrainParamStr)
	cacheFilename = [cacheFilename '_trainparam_' ...
        lower(miRnaDuplexFinderTrainParamStr)];   
end

if exist(sprintf('%s.mat', cacheFilename), 'file')
    
    fprintf('\nLoading from cache...\n');
    load(cacheFilename, 'model');
    
else
    
    % load sample
    load(['output' filesep 'data' filesep 'hairpin_' ...
        miRnaDuplexFinderTrainConfig.sampleName{1}], 'hairpin');
    load(['output' filesep 'data' filesep 'hairpinbracket_' ...
        miRnaDuplexFinderTrainConfig.sampleName{1}], 'hairpinBracket');
    
    % create duplex matrix
    miRnaDuplex = [hairpin.miRnaDuplex5pStrand5pEndPos ...
        hairpin.miRnaDuplex5pStrand3pEndPos ...
        hairpin.miRnaDuplex3pStrand5pEndPos ...
        hairpin.miRnaDuplex3pStrand3pEndPos];
    
    % train finder
    model = miRnaDuplexFinderTrainConfig.trainFun{1} ...
        (hairpin.sequence, hairpinBracket.bracket, miRnaDuplex, ...
        miRnaDuplexFinderTrainConfig);
    
    % save model
    save(cacheFilename, 'model');
    
end

% put model to dataset
miRnaDuplexFinderModel.model{1} = model;

end