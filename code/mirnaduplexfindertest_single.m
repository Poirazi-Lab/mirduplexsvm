function miRnaDuplexEstDs = mirnaduplexfindertest_single(miRnaDuplexEstDs, ...
    miRnaDuplexFindTrainConfig, miRnaDuplexFinderModel, miRnaDuplexTestConfig)

% create cache filename
cacheFilename = ['cache' filesep 'mirnaduplexest_sample_' ...
    lower(miRnaDuplexTestConfig.sampleName{1}) '_trainsample_' ...
    lower(miRnaDuplexFindTrainConfig.sampleName{1}) '_finder_' ...
    lower(miRnaDuplexFindTrainConfig.finderName{1})]; 

miRnaDuplexFinderTrainParamStr = miRnaDuplexFindTrainConfig.trainParam2StrFun{1}...
    (miRnaDuplexFindTrainConfig.trainParam{1});

if ~isempty(miRnaDuplexFinderTrainParamStr)
	cacheFilename = [cacheFilename '_trainparam_' lower(miRnaDuplexFinderTrainParamStr)];   
end

if exist(sprintf('%s.mat', cacheFilename), 'file')
    
    fprintf('\nLoading from cache...\n');

    load(cacheFilename, 'miRnaDuplexEst');
    
else
    
    % load sample
    load(['output' filesep 'data' filesep 'hairpin_' ...
        miRnaDuplexTestConfig.sampleName{1}], 'hairpin');
    load(['output' filesep 'data' filesep 'hairpinbracket_' ...
        miRnaDuplexTestConfig.sampleName{1}], 'hairpinBracket');
    
    % apply finder
    miRnaDuplexEst = miRnaDuplexTestConfig.findFun{1}...
        (miRnaDuplexFinderModel.model{1}, hairpin.sequence, ...
        hairpinBracket.bracket, miRnaDuplexTestConfig.findParam{1});
    
    save(cacheFilename, 'miRnaDuplexEst');
    
end

% put to dataset
miRnaDuplexEstDs.miRnaDuplex5pStrand5pEndPosEst{1} = miRnaDuplexEst(:, 1);
miRnaDuplexEstDs.miRnaDuplex5pStrand3pEndPosEst{1} = miRnaDuplexEst(:, 2);
miRnaDuplexEstDs.miRnaDuplex3pStrand5pEndPosEst{1} = miRnaDuplexEst(:, 3);
miRnaDuplexEstDs.miRnaDuplex3pStrand3pEndPosEst{1} = miRnaDuplexEst(:, 4);

end







