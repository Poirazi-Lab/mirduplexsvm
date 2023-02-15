function miRnaDuplexEst = mirnaduplexfindertest(miRnaDuplexFindTrainConfig, ...
    miRnaDuplexFinderModel, miRnaDuplexTestConfig)
%MIRNADUPLEXFINDERTEST Run miRNA:miRNA*-duplex finder test configurations

% get number of tests
numTests = size(miRnaDuplexTestConfig, 1);

% initialize estimated miRNA:miRNA* duplex dataset
miRnaDuplexEst = dataset( ...
    {cell(numTests, 1), 'miRnaDuplex5pStrand5pEndPosEst'}, ...
    {cell(numTests, 1), 'miRnaDuplex5pStrand3pEndPosEst'}, ...
    {cell(numTests, 1), 'miRnaDuplex3pStrand5pEndPosEst'}, ...
    {cell(numTests, 1), 'miRnaDuplex3pStrand3pEndPosEst'});


parfor i = 1:numTests % for each configuration    
    fprintf('\nPerforming miRNA:miRNA*-duplex finder testing %d of %d', i, numTests);
    
    miRnaDuplexEst(i, :) = mirnaduplexfindertest_single(miRnaDuplexEst(i, :),...
        miRnaDuplexFindTrainConfig(i, :), miRnaDuplexFinderModel(i, :), ...
        miRnaDuplexTestConfig(i, :));
    
end 

 end

