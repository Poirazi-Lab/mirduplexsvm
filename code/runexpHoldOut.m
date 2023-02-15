%% Create parameters
clc, clear all, close all;

crossValK = 1;
holdOutP = 0.3;

% create dataset of miRNA:miRNA* duplex finders
miRnaDuplexFinder = dataset( ...
    {{'SVM'}, 'label'}, ...
    {{@org.mensxmachina.mirna.mirnaduplexsvmfindertrainq}, 'trainFun'}, ...
    {{@mirnaduplexsvmfindertrainparam2str}, 'trainParam2StrFun'}, ...
    {{@org.mensxmachina.mirna.mirnaduplexsvmfindq}, 'findFun'}, ...
    'ObsNames', {'miRnaDuplexSvmFinder'} ...
    );

% create dataset of miRNA:miRNA* duplex strands
miRnaDuplexStrand = dataset( ...
    {{'5'' strand'; '3'' strand'}, 'label'}, ...
    {[5; 3], 'number'}, ...
    'ObsNames', {'strand5p', 'strand3p'} ...
    );

% create dataset of strand ends
strandEnd = dataset( ...
    {{'5'' end'; '3'' end'}, 'label'}, ...
    {[5; 3], 'number'}, ...
    'ObsNames', {'end5p', 'end3p'} ...
    );

save('param', 'miRnaDuplexFinder', 'miRnaDuplexStrand', 'strandEnd', 'crossValK', 'holdOutP');

%% create miRNA:miRNA* finder training configurations
clc, clear all, close all;
flankSeq = 10;
degree = 3;
cost = 100;

load('param');
% create default configuration
Config_default = dataset( ...
    {{'train1'}, 'sampleName'}, ...
    {{'miRnaDuplexSvmFinder'}, 'finderName'}, ...
    {{[]}, 'trainFun'}, ...
    {{struct()}, 'trainParam'}, ...
    {{[]}, 'trainParam2StrFun'} ...
    );

miRnaDuplexFinderTrainConfig = repmat(Config_default, 0, 1);

load(['output' filesep 'data' filesep 'hairpin_train'], 'hairpin');
Config_new = Config_default;
Config_new.sampleName{1} = sprintf('train');
Config_new.finderName{1} = miRnaDuplexFinder.Properties.ObsNames{1};
Config_new.trainFun{1} = miRnaDuplexFinder.trainFun{1};
if strcmp(miRnaDuplexFinder.Properties.ObsNames{1}, 'miRnaDuplexSvmFinder')
    candidateMiRnaDuplexCacheFilename = cellfun(@(name) ['cache' ...
        filesep 'candidatemirnaduplex_hairpin_' name ...
        '_AllTrain_train'], hairpin.Properties.ObsNames, ...
        'UniformOutput', false);
    Config_new.trainParam{1} = struct(...
        'CandidateMiRnaDuplexCaching', true, ...
        'CandidateMiRnaDuplexCacheFilename', {candidateMiRnaDuplexCacheFilename}, ...
        'FlankingSequenceLength', flankSeq, ...
        'Ratio', 100, ...
        'SvmTrainParam', struct('KernelType', 'poly', 'Degree', degree, 'cost', cost), ...
        'Verbose', false ...
        );            
else
    error('Assertion error'); % cannot happen
end

Config_new.trainParam2StrFun{1} = miRnaDuplexFinder.trainParam2StrFun{1};
miRnaDuplexFinderTrainConfig = [miRnaDuplexFinderTrainConfig; Config_new];

% save configurations
save('mirnaduplexfindertrainconfig', 'miRnaDuplexFinderTrainConfig');

%% train finders
clc, clear all, close all;

load('mirnaduplexfindertrainconfig', 'miRnaDuplexFinderTrainConfig');

miRnaDuplexFinderModel = mirnaduplexfindertrain(miRnaDuplexFinderTrainConfig);

save('mirnaduplexfindermodel', 'miRnaDuplexFinderModel');

%% create miRNA:miRNA* finding configurations

clc, clear all, close all;
load('param');

% create default configuration
Config_default = dataset( ...
    {{'test1'}, 'sampleName'}, ...
    {{[]}, 'findFun'}, ...
    {{struct()}, 'findParam'} ...
    );

miRnaDuplexFinderTestConfig = repmat(Config_default, 0, 1);

load(['output' filesep 'data' filesep 'hairpin_test'], 'hairpin');
candidateMiRnaDuplexCacheFilename = cellfun(@(name) ['cache' filesep ...
    'candidatemirnaduplex_hairpin_' name '_HoldOut_test' ...
    ], hairpin.Properties.ObsNames, 'UniformOutput', false);

Config_new = Config_default;
Config_new.sampleName{1} = sprintf('test');
Config_new.findFun{1} = miRnaDuplexFinder.findFun{1};
if strcmp(miRnaDuplexFinder.Properties.ObsNames{1}, 'miRnaDuplexSvmFinder')
    Config_new.findParam{1} = struct(...
        'CandidateMiRnaDuplexCaching', true, ...
        'CandidateMiRnaDuplexCacheFilename', {candidateMiRnaDuplexCacheFilename}, ...
        'Verbose', false);

else
    error('Assertion error'); % cannot happen
end
miRnaDuplexFinderTestConfig = [miRnaDuplexFinderTestConfig; Config_new];

% save configurations
save('mirnaduplexfindertestconfig', 'miRnaDuplexFinderTestConfig');

%% test finders 
clc, clear all, close all;

load('mirnaduplexfindertrainconfig', 'miRnaDuplexFinderTrainConfig');
load('mirnaduplexfindermodel', 'miRnaDuplexFinderModel');
load('mirnaduplexfindertestconfig', 'miRnaDuplexFinderTestConfig');

miRnaDuplexEst = mirnaduplexfindertest(miRnaDuplexFinderTrainConfig, ...
    miRnaDuplexFinderModel, miRnaDuplexFinderTestConfig);

save('mirnaduplexest', 'miRnaDuplexEst');

%% plot errors for all or for every corner separately

clear all, close all;clc;
load('param');
load('mirnaduplexfindertrainconfig', 'miRnaDuplexFinderTrainConfig');
load('mirnaduplexfindertestconfig', 'miRnaDuplexFinderTestConfig');
load('mirnaduplexest', 'miRnaDuplexEst');

% set exp Name 
name = pwd;
dach = strfind(name,'/');
expName = name(dach(end - 1)+1:dach(end)-1)

hold_out = true;

errorPlotter(miRnaDuplexFinderTrainConfig, miRnaDuplexFinderTestConfig, ...
    miRnaDuplexEst, expName, true)

[meanAbsError absErrorStd] = errorPlotterCumFred(miRnaDuplexFinderTrainConfig, ...
    miRnaDuplexFinderTestConfig, miRnaDuplexEst, ...
    expName, true); 



















