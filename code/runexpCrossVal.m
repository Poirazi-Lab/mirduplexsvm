%% create miRNA:miRNA* finder training configurations
clc, clear all, close all;

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

flankSeq = 10;
degrees = [1,2,3];
cost = [0.01, 0.1, 1, 10, 100, 1000];
for i = 1:crossValK % for each cross validation fold
    for d = 1:length(degrees)
        for c = 1:length(cost)
            degree = degrees(d);
            load(['output' filesep 'data' filesep 'hairpin_train-' int2str(i)], 'hairpin');
            Config_new = Config_default;
            Config_new.sampleName{1} = sprintf('train-%d', i);
            Config_new.finderName{1} = miRnaDuplexFinder.Properties.ObsNames{1};
            Config_new.trainFun{1} = miRnaDuplexFinder.trainFun{1};
            if strcmp(miRnaDuplexFinder.Properties.ObsNames{1}, 'miRnaDuplexSvmFinder')
                candidateMiRnaDuplexCacheFilename = cellfun(@(name) ['cache' ...
                    filesep 'candidatemirnaduplex_hairpin_' name ...
                    '_trainsample_train-' int2str(i)], hairpin.Properties.ObsNames, ...
                    'UniformOutput', false);
                Config_new.trainParam{1} = struct(...
                    'CandidateMiRnaDuplexCaching', true, ...
                    'CandidateMiRnaDuplexCacheFilename', {candidateMiRnaDuplexCacheFilename}, ...
                    'FlankingSequenceLength', flankSeq, ...
                    'Ratio', 100, ...
                    'SvmTrainParam', struct('KernelType', 'poly', 'Degree', degree, 'cost', cost(c)), ...
                    'Verbose', false);                   
            else
                error('Assertion error'); % cannot happen
            end
            Config_new.trainParam2StrFun{1} = miRnaDuplexFinder.trainParam2StrFun{1};
            miRnaDuplexFinderTrainConfig = [miRnaDuplexFinderTrainConfig; Config_new];
        end
    end
end

% save configurations
save('mirnaduplexfindertrainconfig_crossval', 'miRnaDuplexFinderTrainConfig');

%% train finders
clc, clear all, close all;

load('mirnaduplexfindertrainconfig_crossval', 'miRnaDuplexFinderTrainConfig');

miRnaDuplexFinderModel = mirnaduplexfindertrain(miRnaDuplexFinderTrainConfig);

save('mirnaduplexfindermodel_crossval', 'miRnaDuplexFinderModel');

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
degrees = [2,3];
cost = [0.001, 0.01];
for i = 1:crossValK % for each cross validation fold
    for d = 1:length(degrees)
        for c = 1:length(cost)
            load(['output' filesep 'data' filesep 'hairpin_test-' int2str(i)], 'hairpin');
            candidateMiRnaDuplexCacheFilename = cellfun(@(name) ['cache' filesep ...
                'candidatemirnaduplex_hairpin_' name '_trainsample_train-' ...
                int2str(i)], hairpin.Properties.ObsNames, 'UniformOutput', false);
            Config_new = Config_default;
            Config_new.sampleName{1} = sprintf('test-%d', i);
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
        end
    end
end

% save configurations
save('mirnaduplexfindertestconfig_crossval', 'miRnaDuplexFinderTestConfig');

%% test finders 
clc, clear all, close all;

load('mirnaduplexfindertrainconfig_crossval', 'miRnaDuplexFinderTrainConfig');
load('mirnaduplexfindermodel_crossval', 'miRnaDuplexFinderModel');
load('mirnaduplexfindertestconfig_crossval', 'miRnaDuplexFinderTestConfig');

miRnaDuplexEst = mirnaduplexfindertest(miRnaDuplexFinderTrainConfig, ...
    miRnaDuplexFinderModel, miRnaDuplexFinderTestConfig);

save('mirnaduplexest_crossval', 'miRnaDuplexEst');

%% plot errors for all or for every corner separately

clear all, close all;clc;
load('param');
load('mirnaduplexfindertrainconfig_crossval', 'miRnaDuplexFinderTrainConfig');
load('mirnaduplexfindertestconfig_crossval', 'miRnaDuplexFinderTestConfig');
load('mirnaduplexest_crossval', 'miRnaDuplexEst');

% Create plots for each parameter combination 
num_dif_comb = sum(strcmp('train-1', miRnaDuplexFinderTrainConfig.sampleName));

for i = 1:num_dif_comb
    inds = i:num_dif_comb:(crossValK*num_dif_comb);
    miRnaDuplexFinderTrainConfig_temp = miRnaDuplexFinderTrainConfig(inds,:);
    miRnaDuplexFinderTestConfig_temp = miRnaDuplexFinderTestConfig(inds,:);
    miRnaDuplexEst_temp = miRnaDuplexEst(inds,:);
    % set exp Name 
    name = pwd;
    dach = strfind(name,'/');
    %expName = name(dach(end - 3)+1:dach(end - 2)-1)
    d = miRnaDuplexFinderTrainConfig_temp.trainParam{1}.SvmTrainParam.Degree;
    c = miRnaDuplexFinderTrainConfig_temp.trainParam{1}.SvmTrainParam.cost;
    expName = sprintf('mirnasSEQdegree%dcost%d', d, c);
    filter = expName == '.';
    expName = expName(~filter)
    crossVall = 5;
   
    errorPlotter(miRnaDuplexFinderTrainConfig_temp, miRnaDuplexFinderTestConfig_temp, ...
        miRnaDuplexEst_temp, expName, false)
    [meanAbsError absErrorStd] = errorPlotterCumFred(miRnaDuplexFinderTrainConfig_temp, ...
        miRnaDuplexFinderTestConfig_temp, miRnaDuplexEst_temp, ...
        expName, false); 
    
end

