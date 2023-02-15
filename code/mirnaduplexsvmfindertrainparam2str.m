function str = mirnaduplexsvmfindertrainparam2str(Param)
%MIRNADUPLEXSVMFINDERTRAINPARAM2STR Convert miRNA:miRNA*-duplex SVM finder parameters to string
%% input
% clear all;
% clc
% load('mirnaduplexfindertrainconfig_crossval', 'miRnaDuplexFinderTrainConfig');
% Param = miRnaDuplexFinderTrainConfig.trainParam{1};
%%
str = '';

if Param.Ratio ~= 100
    str = [str sprintf('ratio-%d', Param.Ratio)];
end

% if Param.FlankingSequenceLength ~= 12
    if ~isempty(str), str = [str '-']; end
    str = [str sprintf('flankseqlen-%d', Param.FlankingSequenceLength)];
% end

if ~strcmp(Param.SvmTrainParam.KernelType, 'radbas')
    if ~isempty(str), str = [str '-']; end
    str = [str sprintf('svmkerneltype-%s', Param.SvmTrainParam.KernelType)];
end

% if strcmp(Param.SvmTrainParam.KernelType, 'poly')
    if ~isempty(str), str = [str '-']; end
    str = [str sprintf('svmdegree-%d', Param.SvmTrainParam.Degree)];
% end

if ~isempty(str), str = [str '-']; end
str = [str sprintf('svmcost-%d', Param.SvmTrainParam.cost)];

% encode more params here

end