%% Add necessary paths
clc, clear all;
addpath('mirna')
addpath('libsvm-3.11/matlab')

%% Create parameters
clc, clear all, close all;

crossValK = 5;
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

% create dataset of media
medium = dataset( ...
    {{'Paper'; 'Presentation'}, 'label'}, ...
    {{{'Colormap', [1 1 1; 0.8 0.8 0.8]}; {'PaperSize', [10 8], 'PaperPosition', [0 0 10 8]}}, 'figureProperties'}, ...
    'ObsNames', {'paper', 'presentation'} ...
    );

save('param', 'miRnaDuplexFinder', 'miRnaDuplexStrand', 'strandEnd', 'medium', 'crossValK', 'holdOutP');

%% Read miRBase EMBL file

clc, clear all, close all;

% read miRBase EMBL file
Embl = emblread(['input' filesep 'data' filesep 'miRNA.dat']);

save(['output' filesep 'data' filesep 'embl'], 'Embl');

%% Read miRBase EMBL struct into dataset

clc, clear all, close all;

load(['output' filesep 'data' filesep 'embl'], 'Embl');

hairpin = hairpinembl2dataset(Embl);

save(['output' filesep 'data' filesep 'hairpin'], 'hairpin');

%% Filter hairpins 

clc, clear all, close all;

load(['output' filesep 'data' filesep 'hairpin'], 'hairpin');

% select hairpins with known duplex and not containing bases other than ACGU
filter = hairpin.miRnaDuplexKnown & hairpin.otherBaseCount == 0;
fprintf('\nSelected %d out of %d hairpins.\n', sum(filter), size(hairpin, 1));
hairpin = hairpin(filter, :);

save(['output' filesep 'data' filesep 'hairpin_filtered'], 'hairpin');

%% Select only human and mouse hairpins 

clc, clear all, close all;

load(['output' filesep 'data' filesep 'hairpin_filtered'], 'hairpin');

filter = ismember(hairpin.speciesName, {'hsa', 'mmu'});
fprintf('\nSelected %d out of %d hairpins.\n', sum(filter), size(hairpin, 1));
hairpin = hairpin(filter, :);

save(['output' filesep 'data' filesep 'hairpin_hsa-mmu'], 'hairpin');

%% Fold human and mouse hairpins using RNAfold

clc, clear all, close all;

load(['output' filesep 'data' filesep 'hairpin_hsa-mmu'], 'hairpin');

[hairpinBracket] = hairpinfold(hairpin);

save(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu'], 'hairpinBracket');

%% Calculate human and mouse hairpin 2D features
clc, clear all, close all;

load(['output' filesep 'data' filesep 'hairpin_hsa-mmu'], 'hairpin');
load(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu'], 'hairpinBracket');

hairpin2dFeatures = hairpin2dfeatures(hairpin, hairpinBracket);

save(['output' filesep 'data' filesep 'hairpin2dfeatures_hsa-mmu'], 'hairpin2dFeatures');

%% Filter human and mouse hairpins according to 2D features

clc, clear all, close all;

load(['output' filesep 'data' filesep 'hairpin_hsa-mmu'], 'hairpin');
load(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu'], 'hairpinBracket');
load(['output' filesep 'data' filesep 'hairpin2dfeatures_hsa-mmu'], 'hairpin2dFeatures');

% create filter
filter_folds = hairpin2dFeatures.folds;
fprintf('\n%d out of %d hairpins are unfoldable.\n', size(hairpin, 1) - ...
    sum(filter_folds), size(hairpin, 1));

filter_multibranched = ~hairpin2dFeatures.multibranched;

fprintf('\n%d out of %d hairpins are multibranch.\n', size(hairpin, 1) - ...
    sum(filter_multibranched), size(hairpin, 1));

filter = filter_folds & filter_multibranched;

% filter hairpins
hairpin = hairpin(filter, :);
hairpinBracket = hairpinBracket(filter, :);
hairpin2dFeatures = hairpin2dFeatures(filter, :);

save(['output' filesep 'data' filesep 'hairpin_hsa-mmu-filtered'], 'hairpin');
save(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu-filtered'], 'hairpinBracket');
save(['output' filesep 'data' filesep 'hairpin2dfeatures_hsa-mmu-filtered'], 'hairpin2dFeatures');

%% Split human and mouse hairpins into training set and test set

clc, clear all, close all;

load('param');

load(['output' filesep 'data' filesep 'hairpin_hsa-mmu-filtered'], 'hairpin');
load(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu-filtered'], 'hairpinBracket');
load(['output' filesep 'data' filesep 'hairpin2dfeatures_hsa-mmu-filtered'], 'hairpin2dFeatures');

filter = strcmp(hairpin.speciesName, 'hsa'); % find human hairpins

% select human hairpins
hairpin_hsa = hairpin(filter, :);
hairpinBracket_hsa = hairpinBracket(filter, :);
hairpin2dFeatures_hsa = hairpin2dFeatures(filter, :);

% select mouse hairpins
hairpin_mmu = hairpin(~filter, :);
hairpinBracket_mmu = hairpinBracket(~filter, :);
hairpin2dFeatures_mmu = hairpin2dFeatures(~filter, :);

% reset default random stream
stream = RandStream.getDefaultStream;
stream.reset;

% find training and test set indices for each species seperatelly
[hairpin_hsa_train_ind, hairpin_hsa_test_ind] = ...
    crossvalind('HoldOut', size(hairpin_hsa, 1), holdOutP);

[hairpin_mmu_train_ind, hairpin_mmu_test_ind] = ...
    crossvalind('HoldOut', size(hairpin_mmu, 1), holdOutP);

% create training set
hairpin = [hairpin_hsa(hairpin_hsa_train_ind, :); hairpin_mmu(hairpin_mmu_train_ind, :)];
hairpinBracket = [hairpinBracket_hsa(hairpin_hsa_train_ind, :); ...
    hairpinBracket_mmu(hairpin_mmu_train_ind, :)];
hairpin2dFeatures = [hairpin2dFeatures_hsa(hairpin_hsa_train_ind, :); ...
    hairpin2dFeatures_mmu(hairpin_mmu_train_ind, :)];

% save training set
save(['output' filesep 'data' filesep 'hairpin_train'], 'hairpin');
save(['output' filesep 'data' filesep 'hairpinbracket_train'], 'hairpinBracket');
save(['output' filesep 'data' filesep 'hairpin2dfeatures_train'], 'hairpin2dFeatures');

% create test set
hairpin = [hairpin_hsa(hairpin_hsa_test_ind, :); hairpin_mmu(hairpin_mmu_test_ind, :)];
hairpinBracket = [hairpinBracket_hsa(hairpin_hsa_test_ind, :); ...
    hairpinBracket_mmu(hairpin_mmu_test_ind, :)];
hairpin2dFeatures = [hairpin2dFeatures_hsa(hairpin_hsa_test_ind, :); ...
    hairpin2dFeatures_mmu(hairpin_mmu_test_ind, :)];

% save test set
save(['output' filesep 'data' filesep 'hairpin_test'], 'hairpin');
save(['output' filesep 'data' filesep 'hairpinbracket_test'], 'hairpinBracket');
save(['output' filesep 'data' filesep 'hairpin2dfeatures_test'], 'hairpin2dFeatures');

%% Remove outliers

clc, clear all, close all;

load(['output' filesep 'data' filesep 'hairpin_train'], 'hairpin');
load(['output' filesep 'data' filesep 'hairpinbracket_train'], 'hairpinBracket');
load(['output' filesep 'data' filesep 'hairpin2dfeatures_train'], 'hairpin2dFeatures');

% create filter

% find 5'' strand overhang outliers 
mu = mean(hairpin2dFeatures.miRnaDuplex5pStrandOverhangSequenceLength);
sigma = std(hairpin2dFeatures.miRnaDuplex5pStrandOverhangSequenceLength);
filter_duplex5pStrandOverhang = abs(hairpin2dFeatures.miRnaDuplex5pStrandOverhangSequenceLength - mu) <= 3*sigma;
fprintf('\n%d outliers in 5'' strand overhang (mu=%.2f, sigma=%.2f).\n', size(hairpin, 1) - sum(filter_duplex5pStrandOverhang), mu, sigma);

% find 3'' strand overhang outliers 
mu = mean(hairpin2dFeatures.miRnaDuplex3pStrandOverhangSequenceLength);
sigma = std(hairpin2dFeatures.miRnaDuplex3pStrandOverhangSequenceLength);
filter_duplex3pStrandOverhang = abs(hairpin2dFeatures.miRnaDuplex3pStrandOverhangSequenceLength - mu) <= 3*sigma;
fprintf('\n%d outliers in 3'' strand overhang (mu=%.2f, sigma=%.2f).\n', size(hairpin, 1) - sum(filter_duplex3pStrandOverhang), mu, sigma);


% find 5'strand 5p distance from tip outliers
mu = mean(hairpin2dFeatures.miRnaDuplex5pStrand5pEndSeqDistFromTip);
sigma = std(hairpin2dFeatures.miRnaDuplex5pStrand5pEndSeqDistFromTip);
filter_duplex5stard5pEnd = ...
    abs(hairpin2dFeatures.miRnaDuplex5pStrand5pEndSeqDistFromTip - mu) <= 3*sigma;
fprintf('\n%d outliers in 5'' strand 5p distance from Tip (mu=%.2f, sigma=%.2f).\n', ...
    size(hairpin, 1) - sum(filter_duplex5stard5pEnd), mu, sigma);

% find 5'strand 3p distance from tip outliers
mu = mean(hairpin2dFeatures.miRnaDuplex5pStrand3pEndSeqDistFromTip);
sigma = std(hairpin2dFeatures.miRnaDuplex5pStrand3pEndSeqDistFromTip);
filter_duplex5stard3pEnd = ...
    abs(hairpin2dFeatures.miRnaDuplex5pStrand3pEndSeqDistFromTip - mu) <= 3*sigma;
fprintf('\n%d outliers in 5'' strand 3p distance from Tip (mu=%.2f, sigma=%.2f).\n', ...
    size(hairpin, 1) - sum(filter_duplex5stard3pEnd), mu, sigma);

% find 3'strand 5p distance from tip outliers
mu = mean(hairpin2dFeatures.miRnaDuplex3pStrand5pEndSeqDistFromTip);
sigma = std(hairpin2dFeatures.miRnaDuplex3pStrand5pEndSeqDistFromTip);
filter_duplex3stard5pEnd = ...
    abs(hairpin2dFeatures.miRnaDuplex3pStrand5pEndSeqDistFromTip - mu) <= 3*sigma;
fprintf('\n%d outliers in 3'' strand 5p distance from Tip (mu=%.2f, sigma=%.2f).\n', ...
    size(hairpin, 1) - sum(filter_duplex3stard5pEnd), mu, sigma);

% find 3'strand 3p distance from tip outliers
mu = mean(hairpin2dFeatures.miRnaDuplex3pStrand3pEndSeqDistFromTip);
sigma = std(hairpin2dFeatures.miRnaDuplex3pStrand3pEndSeqDistFromTip);
filter_duplex3stard3pEnd = ...
    abs(hairpin2dFeatures.miRnaDuplex3pStrand3pEndSeqDistFromTip - mu) <= 3*sigma;
fprintf('\n%d outliers in 3'' strand 3p distance from Tip (mu=%.2f, sigma=%.2f).\n', ...
    size(hairpin, 1) - sum(filter_duplex3stard3pEnd), mu, sigma);

filter = filter_duplex5pStrandOverhang & filter_duplex3pStrandOverhang ...
    & filter_duplex5stard5pEnd ...
    & filter_duplex5stard3pEnd ...
    & filter_duplex3stard5pEnd ...
    & filter_duplex3stard3pEnd;

% filter hairpins
hairpin = hairpin(filter, :);
hairpinBracket = hairpinBracket(filter, :);
hairpin2dFeatures = hairpin2dFeatures(filter, :);

% save
save(['output' filesep 'data' filesep 'hairpin_hsa-mmu-filtered-nooutliers'], 'hairpin');
save(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu-filtered-nooutliers'], 'hairpinBracket');
save(['output' filesep 'data' filesep 'hairpin2dfeatures_hsa-mmu-filtered-nooutliers'], 'hairpin2dFeatures');

%% Split training set for cross validation

clc, clear all, close all;

load('param');

load(['output' filesep 'data' filesep 'hairpin_hsa-mmu-filtered-nooutliers'], 'hairpin');
load(['output' filesep 'data' filesep 'hairpinbracket_hsa-mmu-filtered-nooutliers'], 'hairpinBracket');
load(['output' filesep 'data' filesep 'hairpin2dfeatures_hsa-mmu-filtered-nooutliers'], 'hairpin2dFeatures');

filter = strcmp(hairpin.speciesName, 'hsa'); % find human hairpins

% select human hairpins
hairpin_hsa = hairpin(filter, :);
hairpinBracket_hsa = hairpinBracket(filter, :);
hairpin2dFeatures_hsa = hairpin2dFeatures(filter, :);

% select mouse hairpins
hairpin_mmu = hairpin(~filter, :);
hairpinBracket_mmu = hairpinBracket(~filter, :);
hairpin2dFeatures_mmu = hairpin2dFeatures(~filter, :);

% reset default random stream
stream = RandStream.getDefaultStream;
stream.reset;

hairpin_hsa_test_ind = crossvalind('Kfold', size(hairpin_hsa, 1), crossValK);
hairpin_mmu_test_ind = crossvalind('Kfold', size(hairpin_mmu, 1), crossValK);

for i = 1:crossValK % for each cross validation fold
    
    i_hairpin_hsa_test_ind = (hairpin_hsa_test_ind == i);
    i_hairpin_hsa_train_ind = ~i_hairpin_hsa_test_ind;
    
    i_hairpin_mmu_test_ind = (hairpin_mmu_test_ind == i);
    i_hairpin_mmu_train_ind = ~i_hairpin_mmu_test_ind;
    
    % create i-th training set
    hairpin = [hairpin_hsa(i_hairpin_hsa_train_ind, :); hairpin_mmu(i_hairpin_mmu_train_ind, :)];
    hairpinBracket = [hairpinBracket_hsa(i_hairpin_hsa_train_ind, :); hairpinBracket_mmu(i_hairpin_mmu_train_ind, :)];
    hairpin2dFeatures = [hairpin2dFeatures_hsa(i_hairpin_hsa_train_ind, :); hairpin2dFeatures_mmu(i_hairpin_mmu_train_ind, :)];
    
    % save i-th training set
    save(['output' filesep 'data' filesep 'hairpin_train-' int2str(i)], 'hairpin');
    save(['output' filesep 'data' filesep 'hairpinbracket_train-' int2str(i)], 'hairpinBracket');
    save(['output' filesep 'data' filesep 'hairpin2dfeatures_train-' int2str(i)], 'hairpin2dFeatures');
    
    % create i-th test set
    hairpin = [hairpin_hsa(i_hairpin_hsa_test_ind, :); hairpin_mmu(i_hairpin_mmu_test_ind, :)];
    hairpinBracket = [hairpinBracket_hsa(i_hairpin_hsa_test_ind, :); hairpinBracket_mmu(i_hairpin_mmu_test_ind, :)];
    hairpin2dFeatures = [hairpin2dFeatures_hsa(i_hairpin_hsa_test_ind, :); ...
        hairpin2dFeatures_mmu(i_hairpin_mmu_test_ind, :)];
        
    % save i-th test set
    save(['output' filesep 'data' filesep 'hairpin_test-' int2str(i)], 'hairpin');
    save(['output' filesep 'data' filesep 'hairpinbracket_test-' int2str(i)], 'hairpinBracket');
    save(['output' filesep 'data' filesep 'hairpin2dfeatures_test-' int2str(i)], 'hairpin2dFeatures');    
end

%% Create figures

clc, clear all, close all;

load('param');

load(['output' filesep 'data' filesep 'hairpin_train'], 'hairpin');
load(['output' filesep 'data' filesep 'hairpinbracket_train'], 'hairpinBracket');
load(['output' filesep 'data' filesep 'hairpin2dfeatures_train'], 'hairpin2dFeatures');

table1 = hairpin2dFeatures.miRnaDuplex5pStrand5pEndSeqDistFromTip;
v = min(table1):1:max(table1);
n = hist(table1,v);
ndata = length(table1);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Distance from Tip Position');
xlabel('5p starnd 5p end');
ylabel('Frequency');

table2 = hairpin2dFeatures.miRnaDuplex5pStrand3pEndSeqDistFromTip;
v = min(table2):1:max(table2);
n = hist(table2,v);
ndata = length(table2);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Distance from Tip Position');
xlabel('5p starnd 3p end');
ylabel('Frequency');

table3 = hairpin2dFeatures.miRnaDuplex3pStrand5pEndSeqDistFromTip;
v = min(table3):1:max(table3);
n = hist(table3,v);
ndata = length(table3);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Distance from Tip Position');
xlabel('3p starnd 5p end');
ylabel('Frequency');

table4 = hairpin2dFeatures.miRnaDuplex3pStrand3pEndSeqDistFromTip;
v = min(table4):1:max(table4);
n = hist(table4,v);
ndata = length(table4);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Distance from Tip Position');
xlabel('3p starnd 3p end');
ylabel('Frequency');

%%
table5 = hairpin2dFeatures.miRnaDuplex5pStrandOverhangSequenceLength;
v = min(table5):1:max(table5);
n = hist(table5,v);
ndata = length(table5);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Overhang 5p strand');
%xlabel('');
ylabel('Frequency');

table6 = hairpin2dFeatures.miRnaDuplex3pStrandOverhangSequenceLength;
v = min(table6):1:max(table6);
n = hist(table6,v);
ndata = length(table6);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Overhang 3p strand');
%xlabel('');
ylabel('Frequency');
%% plot distance 5'strand 5p from hairpin start
table6 = hairpin. miRnaDuplex5pStrand5pEndPos;
v = min(table6):1:max(table6);
n = hist(table6,v);
ndata = length(table6);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Distance 5p strand 5p end from hairpin start');
ylabel('Frequency');

%% plot distance 3'strand 3p from hairpin end
table6 = hairpin.sequenceLength - hairpin.miRnaDuplex3pStrand3pEndPos
v = min(table6):1:max(table6);
n = hist(table6,v);
ndata = length(table6);
nnorm = n./ndata;
figure; bar(v, nnorm);
title('Distance 3p strand 3p end from hairpin end');
ylabel('Frequency');
%%
numHairpins = size(hairpin, 1);

% initialize figures
miRnaDuplexStrandSequenceLengthRelFreqFigure = figure();
miRnaDuplexStrandEndSeqDistFromTipRelFreqFigure = figure();

miRnaDuplexStrandSequenceLengthRelFreqAxes = zeros(1, 2);
miRnaDuplexStrandSequenceLengthRelFreqLim2 = zeros(1, 2);

miRnaDuplexStrandEndSeqDistFromTipRelFreqAxes = zeros(2, 2);
miRnaDuplexStrandEndSeqDistFromTipRelFreqYLim2 = zeros(2, 2);
miRnaDuplexStrandEndSeqDistFromTipRelFreqXLim2 = zeros(2, 2);

% create figures

for i = 1:2 % for each strand
%i
    iMiRnaDuplexStrandSequenceLength = hairpin.(sprintf('miRnaDuplex%dpStrandSequenceLength', miRnaDuplexStrand.number(i)));
    iMiRnaDuplexStrandSequenceLengthValues = min(iMiRnaDuplexStrandSequenceLength):max(iMiRnaDuplexStrandSequenceLength);
    iMiRnaDuplexStrandSequenceLengthRelFreq = histc(iMiRnaDuplexStrandSequenceLength, iMiRnaDuplexStrandSequenceLengthValues)/numHairpins;
    
    % -- plot relative frequency of strand length -------------------------

    ah = subplot(1, 2, i, 'Parent', miRnaDuplexStrandSequenceLengthRelFreqFigure);
    miRnaDuplexStrandSequenceLengthRelFreqAxes(i) = ah;
    
    bsh = bar(ah, iMiRnaDuplexStrandSequenceLengthValues, iMiRnaDuplexStrandSequenceLengthRelFreq');

    axis(ah, 'tight');
    iMiRnaDuplexStrandSequenceLengthRelFreqLim = ylim(ah);
    miRnaDuplexStrandSequenceLengthRelFreqLim2(i) = iMiRnaDuplexStrandSequenceLengthRelFreqLim(2);
    ylim(ah, [0 miRnaDuplexStrandSequenceLengthRelFreqLim2(i)]);
    xlim(ah, [iMiRnaDuplexStrandSequenceLengthValues(1) iMiRnaDuplexStrandSequenceLengthValues(end)]);
    title(ah, miRnaDuplexStrand.label{i});
    xlabel(ah, 'Length (nt)');
    ylabel(ah, 'Relative frequency');
    
    for j = 1:2 % for each end

        ijMiRnaDuplexStrandEndSeqDistFromTip = hairpin2dFeatures.(sprintf('miRnaDuplex%dpStrand%dpEndSeqDistFromTip', miRnaDuplexStrand.number(i), strandEnd.number(j)));
        ijMiRnaDuplexStrandEndSeqDistFromTipValues = min(ijMiRnaDuplexStrandEndSeqDistFromTip):max(ijMiRnaDuplexStrandEndSeqDistFromTip);
        ijMiRnaDuplexStrandEndSeqDistFromTipRelFreq = histc(ijMiRnaDuplexStrandEndSeqDistFromTip, ijMiRnaDuplexStrandEndSeqDistFromTipValues)/numHairpins;

        % -- plot relative frequency of strand length ---------------------

        ah = subplot(2, 2, (i - 1)*2 + j, 'Parent', miRnaDuplexStrandEndSeqDistFromTipRelFreqFigure);
        miRnaDuplexStrandEndSeqDistFromTipRelFreqAxes(i, j) = ah;

        bsh = bar(ah, ijMiRnaDuplexStrandEndSeqDistFromTipValues, ijMiRnaDuplexStrandEndSeqDistFromTipRelFreq');

        axis(ah, 'tight');
        ijMiRnaDuplexStrandEndSeqDistFromTipRelFreqYLim = ylim(ah);
        miRnaDuplexStrandEndSeqDistFromTipRelFreqYLim2(i, j) = ijMiRnaDuplexStrandEndSeqDistFromTipRelFreqYLim(2);
        
        ijMiRnaDuplexStrandEndSeqDistFromTipRelFreqXLim = xlim(ah);
        miRnaDuplexStrandEndSeqDistFromTipRelFreqXLim2(i, j) = ijMiRnaDuplexStrandEndSeqDistFromTipRelFreqXLim(2);
  
        xlim(ah, [ijMiRnaDuplexStrandEndSeqDistFromTipValues(1) ijMiRnaDuplexStrandEndSeqDistFromTipValues(end)]);
        title(ah, [miRnaDuplexStrand.label{i} ' ' strandEnd.label{j}]);
        if i == 2, xlabel(ah, 'Distance (nt)'); end
        if j == 1, ylabel(ah, 'Relative frequency'); end
    
    end

end

set(miRnaDuplexStrandSequenceLengthRelFreqAxes, 'YLim', [0 max(miRnaDuplexStrandSequenceLengthRelFreqLim2)]);
print(sprintf('-f%d', miRnaDuplexStrandSequenceLengthRelFreqFigure), '-r300', '-dpdf', ['output' filesep 'figures' filesep 'mirnaduplexstrandseqlengthrelfreq_foreach_mirnaduplexstrand']);

set(miRnaDuplexStrandEndSeqDistFromTipRelFreqAxes(:), 'XLim', [0 max(miRnaDuplexStrandEndSeqDistFromTipRelFreqXLim2(:))]);
set(miRnaDuplexStrandEndSeqDistFromTipRelFreqAxes(:), 'YLim', [0 max(miRnaDuplexStrandEndSeqDistFromTipRelFreqYLim2(:))]);
print(sprintf('-f%d', miRnaDuplexStrandEndSeqDistFromTipRelFreqFigure), '-r300', '-dpdf', ['output' filesep 'figures' filesep 'mirnaduplexstrandendseqdistfromtiprelfreq_foreach_mirnaduplexstrand_strandend'])



%runexp


