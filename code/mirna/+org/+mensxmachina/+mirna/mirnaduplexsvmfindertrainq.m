function model = mirnaduplexsvmfindertrainq(hairpinSeq, hairpinBracket, ...
    miRnaDuplex, miRnaDuplexFinderTrainConfig)
%MIRNADUPLEXSVMFINDERTRAINQ Train miRNA:miRNA*-duplex SVM finder

Param = miRnaDuplexFinderTrainConfig.trainParam{1};
import org.mensxmachina.mirna.*;

% calculate miRNA:miRNA*-duplex sequence lengths
miRnaDuplex5pStrandSeqLength = miRnaDuplex(:, 2) - miRnaDuplex(:, 1) + 1;
miRnaDuplex3pStrandSeqLength = miRnaDuplex(:, 4) - miRnaDuplex(:, 3) + 1;

% get number of hairpins
numHairpins = length(hairpinSeq);

% initialize overhang sequence lengths
miRnaDuplexOverhang = zeros(numHairpins, 2);

% initialize hairpin secondary structure features
hairpin5pArmMatchPos = cell(numHairpins, 1);
hairpin3pArmMatchPos = cell(numHairpins, 1);
hairpinLoop5pEndPos = zeros(numHairpins, 1);
hairpinLoop3pEndPos = zeros(numHairpins, 1);
hairpinLoopSeqLength = zeros(numHairpins, 1);
hairpinTipPos = zeros(numHairpins, 1);

fprintf('\nFolding hairpins...\n');

for i = 1:numHairpins % for each hairpin
    
    % find matches per arm
    hairpin5pArmMatchPos{i} = strfind(hairpinBracket{i}, '('); 
    hairpin3pArmMatchPos{i} = strfind(hairpinBracket{i}, ')');

    % loop start = just after last 5' match
     hairpinLoop5pEndPos(i) = hairpin5pArmMatchPos{i}(end) + 1;

    % loop end = just before first 3' match
     hairpinLoop3pEndPos(i) = hairpin3pArmMatchPos{i}(1) - 1;

    % calculate loop sequence length
    hairpinLoopSeqLength(i) = hairpinLoop3pEndPos(i) - hairpinLoop5pEndPos(i) + 1; 
    
    % find tip position
    hairpinTipPos(i) = hairpinLoop5pEndPos(i) - 1 + ceil(hairpinLoopSeqLength(i)/2);
    
    % calculate overhang lengths
    miRnaDuplexOverhang(i, 1) = org.mensxmachina.mirna.duplex5pstrandoverhangseqlengthq(hairpinBracket{i},...
        hairpin5pArmMatchPos{i}, hairpin3pArmMatchPos{i}, miRnaDuplex(i, :));
    miRnaDuplexOverhang(i, 2) = org.mensxmachina.mirna.duplex3pstrandoverhangseqlengthq(hairpinBracket{i},...
        hairpin5pArmMatchPos{i}, hairpin3pArmMatchPos{i}, miRnaDuplex(i, :));

end

% calculate miRNA:miRNA*-duplex sequence distance from hairpin tip
miRnaDuplex5pStrand5pEndSeqDistFromTip = hairpinTipPos - miRnaDuplex(:, 1);
miRnaDuplex5pStrand3pEndSeqDistFromTip = hairpinTipPos - miRnaDuplex(:, 2);
miRnaDuplex3pStrand5pEndSeqDistFromTip = miRnaDuplex(:, 3) - hairpinTipPos;
miRnaDuplex3pStrand3pEndSeqDistFromTip = miRnaDuplex(:, 4) - hairpinTipPos;

candidateMiRnaDuplexModel = struct(...
    'strand5pSeqLengthLim', [min(miRnaDuplex5pStrandSeqLength) ...
    max(miRnaDuplex5pStrandSeqLength)], ...
    'strand3pSeqLengthLim', [min(miRnaDuplex3pStrandSeqLength) ...
    max(miRnaDuplex3pStrandSeqLength)], ...
    'strand5pOverhangSeqLengthLim', [min(miRnaDuplexOverhang(:, 1)) ...
    max(miRnaDuplexOverhang(:, 1))], ...
    'strand3pOverhangSeqLengthLim', [min(miRnaDuplexOverhang(:, 2)) ...
    max(miRnaDuplexOverhang(:, 2))], ...
    'Strand5pEnd5pSeqDistFromTip', [min(miRnaDuplex5pStrand5pEndSeqDistFromTip)...
    max(miRnaDuplex5pStrand5pEndSeqDistFromTip)], ...
    'Strand5pEnd3pSeqDistFromTip', [min(miRnaDuplex5pStrand3pEndSeqDistFromTip)...
    max(miRnaDuplex5pStrand3pEndSeqDistFromTip)], ...
    'Strand3pEnd5pSeqDistFromTip', [min(miRnaDuplex3pStrand5pEndSeqDistFromTip)...
    max(miRnaDuplex3pStrand5pEndSeqDistFromTip)], ...
    'Strand3pEnd3pSeqDistFromTip', [min(miRnaDuplex3pStrand3pEndSeqDistFromTip)...
    max(miRnaDuplex3pStrand3pEndSeqDistFromTip)]);

candidateMiRnaDuplexParam = struct('Verbose', Param.Verbose);

trainCandidateMiRnaDuplexParam = struct('Ratio', Param.Ratio, 'Verbose',...
    Param.Verbose);

hairpinTrainCandidateMiRnaDuplex = cell(numHairpins, 1);
hairpinTrainCandidateMiRnaDuplexOverhang = cell(numHairpins, 1);
hairpinTrainCandidateMiRnaDuplexIsMiRnaDuplex = cell(numHairpins, 1);
fprintf('\nGenerating and selecting candidate miRNA:miRNA* duplexes...\n');

    candidateMiRnaDuplexCacheFilename = Param.CandidateMiRnaDuplexCacheFilename;
    found = zeros(size(candidateMiRnaDuplexCacheFilename,1),1) ;
    

    for i = 1:numHairpins % for each hairpin
 
        [candidateMiRnaDuplex candidateMiRnaDuplexOverhang] = ...
            candidatemirnaduplexwrapper(candidateMiRnaDuplexModel, hairpinBracket{i}, ...
            hairpin5pArmMatchPos{i}, hairpin3pArmMatchPos{i}, hairpinTipPos(i), ...
            candidateMiRnaDuplexParam, candidateMiRnaDuplexCacheFilename{i});


        [hairpinTrainCandidateMiRnaDuplex{i} hairpinTrainCandidateMiRnaDuplexOverhang{i}...
            hairpinTrainCandidateMiRnaDuplexIsMiRnaDuplex{i}, FOUND] = ...
            traincandidatemirnaduplexq(candidateMiRnaDuplex, candidateMiRnaDuplexOverhang,...
            miRnaDuplex(i, :), miRnaDuplexOverhang(i, :), trainCandidateMiRnaDuplexParam);
                        
        if FOUND == 1
            found(i,1) = 1;
        else
            found(i,1) = 1000;
        end 
        
    end
filename = strcat(miRnaDuplexFinderTrainConfig.sampleName{1},'found');
save(['output' filesep 'data' filesep filename], 'found');


% merge duplexes
numTrainCandidateMiRnaDuplexes = numHairpins*(1 + Param.Ratio);
trainCandidateMiRnaHairpinSeq = cell(numTrainCandidateMiRnaDuplexes, 1);
trainCandidateMiRnaHairpinTipPos = zeros(numTrainCandidateMiRnaDuplexes, 1);
trainCandidateMiRnaHairpinLoopSeqLength = zeros(numTrainCandidateMiRnaDuplexes, 1);
trainCandidateMiRnaDuplex = zeros(numTrainCandidateMiRnaDuplexes, 4);
trainCandidateMiRnaDuplexOverhang = zeros(numTrainCandidateMiRnaDuplexes, 2);
trainCandidateMiRnaDuplexSeqDistFromHairpinLoopTip = zeros(numTrainCandidateMiRnaDuplexes, 4); 
trainCandidateMiRnaDuplexIsMiRnaDuplex = false(numTrainCandidateMiRnaDuplexes, 1);

for i = 1:numHairpins % for each hairpin
 
    % create training candidate miRNA:miRNA* duplex indices
    ind = ((i - 1)*(1 + Param.Ratio) + 1):(i*(1 + Param.Ratio));
    
    % replicate hairpin data for all training candidate miRNA:miRNA* duplexes
    trainCandidateMiRnaHairpinSeq(ind) = ...
        repmat(hairpinSeq(i), 1 + Param.Ratio, 1);
    %trainCandidateMiRnaHairpinThermo(ind) = repmat(hairpinThermo(i), 1 + Param.Ratio, 1);
    trainCandidateMiRnaHairpinTipPos(ind) = ...
        repmat(hairpinTipPos(i), 1 + Param.Ratio, 1);
    trainCandidateMiRnaHairpinLoopSeqLength(ind) = ...
        repmat(hairpinLoopSeqLength(i), 1 + Param.Ratio, 1);
    
    % get training candidate miRNA:miRNA* duplexes
    trainCandidateMiRnaDuplex(ind, :) = hairpinTrainCandidateMiRnaDuplex{i};
    trainCandidateMiRnaDuplexOverhang(ind, :) = hairpinTrainCandidateMiRnaDuplexOverhang{i};
    trainCandidateMiRnaDuplexIsMiRnaDuplex(ind) = hairpinTrainCandidateMiRnaDuplexIsMiRnaDuplex{i};

    trainCandidateMiRnaDuplexSeqDistFromHairpinLoopTip(ind, 1) = ...
        hairpinTipPos(i) - hairpinTrainCandidateMiRnaDuplex{i}(:, 1);
    trainCandidateMiRnaDuplexSeqDistFromHairpinLoopTip(ind, 2) = ...
        hairpinTipPos(i) - hairpinTrainCandidateMiRnaDuplex{i}(:, 2);
   trainCandidateMiRnaDuplexSeqDistFromHairpinLoopTip(ind, 3) = ...
        hairpinTrainCandidateMiRnaDuplex{i}(:, 3) - hairpinTipPos(i);    
    trainCandidateMiRnaDuplexSeqDistFromHairpinLoopTip(ind, 4) = ...
        hairpinTrainCandidateMiRnaDuplex{i}(:, 4) - hairpinTipPos(i);    
    
    hairpinTrainCandidateMiRnaDuplex{i} = [];
    hairpinTrainCandidateMiRnaDuplexIsMiRnaDuplex{i} = [];
    
end

clear('hairpinTrainCandidateMiRnaDuplex', 'hairpinTrainCandidateMiRnaDuplexIsMiRnaDuplex');


miRnaDuplexSvmVectorModel = candidateMiRnaDuplexModel; 

miRnaDuplexSvmVectorModel.miRnaDuplex5pStrand5pEndSeqDistFromTipLim = ...
    [min(miRnaDuplex5pStrand5pEndSeqDistFromTip) max(miRnaDuplex5pStrand5pEndSeqDistFromTip)];

miRnaDuplexSvmVectorModel.miRnaDuplex5pStrand3pEndSeqDistFromTipLim = ...
    [min(miRnaDuplex5pStrand3pEndSeqDistFromTip) max(miRnaDuplex5pStrand3pEndSeqDistFromTip)];

miRnaDuplexSvmVectorModel.miRnaDuplex3pStrand5pEndSeqDistFromTipLim = ... 
    [min(miRnaDuplex3pStrand5pEndSeqDistFromTip) max(miRnaDuplex3pStrand5pEndSeqDistFromTip)];

miRnaDuplexSvmVectorModel.miRnaDuplex3pStrand3pEndSeqDistFromTipLim = ... 
    [min(miRnaDuplex3pStrand3pEndSeqDistFromTip) max(miRnaDuplex3pStrand3pEndSeqDistFromTip)];

miRnaDuplexSvmVectorModel.flankSeqLength = Param.FlankingSequenceLength;

fprintf('\nCreating SVM input...\n');

numSeqCols = 4*(4*miRnaDuplexSvmVectorModel.flankSeqLength + ... 
    miRnaDuplexSvmVectorModel.strand5pSeqLengthLim(2) + ...
    miRnaDuplexSvmVectorModel.strand3pSeqLengthLim(2));

numCols = numSeqCols;

x = zeros(numTrainCandidateMiRnaDuplexes, numCols);

for i = 1:numTrainCandidateMiRnaDuplexes 
    
    x(i, :) = mirnaDuplexSVMvectorSEQ_THERMq( ...
        miRnaDuplexSvmVectorModel, ...
        trainCandidateMiRnaHairpinSeq{i}, ...
        trainCandidateMiRnaHairpinTipPos(i), ...
        trainCandidateMiRnaDuplex(i, :));

end
group = mirnaduplexsvmgroupq(trainCandidateMiRnaDuplexIsMiRnaDuplex);

fprintf('\nTraining SVM...\n');
% it works with libSVM 
miRnaDuplexSvmModel = mirnaduplexsvmtrainq(x, group, Param.SvmTrainParam);

% combine models
model = struct( ...
    'candidateMiRnaDuplexModel', candidateMiRnaDuplexModel, ...
    'miRnaDuplexSvmVectorModel', miRnaDuplexSvmVectorModel, ...
    'miRnaDuplexSvmModel', miRnaDuplexSvmModel); 
    
 end

