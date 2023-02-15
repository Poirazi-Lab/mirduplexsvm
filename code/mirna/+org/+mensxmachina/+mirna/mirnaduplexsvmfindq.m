function miRnaDuplex = mirnaduplexsvmfindq(model, hairpinSeq, ...
    hairpinBracket, Param)
%MIRNADUPLEXSVMFINDQ Apply miRNA:miRNA*-miRnaDuplex SVM finder

import org.mensxmachina.mirna.*;

numHairpins = length(hairpinSeq);

miRnaDuplex = zeros(numHairpins, 4);

numSeqCols = 4*(4*model.miRnaDuplexSvmVectorModel.flankSeqLength + ...
    model.miRnaDuplexSvmVectorModel.strand5pSeqLengthLim(2) + ...
    model.miRnaDuplexSvmVectorModel.strand3pSeqLengthLim(2));


numCols = numSeqCols;

candidateMiRnaDuplexModel = model.candidateMiRnaDuplexModel;
Verbose = Param.Verbose;
CandidateMiRnaDuplexCacheFilename = Param.CandidateMiRnaDuplexCacheFilename;
miRnaDuplexSvmVectorModel = model.miRnaDuplexSvmVectorModel;
miRnaDuplexSvmModel = model.miRnaDuplexSvmModel;

for i = 1:numHairpins % for each hairpin

    fprintf('\nFolding hairpin...\n');

    % find matches per arm
    hairpin5pArmMatchPos = strfind(hairpinBracket{i}, '('); 
    hairpin3pArmMatchPos = strfind(hairpinBracket{i}, ')');
    
    % loop start = just after last 5' match
    hairpinLoop5pEndPos = hairpin5pArmMatchPos(end) + 1;

    % loop end = just before first 3' match
    hairpinLoop3pEndPos = hairpin3pArmMatchPos(1) - 1;

    % calculate loop sequence length
    hairpinLoopSeqLength = hairpinLoop3pEndPos - hairpinLoop5pEndPos + 1; 

    % find tip position
    hairpinTipPos = hairpinLoop5pEndPos - 1 + ceil(hairpinLoopSeqLength/2);
    
    candidateMiRnaDuplexParam = struct('Verbose', Verbose);

    fprintf('\nGenerating candidate miRNA:miRNA* duplexes...\n');

    [candidateMiRnaDuplex candidateMiRnaDuplexOverhang] = ...
        candidatemirnaduplexwrapper(candidateMiRnaDuplexModel, ...
        hairpinBracket{i}, hairpin5pArmMatchPos, hairpin3pArmMatchPos, ...
        hairpinTipPos, candidateMiRnaDuplexParam, ...
        CandidateMiRnaDuplexCacheFilename{i});

    numTrainCandidateMiRnaDuplexes = size(candidateMiRnaDuplex, 1);
    
    % calculate candidate miRNA:miRNA* duplex sequence distance from loop Tip
    candidateMiRnaDuplexDistTip = zeros(numTrainCandidateMiRnaDuplexes, 4);
    candidateMiRnaDuplexDistTip(:, 1) = hairpinTipPos - candidateMiRnaDuplex(:, 1);
    candidateMiRnaDuplexDistTip(:, 2) = hairpinTipPos - candidateMiRnaDuplex(:, 2);
    candidateMiRnaDuplexDistTip(:, 3) = candidateMiRnaDuplex(:, 3) - hairpinTipPos;
    candidateMiRnaDuplexDistTip(:, 4) = candidateMiRnaDuplex(:, 4) - hairpinTipPos;

    fprintf('\nCreating SVM input...\n');

    x = zeros(numTrainCandidateMiRnaDuplexes, numCols);

    % for each training candidate miRNA:miRNA* duplex 
    for j = 1:numTrainCandidateMiRnaDuplexes 
        x(j, :) = mirnaDuplexSVMvectorSEQ_THERMq(...
            miRnaDuplexSvmVectorModel,...
            hairpinSeq{i}, ... 
            hairpinTipPos, ...
            candidateMiRnaDuplex(j, :)); 
    end

    fprintf('\nClassifying using SVM...\n');

    [a, candidateMiRnaDuplexScore] = org.mensxmachina.mirna.mirnaduplexsvmclassifyq(miRnaDuplexSvmModel, x);

    % find max score
    [b, ind] = max(candidateMiRnaDuplexScore);

    % get candidate miRNA:miRNA* duplex with max score
    miRnaDuplex(i, :) = candidateMiRnaDuplex(ind, :);
    
end

 end

