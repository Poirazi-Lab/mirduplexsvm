function miRnaDuplex = mirnaduplextrivialfindq(model, hairpinSeq, hairpinBracket)
%MIRNADUPLEXTRIVIALFINDQ Apply trivial miRNA:miRNA*-duplex  finder

import org.mensxmachina.mirna.*;

numHairpins = length(hairpinSeq);

hairpinSeqLength = zeros(numHairpins, 1);
hairpinTipPos = zeros(numHairpins, 1);

for i = 1:numHairpins % for each hairpin

    fprintf('\nFolding hairpin...\n');
    
    % set hairpin length
    hairpinSeqLength(i) = length(hairpinSeq{i});

    % find matches per arm
    hairpin5pArmMatchPos = strfind(hairpinBracket{i}, '('); 
    hairpin3pArmMatchPos = strfind(hairpinBracket{i}, ')');

    % find tip position
    hairpinTipPos(i) = org.mensxmachina.mirna.findhairpintipq(hairpin5pArmMatchPos, hairpin3pArmMatchPos);
    
end

% create duplex with strands of mean length and mean distance from hairpin tip

miRnaDuplex = zeros(numHairpins, 4);

miRnaDuplex(:, 1) = max(hairpinTipPos - model.meanMiRnaDuplex5pStrand5pEndSeqDistFromTip, 1);
miRnaDuplex(:, 2) = max(hairpinTipPos - model.meanMiRnaDuplex5pStrand3pEndSeqDistFromTip, 1);
miRnaDuplex(:, 3) = min(hairpinTipPos + model.meanMiRnaDuplex3pStrand5pEndSeqDistFromTip, hairpinSeqLength);
miRnaDuplex(:, 4) = min(hairpinTipPos + model.meanMiRnaDuplex3pStrand3pEndSeqDistFromTip, hairpinSeqLength);

end