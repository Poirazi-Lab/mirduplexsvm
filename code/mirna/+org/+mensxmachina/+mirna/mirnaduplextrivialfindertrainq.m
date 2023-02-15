function model = mirnaduplextrivialfindertrainq(hairpinSeq, hairpinBracket, miRnaDuplex)
%MIRNADUPLEXTRIVIALFINDERTRAINQ Train trivial miRNA:miRNA*-duplex SVM finder

import org.mensxmachina.mirna.*;

% get number of hairpins
numHairpins = length(hairpinSeq);

% initialize hairpin tip positions
hairpinTipPos = zeros(numHairpins, 1);

fprintf('\nFolding hairpins...\n');

for i = 1:numHairpins % for each hairpin
    
    % find matches per arm
    hairpin5pArmMatchPos = strfind(hairpinBracket{i}, '('); 
    hairpin3pArmMatchPos = strfind(hairpinBracket{i}, ')');

    % find tip position
    hairpinTipPos(i) = org.mensxmachina.mirna.findhairpintipq(hairpin5pArmMatchPos, hairpin3pArmMatchPos);
    
end

% calculate miRNA:miRNA*-duplex sequence distance from hairpin tip
miRnaDuplex5pStrand5pEndSeqDistFromTip = hairpinTipPos - miRnaDuplex(:, 1);
miRnaDuplex5pStrand3pEndSeqDistFromTip = hairpinTipPos - miRnaDuplex(:, 2);
miRnaDuplex3pStrand5pEndSeqDistFromTip = miRnaDuplex(:, 3) - hairpinTipPos;
miRnaDuplex3pStrand3pEndSeqDistFromTip = miRnaDuplex(:, 4) - hairpinTipPos;

% learn miRNA:miRNA*-duplex model
model = struct(...
    'meanMiRnaDuplex5pStrand5pEndSeqDistFromTip', round(mean(miRnaDuplex5pStrand5pEndSeqDistFromTip)), ...
    'meanMiRnaDuplex5pStrand3pEndSeqDistFromTip', round(mean(miRnaDuplex5pStrand3pEndSeqDistFromTip)), ...
    'meanMiRnaDuplex3pStrand5pEndSeqDistFromTip', round(mean(miRnaDuplex3pStrand5pEndSeqDistFromTip)), ...
    'meanMiRnaDuplex3pStrand3pEndSeqDistFromTip', round(mean(miRnaDuplex3pStrand3pEndSeqDistFromTip)) ...
);

end