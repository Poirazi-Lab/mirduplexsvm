function hairpin2dFeatures = hairpin2dfeatures(hairpin, hairpinBracket)
%HAIRPIN2DFEATURES calculate hairpin 2D features

% get number of hairpins
numHairpins = size(hairpin, 1);

% initialize secondary structure features dataset
hairpin2dFeatures = dataset( ...
    {cell(numHairpins, 1), 'arm5pMatchPos'}, ...
    {cell(numHairpins, 1), 'arm3pMatchPos'}, ... 
    ...
    {false(numHairpins, 1), 'folds'}, ...
    {false(numHairpins, 1), 'multibranched'}, ...
    ...
    {zeros(numHairpins, 1), 'loop5pEndPos'}, ...
    {zeros(numHairpins, 1), 'loop3pEndPos'}, ...
    {zeros(numHairpins, 1), 'loopSequenceLength'}, ...
    ...
    {zeros(numHairpins, 1), 'tipPos'}, ...
    ...
    {zeros(numHairpins, 1), 'miRnaDuplex5pStrand5pEndSeqDistFromTip'}, ...
    {zeros(numHairpins, 1), 'miRnaDuplex5pStrand3pEndSeqDistFromTip'}, ...
    {zeros(numHairpins, 1), 'miRnaDuplex3pStrand5pEndSeqDistFromTip'}, ...
    {zeros(numHairpins, 1), 'miRnaDuplex3pStrand3pEndSeqDistFromTip'}, ...
    ...
    {false(numHairpins, 1), 'miRnaDuplex5pStrand5pEndMatches'}, ...
    {false(numHairpins, 1), 'miRnaDuplex3pStrand5pEndMatches'}, ...
    {zeros(numHairpins, 1), 'matureMatches'}, ...
    ...
    {zeros(numHairpins, 1), 'miRnaDuplex5pStrandOverhangSequenceLength'}, ...
    {zeros(numHairpins, 1), 'miRnaDuplex3pStrandOverhangSequenceLength'} ...
    );
%%
parfor i=1:numHairpins % for each hairpin
    
    fprintf('\nCalculating 2D features of hairpin #%d...\n', i); 
    hairpin2dFeatures(i, :) = hairpin2dfeatures_single(hairpin2dFeatures(i, :), hairpin(i, :), hairpinBracket(i, :));
    
end

end






function hairpin2dFeatures = hairpin2dfeatures_single(hairpin2dFeatures, hairpin, hairpinBracket)
%%
% -- hairpin matching positions -------------------------------------------

% find matches per arm
hairpin2dFeatures.arm5pMatchPos{1} = strfind(hairpinBracket.bracket{1}, '('); 
hairpin2dFeatures.arm3pMatchPos{1} = strfind(hairpinBracket.bracket{1}, ')');

if isempty(hairpin2dFeatures.arm5pMatchPos{1}) % no folding
    
    disp('Non-folding hairpin');
    return;
    
end

hairpin2dFeatures.folds = true;
%%
% -- terminal loop --------------------------------------------------------

% loop start = just after last 5' match
loop5pEndPos = hairpin2dFeatures.arm5pMatchPos{1}(end) + 1;

% loop end = just before first 3' match
loop3pEndPos = hairpin2dFeatures.arm3pMatchPos{1}(1) - 1;

if loop5pEndPos > loop3pEndPos

    % multibranched hairpin
    disp('Multibranch hairpin');
    hairpin2dFeatures.multibranched = true;
    return;

end

hairpin2dFeatures.loop5pEndPos = loop5pEndPos;
hairpin2dFeatures.loop3pEndPos = loop3pEndPos;
hairpin2dFeatures.loopSequenceLength = loop3pEndPos - loop5pEndPos + 1; % length of the loop

% set tip position
hairpin2dFeatures.tipPos = hairpin2dFeatures.loop5pEndPos - 1 + ceil(hairpin2dFeatures.loopSequenceLength/2);

hairpin2dFeatures.miRnaDuplex5pStrand5pEndSeqDistFromTip = hairpin2dFeatures.tipPos - hairpin.miRnaDuplex5pStrand5pEndPos;
hairpin2dFeatures.miRnaDuplex5pStrand3pEndSeqDistFromTip = hairpin2dFeatures.tipPos - hairpin.miRnaDuplex5pStrand3pEndPos;
hairpin2dFeatures.miRnaDuplex3pStrand5pEndSeqDistFromTip = hairpin.miRnaDuplex3pStrand5pEndPos - hairpin2dFeatures.tipPos;
hairpin2dFeatures.miRnaDuplex3pStrand3pEndSeqDistFromTip = hairpin.miRnaDuplex3pStrand3pEndPos - hairpin2dFeatures.tipPos;

% -- 5' ends match -------------------------------------------------------- 

hairpin2dFeatures.miRnaDuplex5pStrand5pEndMatches = hairpinBracket.bracket{1}(hairpin.miRnaDuplex5pStrand5pEndPos) ~= '.';
hairpin2dFeatures.miRnaDuplex3pStrand5pEndMatches = hairpinBracket.bracket{1}(hairpin.miRnaDuplex3pStrand5pEndPos) ~= '.';

% -- overhangs ------------------------------------------------------------

duplex = [hairpin.miRnaDuplex5pStrand5pEndPos hairpin.miRnaDuplex5pStrand3pEndPos hairpin.miRnaDuplex3pStrand5pEndPos hairpin.miRnaDuplex3pStrand3pEndPos];

hairpin2dFeatures.miRnaDuplex5pStrandOverhangSequenceLength = org.mensxmachina.mirna.duplex5pstrandoverhangseqlengthq(hairpinBracket.bracket{1}, hairpin2dFeatures.arm5pMatchPos{1}, hairpin2dFeatures.arm3pMatchPos{1}, duplex);
hairpin2dFeatures.miRnaDuplex3pStrandOverhangSequenceLength = org.mensxmachina.mirna.duplex3pstrandoverhangseqlengthq(hairpinBracket.bracket{1}, hairpin2dFeatures.arm5pMatchPos{1}, hairpin2dFeatures.arm3pMatchPos{1}, duplex);

%%
% -- matureMatches----------------------------------------------------------

duplex = [hairpin.miRnaDuplex5pStrand5pEndPos(1) , ...
    hairpin.miRnaDuplex5pStrand3pEndPos(1) , hairpin.miRnaDuplex3pStrand5pEndPos(1),...
    hairpin.miRnaDuplex3pStrand3pEndPos(1)];

hairpin2dFeatures.matureMatches = ...
    org.mensxmachina.mirna.findMatureMatches(hairpinBracket.bracket{1}, ...
    hairpin2dFeatures.arm5pMatchPos{1}, hairpin2dFeatures.arm3pMatchPos{1}, ...
    duplex);





end




