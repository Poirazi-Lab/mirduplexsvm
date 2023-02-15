function hairpin = hairpinembl2dataset(Embl)
%HAIRPINEMBL2DATASET Convert miRBase EMBL struct to dataset

% get number of hairpins
numHairpins = size(Embl, 2);

% initialize hairpin dataset
hairpin = dataset( ...
    {cell(numHairpins, 1), 'speciesName'}, ...
    {cell(numHairpins, 1), 'id'}, ...
    ...
    {cell(numHairpins, 1), 'sequence'}, ...
    {zeros(numHairpins, 1), 'sequenceLength'}, ...
    ...
    {zeros(numHairpins, 1), 'aCount'}, ...
    {zeros(numHairpins, 1), 'cCount'}, ...
    {zeros(numHairpins, 1), 'gCount'}, ...
    {zeros(numHairpins, 1), 'uCount'}, ...
    {zeros(numHairpins, 1), 'otherBaseCount'}, ...
    ...
    {zeros(numHairpins, 1), 'numMatureProductsKnown'}, ...
    ...
    {false(numHairpins, 1), 'miRnaDuplexKnown'}, ...
    ...
    {NaN(numHairpins, 1), 'miRnaDuplex5pStrand5pEndPos'}, ...
    {NaN(numHairpins, 1), 'miRnaDuplex5pStrand3pEndPos'}, ...
    {cell(numHairpins, 1), 'miRnaDuplex5pStrandSequence'}, ...
    {NaN(numHairpins, 1), 'miRnaDuplex5pStrandSequenceLength'}, ...
    ...
    {NaN(numHairpins, 1), 'miRnaDuplex3pStrand5pEndPos'}, ...
    {NaN(numHairpins, 1), 'miRnaDuplex3pStrand3pEndPos'}, ...
    {cell(numHairpins, 1), 'miRnaDuplex3pStrandSequence'}, ...
    {NaN(numHairpins, 1), 'miRnaDuplex3pStrandSequenceLength'}, ...
    'ObsNames', arrayfun(@(i) sprintf('temp%d', i), 1:numHairpins, 'UniformOutput', false) ...
    );

for i = 1:numHairpins % for each hairpin
    
    fprintf('\nProcessing hairpin #%d...\n', i);
    
    % set hairpin name
    hairpin.Properties.ObsNames{i} = Embl(i).Identification.EntryName;
    
    % set hairpin ID
    hairpin.id{i} = Embl(i).Accession;
    
    % set species name
    hairpin.speciesName{i} = lower(Embl(i).Identification.Division);
    
    % set hairpin sequence
    hairpin.sequence{i} = upper(Embl(i).Sequence);
    
    % set hairpin length
    hairpin.sequenceLength(i) = length(Embl(i).Sequence);
    
    % set hairpin base counts
    hairpin.aCount(i) = Embl(i).BaseCount.A;
    assert(sum(hairpin.sequence{i} == 'A') == hairpin.aCount(i));
    hairpin.cCount(i) = Embl(i).BaseCount.C;
    assert(sum(hairpin.sequence{i} == 'C') == hairpin.cCount(i));
    hairpin.gCount(i) = Embl(i).BaseCount.G;
    assert(sum(hairpin.sequence{i} == 'G') == hairpin.gCount(i));
    hairpin.uCount(i) = sum(hairpin.sequence{i} == 'U');
    hairpin.otherBaseCount(i) = sum(~ismember(hairpin.sequence{i}, 'ACGU'));
    assert(hairpin.otherBaseCount(i) - sum(hairpin.sequence{i} == 'T') + hairpin.uCount(i) == Embl(i).BaseCount.Other);
    
    % parse feature table 
    MiRNA = featuresparse(Embl(i).Feature, 'feature', 'miRNA');
    
    hairpin.numMatureProductsKnown(i) = length(MiRNA); % set number of mature products in MiRBase
    
    if hairpin.numMatureProductsKnown(i) == 2 % if there are 2 mature products
        
        % assert first listed product is upstream from second listed product
        assert(MiRNA(1).Indices(1) <= MiRNA(2).Indices(1));
        
        if MiRNA(2).Indices(1) <= MiRNA(1).Indices(2) % product sequences are overlapping
            continue;
        end
        
        % warning: may still be not a duplex!
        
        % set that duplex is known
        hairpin.miRnaDuplexKnown(i) = true;
        
        % set duplex strand indices and lengths
        hairpin.miRnaDuplex5pStrand5pEndPos(i) = MiRNA(1).Indices(1);
        hairpin.miRnaDuplex5pStrand3pEndPos(i) = MiRNA(1).Indices(2);
        hairpin.miRnaDuplex5pStrandSequence{i} = hairpin.sequence{i}(hairpin.miRnaDuplex5pStrand5pEndPos(i):hairpin.miRnaDuplex5pStrand3pEndPos(i));
        hairpin.miRnaDuplex5pStrandSequenceLength(i) = length(hairpin.miRnaDuplex5pStrandSequence{i});

        hairpin.miRnaDuplex3pStrand5pEndPos(i) = MiRNA(2).Indices(1);
        hairpin.miRnaDuplex3pStrand3pEndPos(i) = MiRNA(2).Indices(2);
        hairpin.miRnaDuplex3pStrandSequence{i} = hairpin.sequence{i}(hairpin.miRnaDuplex3pStrand5pEndPos(i):hairpin.miRnaDuplex3pStrand3pEndPos(i));
        hairpin.miRnaDuplex3pStrandSequenceLength(i) = length(hairpin.miRnaDuplex3pStrandSequence{i});

    end
    
end

end