function [hairpinBracket hairpinEntropy] = hairpinfold(hairpin)
%HAIRPINBRACKET fold hairpins

% get number of hairpins
numHairpins = size(hairpin, 1);

% initialize secondary structure dataset
hairpinBracket = dataset( ...
    {cell(numHairpins, 1), 'bracket'} ...
    );

hairpinEntropy = dataset( ...
    {cell(numHairpins, 1), 'thermo'} ...
    );
for i = 1:numHairpins % for each hairpin
    
    fprintf('\nFolding hairpin #%d...\n', i); 
    [hairpinBracket(i, :), hairpinEntropy(i,:)] = ...
        hairpinfold_single(hairpinBracket(i, :),hairpinEntropy(i,:), hairpin(i, :));
    
end

end

function [hairpinBracket, hairpinEntropy] = ...
    hairpinfold_single(hairpinBracket,hairpinEntropy, hairpin)

% for short sequences
cacheFile = ['cache' filesep 'hairpinbracket_' lower(hairpin.Properties.ObsNames{1})];

if exist(sprintf('%s.mat', cacheFile), 'file')

    fprintf('\nLoading secondary structure from cache...\n');

    load(cacheFile, 'bracket');
else

    [bracket] = org.mensxmachina.mirna.hairpinfoldq(hairpin.sequence{1});
    
    save(cacheFile, 'bracket');
    
end
hairpinBracket.bracket{1} = bracket;

end




