function hairpinThermo = hairpinfold_UNAfold(hairpin)
%HAIRPINBRACKET fold hairpins

% get number of hairpins
numHairpins = size(hairpin, 1);

% initialize secondary structure dataset
hairpinThermo = dataset( ...
    {cell(numHairpins, 1), 'thermo'} ...
    );

parfor i = 1:numHairpins % for each hairpin
    
    fprintf('\nFolding hairpin UNAfold #%d...\n', i); 
    hairpinThermo(i, :) = hairpinfold_single(hairpinThermo(i, :), hairpin(i, :));
    
end

end

function hairpinThermo = hairpinfold_single(hairpinThermo, hairpin)

cacheFile = ['cache' filesep 'hairpinThermo_hairpin_' ...
    lower(hairpin.Properties.ObsNames{1})];

if exist(sprintf('%s.mat', cacheFile), 'file')

    fprintf('\nLoading thermo data from cache...\n');

    load(cacheFile, 'thermo');

else

    thermo = org.mensxmachina.mirna.thermodynamics.UNAfoldCaller(hairpin.sequence{1});
    
    save(cacheFile, 'thermo');

end

hairpinThermo.thermo{1} = thermo;

end



