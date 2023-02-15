function [duplex overhang] = candidatemirnaduplexwrapperF(model, hairpinBracket, ...
    hairpin5pArmMatchPos, hairpin3pArmMatchPos, hairpinTipPos, ...
    candidateMiRnaDuplexParam, candidateMiRnaDuplexCacheFilename)

candidateMiRnaDuplexCacheFilename = sprintf(sprintf('%s.mat', candidateMiRnaDuplexCacheFilename));

if exist(candidateMiRnaDuplexCacheFilename, 'file')
    
    fprintf('\nLoading candidate miRNA:miRNA* duplexes from cache...\n'); 
    load(candidateMiRnaDuplexCacheFilename, 'duplex', 'overhang');
    
else
    [duplex overhang] = org.mensxmachina.mirna.candidatemirnaduplexNEWq(model, ...
        hairpinBracket, hairpin5pArmMatchPos, hairpin3pArmMatchPos, hairpinTipPos, ...
        candidateMiRnaDuplexParam);
    
    save(candidateMiRnaDuplexCacheFilename, 'duplex', 'overhang');
end

end