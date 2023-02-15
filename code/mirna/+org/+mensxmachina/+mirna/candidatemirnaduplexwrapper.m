
function [duplex overhang] = candidatemirnaduplexwrapper(candidateMiRnaDuplexModel, ...
    hairpinBracket, hairpin5pArmMatchPos, hairpin3pArmMatchPos, ...
    hairpinTipPos, candidateMiRnaDuplexParam, candidateMiRnaDuplexCacheFilename)
% CHECK IF THE FILE WITH THE CANDIDIATE DUPLEXES ALREADY EXISTS
% IF TRUE LOAD THE FILE


%   %% to run it alone
%   candidateMiRnaDuplexModel; 
%   hairpinBracket = hairpinBracket{1};
%   hairpin5pArmMatchPos = hairpin5pArmMatchPos{1}; 
%   hairpin3pArmMatchPos = hairpin3pArmMatchPos{1};
%   hairpinTipPos = hairpinTipPos(1);
%   candidateMiRnaDuplexParam;
%   candidateMiRnaDuplexCacheFilename = candidateMiRnaDuplexCacheFilename{1};
%  %%

% candidateMiRnaDuplexModel = candidateMiRnaDuplexModel{1};
% candidateMiRnaDuplexParam = candidateMiRnaDuplexParam{1};
  
candidateMiRnaDuplexCacheFilename = sprintf(sprintf('%s.mat', ...
    candidateMiRnaDuplexCacheFilename));

if exist(candidateMiRnaDuplexCacheFilename, 'file')
    
    fprintf('\nLoading candidate miRNA:miRNA* duplexes from cache...\n');
    
    load(candidateMiRnaDuplexCacheFilename, 'duplex', 'overhang');
    
else
    [duplex overhang] = org.mensxmachina.mirna.candidatemirnaduplexNEWq...
        (candidateMiRnaDuplexModel, hairpinBracket, hairpin5pArmMatchPos, ...
        hairpin3pArmMatchPos, hairpinTipPos, candidateMiRnaDuplexParam);
    save(candidateMiRnaDuplexCacheFilename, 'duplex', 'overhang');
end

end
