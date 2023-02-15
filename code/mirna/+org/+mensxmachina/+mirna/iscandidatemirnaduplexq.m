function [tf overhang] = iscandidatemirnaduplexq(model, hairpinBracket, ...
    hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex, Param)
%ISCANDIDATEMIRNADUPLEXQ Summary of this function goes here
%   Detailed explanation goes here

tf = false;
overhang = [NaN NaN];

if ~org.mensxmachina.mirna.isduplexq(hairpinBracket, hairpin5pArmMatchPos, ...
        hairpin3pArmMatchPos, duplex)                  
    if Param.Verbose
        fprintf(' Not a duplex.');
    end
    return;
end

% calculate 5' strand 3' overhang
overhang(1) = org.mensxmachina.mirna.duplex5pstrandoverhangseqlengthq...
    (hairpinBracket, hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex);

if overhang(1) < model.strand5pOverhangSeqLengthLim(1) || ...
        overhang(1) > model.strand5pOverhangSeqLengthLim(2)
    if Param.Verbose
        fprintf(' 5'' strand 3'' overhang %d out of range [%d %d].', ...
            overhang(1), model.strand5pOverhangSeqLengthLim(1), ...
            model.strand5pOverhangSeqLengthLim(2));
    end
    return;
end

% calculate 3' strand 3' overhang
overhang(2) = org.mensxmachina.mirna.duplex3pstrandoverhangseqlengthq...
    (hairpinBracket, hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex);

if overhang(2) < model.strand3pOverhangSeqLengthLim(1) || ...
        overhang(2) > model.strand3pOverhangSeqLengthLim(2)
    if Param.Verbose
        fprintf(' 3'' strand 3'' overhang %d out of range [%d %d].', ...
            overhang(2), model.strand3pOverhangSeqLengthLim(1), ...
            model.strand3pOverhangSeqLengthLim(2));
    end 
    return;
end

tf = true;

end