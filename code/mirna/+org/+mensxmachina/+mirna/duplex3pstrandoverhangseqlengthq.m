function overhang = duplex3pstrandoverhangseqlengthq(hairpinBracket, hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex)
%DUPLEX3PSTRANDOVERHANGSEQLENGTHQ  miRNA:miRNA* duplex 3' strand overhang sequence length

%     V-- 1st upstream matching position 
%   V-- 5' strand 5' end
%   ..((((
% ....))))
% ^-- 3' strand 3' end
%   ^-- opposite position from 5' strand 5' end end
% || 2nt overhang

duplex5pStrand5pEndPos = duplex(1);
duplex3pStrand3pEndPos = duplex(4);

% find offset of the first upstream matching position from duplex 5' strand 5' end

offset = 0; % initialize offset to zero

match = hairpinBracket(duplex5pStrand5pEndPos + offset) ~= '.'; % determine whether there is a match at that offset

if match
    
    oppositePos = org.mensxmachina.mirna.findhairpin3parmmatchq(hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex5pStrand5pEndPos + offset);
    match = oppositePos <= duplex3pStrand3pEndPos;
    
end

while ~match % while there is no match at the current offset
    
    offset = offset + 1; % increment offset
    
    match = hairpinBracket(duplex5pStrand5pEndPos + offset) ~= '.'; % determine whether there is a match at that offset
    
    if match
        
        oppositePos = org.mensxmachina.mirna.findhairpin3parmmatchq(hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex5pStrand5pEndPos + offset);
        match = oppositePos <= duplex3pStrand3pEndPos;
        
    end
    
end

% find opposite position from 5' strand 5' end
duplex5pStrand5pEndOppositePos = oppositePos + offset;

% calculate overhang
overhang = duplex3pStrand3pEndPos - duplex5pStrand5pEndOppositePos;

end