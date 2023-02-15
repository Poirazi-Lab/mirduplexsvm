function overhang = duplex5pstrandoverhangseqlengthq(hairpinBracket, hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex)
%DUPLEX5PSTRANDOVERHANGSEQLENGTHQ miRNA:miRNA* duplex 5' strand overhang sequence length

%        || 2nt overhang
%       V-- opposite position from 3' strand 5' end
%         V- 5' strand 3' end
% (((((....
% )))))..
%       ^-- 3' strand 5' end
%     ^-- 1st upstream matching position

duplex5pStrand3pEndPos = duplex(2);
duplex3pStrand5pEndPos = duplex(3);

% find offset of the first upstream matching position from duplex 3' strand5' end

offset = 0; % initialize offset to zero

match = hairpinBracket(duplex3pStrand5pEndPos + offset) ~= '.'; % determine whether there is a match at that offset

if match
    
    oppositePos = org.mensxmachina.mirna.findhairpin5parmmatchq(hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex3pStrand5pEndPos + offset);
    match = oppositePos <= duplex5pStrand3pEndPos;
    
end

while ~match % while there is no match at the current offset
    
    offset = offset + 1; % increment offset
    
    match = hairpinBracket(duplex3pStrand5pEndPos + offset) ~= '.'; % determine whether there is a match at that offset
    
    if match
    
        oppositePos = org.mensxmachina.mirna.findhairpin5parmmatchq(hairpin5pArmMatchPos, hairpin3pArmMatchPos, duplex3pStrand5pEndPos + offset);
        match = oppositePos <= duplex5pStrand3pEndPos;
    
    end
    
end

% find opposite position from 3' strand 5' end
duplex3pStrand5pEndOppositePos = oppositePos + offset;

% calculate overhang
overhang = duplex5pStrand3pEndPos - duplex3pStrand5pEndOppositePos;

end