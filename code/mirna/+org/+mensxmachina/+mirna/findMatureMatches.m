function matches = findMatureMatches(hairpinBracket, match5pPos, match3pPos, duplex)
%ISDUPLEXQ Summary of this function goes here
%   Detailed explanation goes here
% %% to run it alone
% clear all;
% clc;
% load ('isduplexqINPUT');
%%

duplex5pStrand5pEndPos = duplex(1);
duplex5pStrand3pEndPos = duplex(2);
duplex3pStrand5pEndPos = duplex(3);
duplex3pStrand3pEndPos = duplex(4);
%%
matches = 0;
for i = duplex5pStrand5pEndPos : duplex5pStrand3pEndPos
    
    if hairpinBracket(i) == '(' % match
        
        pos3p = org.mensxmachina.mirna.findhairpin3parmmatchq(match5pPos, ...
            match3pPos, i);
        
        if pos3p >= duplex3pStrand5pEndPos && pos3p <= duplex3pStrand3pEndPos
            
            % OK, there is a match with the mature on the 3' arm
            matches = matches + 1;
%             if matches >= 10
%                 %tf = 1;
%                 return;
%             end
        end
        
    end

end
%matches
%tf = 0;

end