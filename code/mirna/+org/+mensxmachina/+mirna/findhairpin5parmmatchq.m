function pos5p = findhairpin5parmmatchq(match5pArmPos, match3pArmPos, pos3p)
%FINDHAIRPIN5PARMMATCHQ

totalHairpinBasePairs = length(match5pArmPos);    
pos5p = match5pArmPos( totalHairpinBasePairs - find(match3pArmPos == pos3p) + 1 );

end