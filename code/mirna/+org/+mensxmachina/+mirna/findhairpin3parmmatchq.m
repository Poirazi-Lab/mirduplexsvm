function pos3p = findhairpin3parmmatchq(match5pArmPos, match3pArmPos, pos5p)    
%FINDHAIRPIN3PARMMATCHQ
% %% run it alone
% clear all;
% clc;
% load('findhairpin3parmmatchqINPUT');
%%
totalHairpinBasePairs = length(match5pArmPos);
pos3p = match3pArmPos( totalHairpinBasePairs - find(match5pArmPos == pos5p) + 1 );

end