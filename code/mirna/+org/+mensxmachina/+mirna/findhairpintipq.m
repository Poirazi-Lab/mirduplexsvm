function pos = findhairpintipq(hairpin5pArmMatchPos, hairpin3pArmMatchPos)
%FINDHAIRPINTIP Find hairpin tip
    
hairpinLoop5pEndPos = hairpin5pArmMatchPos(end) + 1; % loop start = just after last 5' match
hairpinLoop3pEndPos = hairpin3pArmMatchPos(1) - 1; % loop end = just before first 3' match
hairpinLoopSequenceLength = hairpinLoop3pEndPos - hairpinLoop5pEndPos + 1; % length of the loop

% determine tip position
pos = hairpinLoop5pEndPos - 1 + ceil(hairpinLoopSequenceLength/2);

end