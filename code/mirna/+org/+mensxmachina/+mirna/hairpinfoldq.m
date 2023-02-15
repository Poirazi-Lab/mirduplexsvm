function [bracket] = hairpinfoldq(seq)

[a, result] = system(sprintf('RNAfold -p -d0 --noLP --noPS <<< %s', seq));

% keep only ensemble prediction structure sequence
seqLength = length(seq);
bracket = result(seqLength + 2:2*seqLength + 1); 

end





