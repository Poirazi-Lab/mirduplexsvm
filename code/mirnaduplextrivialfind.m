function miRnaDuplexEst = mirnaduplextrivialfind(model, seq, bracket, thermo, a)
%MIRNADUPLEXTRIVIALFIND Find miRNA:miRNA*-duplex using trivial finder

miRnaDuplexEst = org.mensxmachina.mirna.mirnaduplextrivialfindq(model, seq, bracket);

end