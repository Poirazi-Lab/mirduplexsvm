function model = mirnaduplextrivialfindertrain(seq, bracket, thermo ,miRnaDuplex, a)
%MIRNADUPLEXTRIVIALFINDERTRAIN Train trivial miRNA:miRNA*-duplex finder

model = org.mensxmachina.mirna.mirnaduplextrivialfindertrainq(seq, bracket, miRnaDuplex);

end