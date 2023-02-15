function [group score] = mirnaduplexsvmclassifyq(model, x)
%MIRNADUPLEXSVMCLASSIFYQ Classify candidate miRNA:miRNA* duplex using miRNA:miRNA* duplex SVM

[group, a, score] = svmpredict(zeros(size(x, 1), 1), x, model);

end