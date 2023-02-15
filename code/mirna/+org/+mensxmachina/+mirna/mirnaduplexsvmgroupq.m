function group = mirnaduplexsvmgroupq(tf)
%MIRNADUPLEXSVMGROUPQ Create miRNA:miRNA duplex SVM grouping variable

group = - ones(size(tf));
group(tf) = 1;

end

