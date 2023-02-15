function model = mirnaduplexsvmtrainq(x, group, Param)
%MIRNADUPLEXSVMTRAINQ Train miRNA:miRNA* duplex SVM

% calculate weights

numTrueInstances = sum(group == 1);
numFalseInstances = sum(group == -1);

numInstances = numTrueInstances + numFalseInstances;
Param.cost;
w1 = numInstances/(Param.cost*numTrueInstances);
wm1 = numInstances/(Param.cost*numFalseInstances);

if strcmp(Param.KernelType, 'radbas')
    t = 2;
else % poly
    t = 1;
end

options = sprintf('-t %d -w1 %f -w-1 %f', t, w1, wm1);

if t == 1 % poly
    % add degree in options
    options = [options sprintf(' -d %d', Param.Degree)];
end

% train SVM
options
model = svmtrain(group, x, options);

end

