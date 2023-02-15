 function errorPlotter(miRnaDuplexFinderTrainConfig, miRnaDuplexFinderTestConfig, ...
     miRnaDuplexEst, expName, hold_out)

load('param');

% configuration
meanAbsErrorEdges = 0:0.25:5; % mean absolute error values in plots
meanAbsErrorEdgesDroshaDicer = 0:0.5:5;
meanAbsErrorEdges5pstarts = 0:0.5:5;
numFinders = size(miRnaDuplexFinder, 1);
numMeanAbsErrorEdges = length(meanAbsErrorEdges);
xAxisMeanError = 0:20;
numMeanAbsErrorEdgesDroshaDicer = length(meanAbsErrorEdgesDroshaDicer);
numMeanAbsErrorEdges5pstarts = length(meanAbsErrorEdges5pstarts);
xAxisDroshaDicer = 0:10;

% initialization
meanAbsError = zeros(numFinders, 1); % mean absolute error per finder
meanAbsError5p5p = zeros(numFinders, 1);
meanAbsError5p3p = zeros(numFinders, 1);
meanAbsError3p5p = zeros(numFinders, 1);
meanAbsError3p3p = zeros(numFinders, 1);
meanAbsErrorDrosha = zeros(numFinders, 1);
meanAbsErrorDicer = zeros(numFinders, 1);
meanAbsError5pstarts = zeros(numFinders, 1);


absErrorStd = zeros(numFinders, 1); % absolute error std per finder
absErrorStd5p5p = zeros(numFinders, 1);
absErrorStd5p3p = zeros(numFinders, 1);
absErrorStd3p5p = zeros(numFinders, 1);
absErrorStd3p3p = zeros(numFinders, 1);
absErrorStdDrosha = zeros(numFinders, 1);
absErrorStdDicer = zeros(numFinders, 1);
absErrorStd5pstarts = zeros(numFinders, 1);

% mean absolute error cumulative relative frequencies per finder
meanAbsErrorMeanCumRelFreq = zeros(numMeanAbsErrorEdges, numFinders); 
meanAbsErrorMeanCumRelFreqDrosha = zeros(numMeanAbsErrorEdgesDroshaDicer, numFinders); 
meanAbsErrorMeanCumRelFreqDicer = zeros(numMeanAbsErrorEdgesDroshaDicer, numFinders); 
meanAbsErrorMeanCumRelFreq5pstarts = zeros(numMeanAbsErrorEdges5pstarts, numFinders); 

for i = 1%:numFinders % for each finder

    iMeanAbsError = zeros(crossValK, 1);
    i5p5pMeanAbsError = zeros(crossValK, 1);
    i5p3pMeanAbsError = zeros(crossValK, 1);
    i3p5pMeanAbsError = zeros(crossValK, 1);
    i3p3pMeanAbsError = zeros(crossValK, 1);
    iMeanAbsErrorDrosha = zeros(crossValK, 1);
    iMeanAbsErrorDicer = zeros(crossValK, 1);
    iMeanAbsError5pstarts = zeros(crossValK, 1);
    
    iAbsErrorStd = zeros(crossValK, 1);
    i5p5pAbsErrorStd = zeros(crossValK, 1);
    i5p3pAbsErrorStd = zeros(crossValK, 1);
    i3p5pAbsErrorStd = zeros(crossValK, 1);
    i3p3pAbsErrorStd = zeros(crossValK, 1);
    iAbsErrorStdDrosha = zeros(crossValK, 1);
    iAbsErrorStdDicer = zeros(crossValK, 1);
    iAbsErrorStd5pstarts = zeros(crossValK, 1);
    
    
    iMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdges, crossValK);
    iDroshaMeanCutErrorCumRelFreq = zeros(numMeanAbsErrorEdgesDroshaDicer, crossValK);
    iDicerMeanCutErrorCumRelFreq = zeros(numMeanAbsErrorEdgesDroshaDicer, crossValK);
    i5pstartsCutErrorCumRelFreq = zeros(numMeanAbsErrorEdges5pstarts, crossValK);
    
    
    for j = 1:crossValK % for each cross validation fold
        
        if(hold_out)
            jSampleName = sprintf('test');
        else
            jSampleName = sprintf('test-%d', j);
        end
        
        iFinderName = miRnaDuplexFinder.Properties.ObsNames{i};
        
        % find configuration with i-th finder and j-th test set
        ijConfigInd = find(strcmp(miRnaDuplexFinderTrainConfig.finderName, ...
            iFinderName) & strcmp(miRnaDuplexFinderTestConfig.sampleName, ...
            jSampleName));
        
        % load sample
        load(['output' filesep 'data' filesep 'hairpin_' jSampleName], ...
            'hairpin')
        % get number of hairpinLongSequencess
        jNumhairpinLongSequencess = size(hairpin, 1);
        
        % calculate errors
        ij5pStrand5pError = miRnaDuplexEst.miRnaDuplex5pStrand5pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex5pStrand5pEndPos;
        ij5pStrand3pError = miRnaDuplexEst.miRnaDuplex5pStrand3pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex5pStrand3pEndPos;
        ij3pStrand5pError = miRnaDuplexEst.miRnaDuplex3pStrand5pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex3pStrand5pEndPos;
        ij3pStrand3pError = miRnaDuplexEst.miRnaDuplex3pStrand3pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex3pStrand3pEndPos;

        lengthIN = length(ij5pStrand5pError);
        
        if j == 1 
            here = 'here';
            p1 = 0;
        else
            p1 = p1 + lengthOUT;
            
        end
        
        posStart1 = p1 + 1;
        
        posEnd1 = posStart1 + lengthIN - 1;
               
        allErrors5p5p(posStart1:posEnd1,i) = ij5pStrand5pError;
        allErrors5p3p(posStart1:posEnd1,i) = ij5pStrand3pError;
        allErrors3p5p(posStart1:posEnd1,i) = ij3pStrand5pError;
        allErrors3p3p(posStart1:posEnd1,i) = ij3pStrand3pError;
        
        lengthOUT = lengthIN;
        
        % errors equals to error
        ij5pStrand5pAbsError = abs(ij5pStrand5pError);
        ij5pStrand3pAbsError = abs(ij5pStrand3pError);
        ij3pStrand5pAbsError = abs(ij3pStrand5pError);
        ij3pStrand3pAbsError = abs(ij3pStrand3pError);
        
        % calculate mean (along ends) error
        ijMeanAbsError = mean([ij5pStrand5pAbsError ij5pStrand3pAbsError...
            ij3pStrand5pAbsError ij3pStrand3pAbsError], 2);        
        ijDroshaCutError = mean([ij5pStrand5pAbsError ...
            ij3pStrand3pAbsError], 2);        
        ijDicerCutError = mean([ij5pStrand3pAbsError...
            ij3pStrand5pAbsError], 2);        
        ij5pStartsError = mean([ij5pStrand5pAbsError ...
            ij3pStrand5pAbsError], 2);
        
        iMeanAbsErrorCumRelFreq(:, j) = histc(ijMeanAbsError, ...
            meanAbsErrorEdges)/jNumhairpinLongSequencess';        
        iDroshaMeanCutErrorCumRelFreq(:, j) = histc(ijDroshaCutError, ...
            meanAbsErrorEdgesDroshaDicer)/jNumhairpinLongSequencess';        
        iDicerMeanCutErrorCumRelFreq(:, j) = histc(ijDicerCutError, ...
            meanAbsErrorEdgesDroshaDicer)/jNumhairpinLongSequencess';
        i5pstartsCutErrorCumRelFreq(:, j) = histc(ij5pStartsError, ...
            meanAbsErrorEdges5pstarts)/jNumhairpinLongSequencess';
        
        % calculate mean (along ends and hairpinLongSequencess) error
        iMeanAbsError(j) = mean(ijMeanAbsError);
        i5p5pMeanAbsError(j) = mean(ij5pStrand5pAbsError);
        i5p3pMeanAbsError(j) = mean(ij5pStrand3pAbsError);
        i3p5pMeanAbsError(j) = mean(ij3pStrand5pAbsError);
        i3p3pMeanAbsError(j) = mean(ij3pStrand3pAbsError);
        iMeanAbsErrorDrosha(j) = mean(ijDroshaCutError);
        iMeanAbsErrorDicer(j) = mean(ijDicerCutError);
        iMeanAbsError5pstarts(j) = mean(ij5pStartsError);
        
        % calculate std of mean (along ends) error
        i5p5pAbsErrorStd(j) = std(ij5pStrand5pAbsError);
        i5p3pAbsErrorStd(j) = std(ij5pStrand3pAbsError);
        i3p5pAbsErrorStd(j) = std(ij3pStrand5pAbsError);
        i3p3pAbsErrorStd(j) = std(ij3pStrand3pAbsError);
        iAbsErrorStd(j) = std(ijMeanAbsError);
        iAbsErrorStdDrosha(j) = std(ijDroshaCutError);
        iAbsErrorStdDicer(j) = std(ijDicerCutError);
        iAbsErrorStd5pstarts(j) = std(ij5pStartsError);
    end
    
    % calculate mean (along ends, hairpinLongSequencess and cross validation folds) error
    meanAbsError(i) = mean(iMeanAbsError);
    meanAbsError5p5p(i) = mean(i5p5pMeanAbsError);
    meanAbsError5p3p(i) = mean(i5p3pMeanAbsError);
    meanAbsError3p5p(i) = mean(i3p5pMeanAbsError);
    meanAbsError3p3p(i) = mean(i3p3pMeanAbsError);
    meanAbsErrorDrosha(i) = mean(iAbsErrorStdDrosha);
    meanAbsErrorDicer(i) = mean(iAbsErrorStdDicer);
    meanAbsError5pstarts(i) = mean(iAbsErrorStd5pstarts);
    
    % calculate mean (along cross validation folds) std of mean (along ends) error
    absErrorStd(i) = mean(iAbsErrorStd);
    absErrorStd5p5p(i) = mean(i5p5pAbsErrorStd);
    absErrorStd5p3p(i) = mean(i5p3pAbsErrorStd);
    absErrorStd3p5p(i) = mean(i3p5pAbsErrorStd);
    absErrorStd3p3p(i) = mean(i3p3pAbsErrorStd);
    absErrorStdDrosha(i) = mean(meanAbsErrorDrosha);
    absErrorStdDicer(i) = mean(meanAbsErrorDicer);
    absErrorStd5pstarts(i) = mean(meanAbsError5pstarts);

    % calculate mean (along ends and hairpinLongSequencess) error mean (along 
    % cross validation folds) cumulative relative frequencies
    meanAbsErrorMeanCumRelFreq(:, i) = mean(iMeanAbsErrorCumRelFreq, 2);
    meanAbsErrorMeanCumRelFreqDrosha(:, i) = mean(iDroshaMeanCutErrorCumRelFreq, 2);
    meanAbsErrorMeanCumRelFreqDicer(:, i) = mean(iDicerMeanCutErrorCumRelFreq, 2);
    meanAbsErrorMeanCumRelFreq5pstarts(:, i) = mean(i5pstartsCutErrorCumRelFreq, 2);

end

flankSequece = int2str...
    (miRnaDuplexFinderTrainConfig.trainParam{1,1}.FlankingSequenceLength);

% Calculate errors of start positions all finders
total = length(allErrors5p5p);
allErrorsFromStart = cell(1,numFinders);
wrong = -5:1:5;

for i = 1%:numFinders
    distanceFromStart = zeros(4,length(wrong));
    start5p5p = histc(allErrors5p5p(:,i), wrong)/ total;
    start5p3p = histc(allErrors5p3p(:,i), wrong)/ total;
    start3p5p = histc(allErrors3p5p(:,i), wrong)/ total;
    start3p3p = histc(allErrors3p3p(:,i), wrong)/ total;
    
    distanceFromStart(1,:) = start5p5p';
    distanceFromStart(2,:) = start5p3p';
    distanceFromStart(3,:) = start3p5p';
    distanceFromStart(4,:) = start3p3p';               
    
    allErrorsFromStart{1,i} = distanceFromStart;
end

Y5p5p = zeros(numFinders,length(distanceFromStart));
Y5p3p = zeros(numFinders,length(distanceFromStart));
Y3p5p = zeros(numFinders,length(distanceFromStart));
Y3p3p = zeros(numFinders,length(distanceFromStart));
for i = 1%:numFinders
    Y5p5p(i,:) = allErrorsFromStart{1,i}(1,:);
    Y5p3p(i,:) = allErrorsFromStart{1,i}(2,:);
    Y3p5p(i,:) = allErrorsFromStart{1,i}(3,:);
    Y3p3p(i,:) = allErrorsFromStart{1,i}(4,:);
end

xaxisLength = wrong;
figure
% subplot(2,2,1)
bar(xaxisLength , Y5p5p');
xlabel('Distance from truth nts','fontsize',12)
ylabel('Percentage of hairpins','fontsize',12);
legend(miRnaDuplexFinder.label, 'Location', 'NorthWest','fontsize',12);
title('Distribution of duplex 5p strand 5p end error','fontsize',12);

figure
% subplot(2,2,2)
bar(xaxisLength , Y5p3p');
xlabel('Distance from truth nts','fontsize',12)
ylabel('Percentage of hairpins','fontsize',12);
legend(miRnaDuplexFinder.label, 'Location', 'NorthWest','fontsize',12);
title('Distribution of duplex 5p strand 3p end error','fontsize',12);

figure
% subplot(2,2,3)
bar(xaxisLength , Y3p3p');
xlabel('Distance from truth nts','fontsize',12)
ylabel('Percentage of hairpins','fontsize',12);
legend(miRnaDuplexFinder.label, 'Location', 'NorthWest','fontsize',12);
title('Distribution of duplex 3p strand 5p end error','fontsize',12);

figure
% subplot(2,2,4)
bar(xaxisLength , Y3p5p');
xlabel('Distance from truth nts','fontsize',12)
ylabel('Percentage of hairpins','fontsize',12);
legend(miRnaDuplexFinder.label, 'Location', 'NorthWest','fontsize',12);
title('Distribution of duplex 3p strand 3p end error','fontsize',12);

% Calculate errors of length all finders
wrong = -5:1:5;
allErrorsLength = cell(1,numFinders);

for i = 1%:numFinders
    lengthMistakes = zeros(2,length(wrong));
    rightStart5p5p = find(allErrors5p5p(:,i) == 0);
    rightStart3p5p = find(allErrors3p5p(:,i) == 0);
    
    % 5' strand
    starts5p3p = allErrors5p3p(rightStart5p5p);
    perc5p = histc(starts5p3p, wrong)/length(rightStart5p5p);
    lengthMistakes(1,:) = perc5p;
    
    % 3' strand
    starts3p3p = allErrors3p3p(rightStart3p5p);
    perc3p = histc(starts3p3p, wrong)/length(rightStart3p5p);
    lengthMistakes(2,:) = perc3p;
    
    allErrorsLength{1,i} = lengthMistakes;
end

L5p3p = zeros(numFinders,length(lengthMistakes));
L3p3p = zeros(numFinders,length(lengthMistakes));
for i = 1%:numFinders
    L5p3p(i,:) = allErrorsLength{1,i}(1,:);
    L3p3p(i,:) = allErrorsLength{1,i}(2,:);
end

xaxisLength = wrong;
figure
subplot(2,1,1)
bar(xaxisLength , L5p3p');
xlabel('Distance from truth nts')
ylabel('Relative frequency');
title('Errors in length for the predicted matures 5strand');
legend(miRnaDuplexFinder.label);%, 'Location', 'NorthWest');
subplot(2,1,2)
bar(xaxisLength , L3p3p');
xlabel('Distance from truth nts')
ylabel('Relative frequency');
title('Errors in length for the predicted matures 3strand');
legend(miRnaDuplexFinder.label);%, 'Location', 'NorthWest');

% Calculate CORRELATIONs all finders
correlationTable = dataset(...
    {zeros(2, numFinders), 'Strand5_Corner5p_Coorelation_Strand5_Corner3p'}, ...
    {zeros(2, numFinders), 'Strand5_Corner3p_Coorelation_Strand3_Corner5p'}, ...
    {zeros(2, numFinders), 'Strand3_Corner5p_Coorelation_Strand3_Corner3p'}, ...
    {zeros(2, numFinders), 'Strand5_Corner5p_Coorelation_Strand3_Corner3p'}, ...
    {zeros(2, numFinders), 'Strand5_Corner5p_Coorelation_Strand3_Corner5p'}, ...
    {zeros(2, numFinders), 'Strand5_Corner3p_Coorelation_Strand3_Corner3p'}, ...
    'ObsNames', arrayfun(@(i) sprintf('temp%d', i), 1:2, ...
    'UniformOutput', false) ...
    );

 for i = 1%: numFinders
    correlationTable.Properties.ObsNames{1} = 'Coorelation';
    correlationTable.Properties.ObsNames{2} = 'pValue';
    
    % 5 strand corners
    [astrand5pC5p3p, pvalStrand5C5p3p] = corr(allErrors5p5p(:,i), allErrors5p3p(:,i));
    correlationTable.Strand5_Corner5p_Coorelation_Strand5_Corner3p(1,i) = ...
        astrand5pC5p3p;
    correlationTable.Strand5_Corner5p_Coorelation_Strand5_Corner3p(2,i) = ...
        pvalStrand5C5p3p;
   
    % Dicer cut corners          
    [astrand53C3p5p, pvalstrand53C3p5p] = corr(allErrors5p3p(:,i), allErrors3p5p(:,i));
    correlationTable.Strand5_Corner3p_Coorelation_Strand3_Corner5p(1,i) = ...
        astrand53C3p5p;
    correlationTable.Strand5_Corner3p_Coorelation_Strand3_Corner5p(2,i) = ...
        pvalstrand53C3p5p;
    
    % 3 strand corners
    [astrand3C5p3p, pvalstrand3C5p3p] = corr(allErrors3p5p(:,i), allErrors3p3p(:,i));
    correlationTable.Strand3_Corner5p_Coorelation_Strand3_Corner3p(1,i) = astrand3C5p3p;
    correlationTable.Strand3_Corner5p_Coorelation_Strand3_Corner3p(2,i) = pvalstrand3C5p3p;            
    
    % Drosha cut corners
    [astrand53C5p3p, pvalstrand53C5p3p] = corr(allErrors5p5p(:,i), allErrors3p3p(:,i));
    correlationTable.Strand5_Corner5p_Coorelation_Strand3_Corner3p(1,i) = astrand53C5p3p;
    correlationTable.Strand5_Corner5p_Coorelation_Strand3_Corner3p(2,i) = pvalstrand53C5p3p;
    
    % 5p starts corners
    [astrand53C5p5p, pvalstrand535p5p] = corr(allErrors5p5p(:,i), allErrors3p5p(:,i));
            
    correlationTable.Strand5_Corner5p_Coorelation_Strand3_Corner5p(1,i) = astrand53C5p5p;
    correlationTable.Strand5_Corner5p_Coorelation_Strand3_Corner5p(2,i) = pvalstrand535p5p;
    
    % 3p start corners
    [astrand53C3p3p, pvalstrand533p3p] = corr(allErrors5p3p(:,i), allErrors3p3p(:,i));
    correlationTable.Strand5_Corner3p_Coorelation_Strand3_Corner3p(1,i) = astrand53C3p3p;
    correlationTable.Strand5_Corner3p_Coorelation_Strand3_Corner3p(2,i) = pvalstrand533p3p;

end

% ----- plot relative frequency of absolute mean error --------------
figure;
bar(xAxisMeanError, meanAbsErrorMeanCumRelFreq);
xlabel('Exact Nucleotides Error for 4 corners');
ylabel('Mean cumulative relative frequency');
legend(miRnaDuplexFinder.label);%, 'Location', 'NorthWest');
title('Sequence')
    
figure;% = figure(iFigureProperties{:});
subplot(2,2,1)
bar(xAxisDroshaDicer, meanAbsErrorMeanCumRelFreqDrosha);
xlabel('Exact Nucleotides Error for DROSHA cut');
ylabel('Mean cumulative relative frequency');
legend(miRnaDuplexFinder.label);%, 'Location', 'NorthWest');
title('Sequence')
subplot(2,2,2)
bar(xAxisDroshaDicer, meanAbsErrorMeanCumRelFreqDicer);
xlabel('Exact Nucleotides Error for DICER cut');
ylabel('Mean cumulative relative frequency');
legend(miRnaDuplexFinder.label);%, 'Location', 'NorthWest');
title('Sequence')
subplot(2,2,3)
bar(xAxisDroshaDicer, meanAbsErrorMeanCumRelFreq5pstarts);
xlabel('Exact Nucleotides Error for 5pstarts cut');
ylabel('Mean cumulative relative frequency');
legend(miRnaDuplexFinder.label);%, 'Location', 'NorthWest');
title('Sequence')

if (hold_out)
    save(['Results' filesep expName, '_hold_out_',flankSequece], 'meanAbsErrorMeanCumRelFreq',...
        'meanAbsErrorMeanCumRelFreqDrosha','meanAbsErrorMeanCumRelFreqDicer',...
        'meanAbsErrorMeanCumRelFreq5pstarts', ...
        'correlationTable', 'allErrorsFromStart' , ...
        'allErrorsLength');
else
    save(['Results' filesep expName, '_CrossVal_',flankSequece], 'meanAbsErrorMeanCumRelFreq',...
        'meanAbsErrorMeanCumRelFreqDrosha','meanAbsErrorMeanCumRelFreqDicer',...
        'meanAbsErrorMeanCumRelFreq5pstarts', ...
        'correlationTable', 'allErrorsFromStart' , ...
        'allErrorsLength');
end


 end
  
