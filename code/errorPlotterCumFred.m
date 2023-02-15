 function [meanAbsError absErrorStd] = errorPlotterCumFred(miRnaDuplexFinderTrainConfig, miRnaDuplexFinderTestConfig, ...
     miRnaDuplexEst, expName, hold_out)

load('param');

% configuration
meanAbsErrorEdges = 0:0.25:5; % mean absolute error values in plots
meanAbsErrorEdgesCorner = 0:5;
meanAbsErrorEdgesDroshaDicer = 0:0.5:5;
meanAbsErrorEdges5pstarts = 0:0.5:5;

numFinders = size(miRnaDuplexFinder, 1);
% numMedia = size(medium, 1);
numMeanAbsErrorEdges = length(meanAbsErrorEdges);
xAxisMeanError = 0:numMeanAbsErrorEdges-1;
numMeanAbsErrorEdgesCorner = length(meanAbsErrorEdgesCorner);
numMeanAbsErrorEdgesDroshaDicer = length(meanAbsErrorEdgesDroshaDicer);
numMeanAbsErrorEdges5pstarts = length(meanAbsErrorEdges5pstarts);

xAxisDroshaDicer = 0:numMeanAbsErrorEdgesDroshaDicer-1;

% initialization
meanAbsError = zeros(numFinders, 1); % mean absolute error per finder
meanAbsError5p5p = zeros(numFinders, 1);
meanAbsError5p3p = zeros(numFinders, 1);
meanAbsError3p5p = zeros(numFinders, 1);
meanAbsError3p3p = zeros(numFinders, 1);
meanAbsErrorDrosha = zeros(numFinders, 1);
meanAbsErrorDicer = zeros(numFinders, 1);
meanAbsError5pstarts = zeros(numFinders, 1);

% absolute error std per finder
absErrorStd = zeros(numFinders, 1); 
absErrorStd5p5p = zeros(numFinders, 1);
absErrorStd5p3p = zeros(numFinders, 1);
absErrorStd3p5p = zeros(numFinders, 1);
absErrorStd3p3p = zeros(numFinders, 1);
absErrorStdDrosha = zeros(numFinders, 1);
absErrorStdDicer = zeros(numFinders, 1);
absErrorStd5pstarts = zeros(numFinders, 1);

% mean absolute error cumulative relative frequencies per finder
meanAbsErrorMeanCumRelFreq = zeros(numMeanAbsErrorEdges, numFinders); 
f5p5pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, numFinders);
f5p3pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, numFinders);
f3p5pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, numFinders);
f3p3pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, numFinders);
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
    i5p5pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, crossValK);
    i5p3pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, crossValK);
    i3p5pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, crossValK);
    i3p3pMeanAbsErrorCumRelFreq = zeros(numMeanAbsErrorEdgesCorner, crossValK);     
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
            iFinderName) & strcmp(miRnaDuplexFinderTestConfig.sampleName, jSampleName));
       
        %load sample
        load(['output' filesep 'data' filesep 'hairpin_' jSampleName], ...
            'hairpin')

        % get number of hairpins
        jNumhairpins = size(hairpin, 1);
        
        % calculate errors
        ij5pStrand5pError = miRnaDuplexEst.miRnaDuplex5pStrand5pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex5pStrand5pEndPos;
        ij5pStrand3pError = miRnaDuplexEst.miRnaDuplex5pStrand3pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex5pStrand3pEndPos;
        ij3pStrand5pError = miRnaDuplexEst.miRnaDuplex3pStrand5pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex3pStrand5pEndPos;
        ij3pStrand3pError = miRnaDuplexEst.miRnaDuplex3pStrand3pEndPosEst{ijConfigInd} ...
            - hairpin.miRnaDuplex3pStrand3pEndPos;
        
        % calculate absolute errors
        ij5pStrand5pAbsError = abs(ij5pStrand5pError);
        ij5pStrand3pAbsError = abs(ij5pStrand3pError);
        ij3pStrand5pAbsError = abs(ij3pStrand5pError);
        ij3pStrand3pAbsError = abs(ij3pStrand3pError);
        
        % calculate mean (along ends) absolute error
        ijMeanAbsError = mean([ij5pStrand5pAbsError ij5pStrand3pAbsError...
            ij3pStrand5pAbsError ij3pStrand3pAbsError], 2);
        
        ijDroshaCutError = mean([ij5pStrand5pAbsError ...
            ij3pStrand3pAbsError], 2);
        
        ijDicerCutError = mean([ij5pStrand3pAbsError...
            ij3pStrand5pAbsError], 2);
        
        ij5pStartsError = mean([ij5pStrand5pAbsError ...
            ij3pStrand5pAbsError], 2);

        % calculate its cumulative relative distribution
        iMeanAbsErrorCumRelFreq(:, j) = cumsum(histc(ijMeanAbsError, ...
            meanAbsErrorEdges)/jNumhairpins)';        
        i5p5pMeanAbsErrorCumRelFreq(:, j) = cumsum(histc(ij5pStrand5pAbsError, ...
            meanAbsErrorEdgesCorner)/jNumhairpins)';
        i5p3pMeanAbsErrorCumRelFreq(:, j) = cumsum(histc(ij5pStrand3pAbsError, ...
            meanAbsErrorEdgesCorner)/jNumhairpins)';
        i3p5pMeanAbsErrorCumRelFreq(:, j) = cumsum(histc(ij3pStrand5pAbsError, ...
            meanAbsErrorEdgesCorner)/jNumhairpins)';
        i3p3pMeanAbsErrorCumRelFreq(:, j) = cumsum(histc(ij3pStrand3pAbsError, ...
            meanAbsErrorEdgesCorner)/jNumhairpins)';
        iDroshaMeanCutErrorCumRelFreq(:, j) = cumsum(histc(ijDroshaCutError, ...
            meanAbsErrorEdgesDroshaDicer)/jNumhairpins)';        
        iDicerMeanCutErrorCumRelFreq(:, j) = cumsum(histc(ijDicerCutError, ...
            meanAbsErrorEdgesDroshaDicer)/jNumhairpins)';
        i5pstartsCutErrorCumRelFreq(:, j) = cumsum(histc(ij5pStartsError, ...
            meanAbsErrorEdges5pstarts)/jNumhairpins)';
        
        
        % calculate mean (along ends and hairpins) absolute error
        iMeanAbsError(j) = mean(ijMeanAbsError);
        i5p5pMeanAbsError(j) = mean(ij5pStrand5pAbsError);
        i5p3pMeanAbsError(j) = mean(ij5pStrand3pAbsError);
        i3p5pMeanAbsError(j) = mean(ij3pStrand5pAbsError);
        i3p3pMeanAbsError(j) = mean(ij3pStrand3pAbsError);
        iMeanAbsErrorDrosha(j) = mean(ijDroshaCutError);
        iMeanAbsErrorDicer(j) = mean(ijDicerCutError);
        iMeanAbsError5pstarts(j) = mean(ij5pStartsError);
        
        % calculate std of mean (along ends) absolute error
        i5p5pAbsErrorStd(j) = std(ij5pStrand5pAbsError);
        i5p3pAbsErrorStd(j) = std(ij5pStrand3pAbsError);
        i3p5pAbsErrorStd(j) = std(ij3pStrand5pAbsError);
        i3p3pAbsErrorStd(j) = std(ij3pStrand3pAbsError);
        iAbsErrorStd(j) = std(ijMeanAbsError);
        iAbsErrorStdDrosha(j) = std(ijDroshaCutError);
        iAbsErrorStdDicer(j) = std(ijDicerCutError);
        iAbsErrorStd5pstarts(j) = std(ij5pStartsError);
    end
    
    % calculate mean (along ends, hairpins and cross validation folds) absolute error
    meanAbsError(i) = mean(iMeanAbsError);
    meanAbsError5p5p(i) = mean(i5p5pMeanAbsError);
    meanAbsError5p3p(i) = mean(i5p3pMeanAbsError);
    meanAbsError3p5p(i) = mean(i3p5pMeanAbsError);
    meanAbsError3p3p(i) = mean(i3p3pMeanAbsError);
    meanAbsErrorDrosha(i) = mean(iAbsErrorStdDrosha);
    meanAbsErrorDicer(i) = mean(iAbsErrorStdDicer);
    meanAbsError5pstarts(i) = mean(iAbsErrorStd5pstarts);
    
    % calculate mean (along cross validation folds) std of mean (along ends) absolute error
    absErrorStd(i) = mean(iAbsErrorStd);
    absErrorStd5p5p(i) = mean(i5p5pAbsErrorStd);
    absErrorStd5p3p(i) = mean(i5p3pAbsErrorStd);
    absErrorStd3p5p(i) = mean(i3p5pAbsErrorStd);
    absErrorStd3p3p(i) = mean(i3p3pAbsErrorStd);
    absErrorStdDrosha(i) = mean(meanAbsErrorDrosha);
    absErrorStdDicer(i) = mean(meanAbsErrorDicer);
    absErrorStd5pstarts(i) = mean(meanAbsError5pstarts);

    % calculate mean (along ends and hairpins) absolute error mean (along 
    %cross validation folds) cumulative relative frequencies
    meanAbsErrorMeanCumRelFreq(:, i) = mean(iMeanAbsErrorCumRelFreq, 2);
    f5p5pMeanAbsErrorCumRelFreq(:, i) = mean(i5p5pMeanAbsErrorCumRelFreq, 2);
    f5p3pMeanAbsErrorCumRelFreq(:, i) = mean(i5p3pMeanAbsErrorCumRelFreq, 2);
    f3p5pMeanAbsErrorCumRelFreq(:, i) = mean(i3p5pMeanAbsErrorCumRelFreq, 2);
    f3p3pMeanAbsErrorCumRelFreq(:, i) = mean(i3p3pMeanAbsErrorCumRelFreq, 2);
    meanAbsErrorMeanCumRelFreqDrosha(:, i) = mean(iDroshaMeanCutErrorCumRelFreq, 2);
    meanAbsErrorMeanCumRelFreqDicer(:, i) = mean(iDicerMeanCutErrorCumRelFreq, 2);
    meanAbsErrorMeanCumRelFreq5pstarts(:, i) = mean(i5pstartsCutErrorCumRelFreq, 2);

end


 flankSequece = int2str...
    (miRnaDuplexFinderTrainConfig.trainParam{1,1}.FlankingSequenceLength);

    % -- plot cumulative relative frequency of mean absolute error --------    
    figure;
    bar(xAxisMeanError, meanAbsErrorMeanCumRelFreq);
    xlabel('MAE(nts)', 'fontsize',20);
    ylabel('Percentage of hairpins','fontsize',20);
    legend(miRnaDuplexFinder.label, 'Location', 'NorthWest','fontsize',20);    
    title('Cumulative distribution of mean absolute error', 'fontsize',20)
    print(sprintf('-f%d', gcf), '-r300', '-djpeg', ['Results' filesep ...
        expName, 'FlankSeq =',flankSequece, ' All Corner Errorsf5p5pMeanAbsErrorCumRelFreq']);
    
    figure;
    subplot(2,2,1)
    bar(xAxisDroshaDicer, meanAbsErrorMeanCumRelFreqDrosha);
    xlabel('Exact Nucleotides Error for DROSHA cut');
    ylabel('Percentage of hairpins');
    legend(miRnaDuplexFinder.label, 'Location', 'NorthWest');
    title('Cumulative distribution of mean absolute error')
    subplot(2,2,2)
    bar(xAxisDroshaDicer, meanAbsErrorMeanCumRelFreqDicer);
    xlabel('Exact Nucleotides Error for DICER cut');
    ylabel('Mean cumulative relative frequency');
    legend(miRnaDuplexFinder.label, 'Location', 'NorthWest');
    title('Sequence')
    subplot(2,2,3)
    bar(xAxisDroshaDicer, meanAbsErrorMeanCumRelFreq5pstarts);
    xlabel('Exact Nucleotides Error for 5pstarts cut');
    ylabel('Mean cumulative relative frequency');
    legend(miRnaDuplexFinder.label, 'Location', 'NorthWest');
    title('Sequence')
    print(sprintf('-f%d', gcf), '-r300', '-djpeg', ['Results' filesep ...
        expName, 'FlankSeq =',flankSequece, '5pstarts - Drosha-Dicer cuts Errors Cum-Freq']);

    figure;
    subplot(2,2,1)
    bar(meanAbsErrorEdgesCorner, f5p5pMeanAbsErrorCumRelFreq);
    xlabel('Exact Nucleotides Error: 5pstrand 5pend');
    ylabel('Mean cumulative relative frequency');
    legend(miRnaDuplexFinder.label, 'Location', 'SouthEast');
    title('Sequence')
    subplot(2,2,2)
    bar(meanAbsErrorEdgesCorner, f5p3pMeanAbsErrorCumRelFreq);
    xlabel('Exact Nucleotides Error: 5pstrand 3pend');
    ylabel('Mean cumulative relative frequency');
    legend(miRnaDuplexFinder.label, 'Location', 'SouthEast');
    title('Sequence')
    subplot(2,2,3)
    bar(meanAbsErrorEdgesCorner, f3p5pMeanAbsErrorCumRelFreq);
    xlabel('Exact Nucleotides Error: 3pstrand 5pend');
    ylabel('Mean cumulative relative frequency');
    legend(miRnaDuplexFinder.label, 'Location', 'SouthEast');
    title('Sequence')
    subplot(2,2,4)
    bar(meanAbsErrorEdgesCorner, f3p3pMeanAbsErrorCumRelFreq);
    xlabel('Exact Nucleotides Error: 3pstrand 3pend');
    ylabel('Mean cumulative relative frequency');
    legend(miRnaDuplexFinder.label, 'Location', 'SouthEast');
    title('Sequence')
    print(sprintf('-f%d', gcf), '-r300', '-djpeg', ['Results' filesep ...
        expName, 'FlankSeq =',flankSequece, ' Errors in each Corner Cum-Freq']);

    if(hold_out)
        save(['Results' filesep expName,'_hold_out_CumFreq_',flankSequece], 'meanAbsErrorMeanCumRelFreq',...
            'meanAbsErrorMeanCumRelFreqDrosha','meanAbsErrorMeanCumRelFreqDicer',...
            'meanAbsErrorMeanCumRelFreq5pstarts', ...
            'f5p5pMeanAbsErrorCumRelFreq','f5p3pMeanAbsErrorCumRelFreq', ...
            'f3p5pMeanAbsErrorCumRelFreq',...
            'f3p3pMeanAbsErrorCumRelFreq', 'iMeanAbsErrorCumRelFreq', ...
            'i5p5pMeanAbsErrorCumRelFreq', 'i3p5pMeanAbsErrorCumRelFreq',...
            'i5pstartsCutErrorCumRelFreq', 'iMeanAbsError');
 
    else
        save(['Results' filesep expName,'_CrossVal_CumFreq_',flankSequece], 'meanAbsErrorMeanCumRelFreq',...
            'meanAbsErrorMeanCumRelFreqDrosha','meanAbsErrorMeanCumRelFreqDicer',...
            'meanAbsErrorMeanCumRelFreq5pstarts', ...
            'f5p5pMeanAbsErrorCumRelFreq','f5p3pMeanAbsErrorCumRelFreq', ...
            'f3p5pMeanAbsErrorCumRelFreq',...
            'f3p3pMeanAbsErrorCumRelFreq', 'iMeanAbsErrorCumRelFreq', ...
            'i5p5pMeanAbsErrorCumRelFreq', 'i3p5pMeanAbsErrorCumRelFreq',...
            'i5pstartsCutErrorCumRelFreq', 'iMeanAbsError');
        
    end
    
    
    



 end
  
