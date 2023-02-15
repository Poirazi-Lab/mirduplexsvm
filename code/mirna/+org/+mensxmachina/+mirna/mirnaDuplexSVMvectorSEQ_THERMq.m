function x = mirnaDuplexSVMvectorSEQ_THERMq(model, hairpinSeq, ...
    hairpinTipPos, candidateMiRnaDuplex)
%MIRNADUPLEXSVMVECTORQ Create miRNA:miRNA duplex SVM vector

strand = [5 3];

numSeqCols = 4*(4*model.flankSeqLength + model.strand5pSeqLengthLim(2) + ...
    model.strand3pSeqLengthLim(2));

% numThermoCols = numSeqCols/4;

numCols = numSeqCols;
x = zeros(1, numCols);

    
i = 1; % position counter se poio noykleotidio briskomai

for k = 1:2 % for each strand

    kStrandId = strand(k);

    % strand positions
    kStrand5pEndPos = candidateMiRnaDuplex((k - 1)*2 + 1);
    kStrand3pEndPos = candidateMiRnaDuplex((k - 1)*2 + 2);
    
    % strand sequence length 
    kStrandSeqLength = kStrand3pEndPos - kStrand5pEndPos + 1;

    % strand sequence length limits
    kStrandSeqLengthLim = model.(sprintf('strand%dpSeqLengthLim', kStrandId));
%     kStrandMinSeqLength = kStrandSeqLengthLim(1);
    kStrandMaxSeqLength = kStrandSeqLengthLim(2);
    
    % strand 5' and 3' sequence lengths
    kStrand5pEndSeqLength = ceil(kStrandSeqLength/2);
    kStrand3pEndSeqLength = floor(kStrandSeqLength/2); 
    
    % max strand 5' and 3' sequence lengths
    kStrandMax5pEndSeqLength = ceil(kStrandMaxSeqLength/2);
    kStrandMax3pEndSeqLength = floor(kStrandMaxSeqLength/2);


     % 5' end flanking sequence columns 
    for offset = model.flankSeqLength : -1 : 1

        pos = kStrand5pEndPos - offset;

        if (kStrandId == 5 && pos >= 1) || (kStrandId == 3 && pos >= hairpinTipPos + 1 )
            x((i - 1)*4 + double(nt2int(hairpinSeq(pos)))) = 1;
%             x(numSeqCols + i) = hairpinTherm(pos);
            
            
        end

        i = i + 1;

    end

    % strand 5' end columns
    for offset = 1 : kStrandMax5pEndSeqLength

        pos = kStrand5pEndPos + offset - 1;

        if offset <= kStrand5pEndSeqLength
            x((i - 1)*4 + double(nt2int(hairpinSeq(pos)))) = 1;
%             x(numSeqCols + i) = hairpinTherm(pos);
                  
        end

        i = i + 1;

    end

    % strand 3' end columns
    
    for offset = kStrandMax3pEndSeqLength : -1 : 1

        pos = kStrand3pEndPos - offset + 1;
        if offset <= kStrand3pEndSeqLength
            x((i - 1)*4 + double(nt2int(hairpinSeq(pos)))) = 1;
%             x(numSeqCols + i) = hairpinTherm(pos);
            %x(numSeqCols + i) = ~strcmp(hairpinBracket(pos), '.');
            %x(i) = hairpinTherm(pos);
            
        end

        i = i + 1;

    end

    % 3' end flanking sequence columns

    for offset = 1 : model.flankSeqLength

        pos = kStrand3pEndPos + offset;

        if (kStrandId == 5 && pos <= hairpinTipPos ) || (kStrandId == 3 && ...
                pos <= length(hairpinSeq) ) 
            x((i - 1)*4 + double(nt2int(hairpinSeq(pos)))) = 1;
%             x(numSeqCols + i) = hairpinTherm(pos);
            %x(numSeqCols + i) = ~strcmp(hairpinBracket(pos), '.');
            %x(i) = hairpinTherm(pos);
            
        end

        i = i + 1;

    end
    

    
end  % for each strand


 end







