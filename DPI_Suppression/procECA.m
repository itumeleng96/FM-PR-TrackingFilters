function [SurvData] = procECA(RefData, SurvData, proc)
fprintf('Performing ECA cancellation\n')

%%%%
% USAGE:
%
% Input:
%	RefData: Reference channel time data
%	SurvData: Surveillance channel time data
%
%	proc.cancellationMaxRange_m = 12650; 
%	proc.cancellationMaxDoppler_Hz = 4;
%	proc.TxToRefRxDistance_m = 12600;
%	proc.nSegments = 16;
%	proc.nIterations = 30;
%	proc.Fs = Fs;
%	proc.alpha = 0;
%	proc.initialAlpha = 0;
%
% Output:
%	SurvData: Cancelled surveillance channel time data
%%%%

%The number of samples per segment
segmentSize_nSamp = floor(max(size(SurvData)) /  proc.nSegments);

%The number of range and Doppler bins over which cancellation will be applied
nRangeBins = ceil((proc.cancellationMaxRange_m - proc.TxToRefRxDistance_m)/(3e8/proc.Fs));
nDopplerBinds = ceil(proc.cancellationMaxDoppler_Hz/(proc.Fs/segmentSize_nSamp)) * 2 + 1;

%The delay in seconds to each sample from the begining of the CPI
sampleTimes_s = 0:1/proc.Fs:(segmentSize_nSamp - 1)/proc.Fs;

for segmentNo = 0:proc.nSegments - 1
    
%     fprintf('\t\t* ECA: Starting cancellation on segment %i of %i:\n', segmentNo, proc.nSegments)
%     fprintf('\t\t* Range:   %g to %g m, %g bins\n', proc.TxToRefRxDistance_m, proc.cancellationMaxRange_m, nRangeBins)
%     fprintf('\t\t* Doppler: %g to %g Hz, %g bins\n', -proc.cancellationMaxDoppler_Hz, proc.cancellationMaxDoppler_Hz, nDopplerBinds)
%     tic
    
    %Start and end sample numbers for this segment
    segStartSampleNo = segmentNo * segmentSize_nSamp + 1;
    segStopSampleNo = segStartSampleNo + segmentSize_nSamp - 1;
    
    %A matrix (As in Ax = b) with only the zero Doppler rows
    ZeroDopplerA=zeros(segmentSize_nSamp, nRangeBins, 'single');
    
    %Create ZeroDoppler A matrix from RCF object.
    %Each row is the surveillance channel left zero padded 
    for i=1:nRangeBins
        ZeroDopplerA(i:segmentSize_nSamp, i)=RefData((segStartSampleNo:segStopSampleNo - (i - 1)));
    end
    
    %The complete A matrix
    A = zeros(segmentSize_nSamp, nRangeBins * nDopplerBinds, 'single');
    K = nRangeBins;
    Pos = 0;
    
    %Create Doppler shifted versions of ZeroDopplerA in  A
    for i=-floor(nDopplerBinds / 2):floor(nDopplerBinds / 2)
        for l=1:K
            A(:,l+Pos)=ZeroDopplerA(:,l).* exp(i * 1j * 2 * pi * sampleTimes_s');
        end
        Pos = Pos + K;
    end
    
    clear ZeroDopplerA
    
    alpha = (A'*A)\A'*SurvData(segStartSampleNo:segStopSampleNo);
    
    SurvData(segStartSampleNo:segStopSampleNo) = SurvData(segStartSampleNo:segStopSampleNo) - (A * alpha);
    
    clear d A Pos K x i l alpha
%     fprintf('\t\tCompleted. ')
%     toc
%     fprintf('\n')

end

clear sampleTimes_s


