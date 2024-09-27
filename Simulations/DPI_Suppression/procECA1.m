function [SurvData] = procECA1(RefData, SurvData, proc)
fprintf('Performing ECA cancellation with GPU acceleration\n')

% The number of samples per segment
segmentSize_nSamp = floor(max(size(SurvData)) / proc.nSegments);

% The number of range and Doppler bins over which cancellation will be applied
nRangeBins = ceil((proc.cancellationMaxRange_m - proc.TxToRefRxDistance_m)/(3e8/proc.Fs));
nDopplerBinds = ceil(proc.cancellationMaxDoppler_Hz/(proc.Fs/segmentSize_nSamp)) * 2 + 1;

% Precompute Doppler shifts to save time in the loop
sampleTimes_s = gpuArray(0:1/proc.Fs:(segmentSize_nSamp - 1)/proc.Fs);
dopplerShifts = exp((-floor(nDopplerBinds/2):floor(nDopplerBinds/2))' * 1j * 2 * pi * sampleTimes_s);

% Move surveillance data to the GPU before loop
SurvData_gpu = gpuArray(SurvData);

% Loop over segments
for segmentNo = 0:proc.nSegments - 1
    % Start and end sample numbers for this segment
    segStartSampleNo = segmentNo * segmentSize_nSamp + 1;
    segStopSampleNo = segStartSampleNo + segmentSize_nSamp - 1;

    % Initialize ZeroDopplerA directly on GPU
    ZeroDopplerA_gpu = gpuArray.zeros(segmentSize_nSamp, nRangeBins, 'single');

    % Create ZeroDoppler A matrix from the reference data (now on GPU)
    for i = 1:nRangeBins
        ZeroDopplerA_gpu(i:segmentSize_nSamp, i) = RefData((segStartSampleNo:segStopSampleNo - (i - 1)));
    end

    % Initialize the full A matrix on the GPU
    A_gpu = gpuArray.zeros(segmentSize_nSamp, nRangeBins * nDopplerBinds, 'single');
    Pos = 0;

    % Create Doppler-shifted versions of ZeroDopplerA on the GPU
    for i = -floor(nDopplerBinds / 2):floor(nDopplerBinds / 2)
        for l = 1:nRangeBins
            A_gpu(:,l + Pos) = ZeroDopplerA_gpu(:,l) .* dopplerShifts(i + floor(nDopplerBinds / 2) + 1, :)';
        end
        Pos = Pos + nRangeBins;
    end

    % Perform least-squares calculation on the GPU
    SurvData_seg_gpu = SurvData_gpu(segStartSampleNo:segStopSampleNo);
    alpha_gpu = (A_gpu' * A_gpu) \ (A_gpu' * SurvData_seg_gpu);

    % Subtract the result from the surveillance data (on GPU)
    SurvData_gpu(segStartSampleNo:segStopSampleNo) = SurvData_seg_gpu - (A_gpu * alpha_gpu);

    % Clear unnecessary variables to free up GPU memory
    clear A_gpu ZeroDopplerA_gpu alpha_gpu
end

% Retrieve the modified surveillance data from the GPU back to the CPU
SurvData = gather(SurvData_gpu);

clear sampleTimes_s dopplerShifts

end
