function [SurvData] = procECA2(RefData, SurvData, proc)
fprintf('Performing ECA cancellation with GPU acceleration\n')

% Precompute Doppler shifts on the GPU
segmentSize_nSamp = floor(max(size(SurvData)) / proc.nSegments);
nRangeBins = ceil((proc.cancellationMaxRange_m - proc.TxToRefRxDistance_m)/(3e8/proc.Fs));
nDopplerBinds = ceil(proc.cancellationMaxDoppler_Hz/(proc.Fs/segmentSize_nSamp)) * 2 + 1;
sampleTimes_s = gpuArray(0:1/proc.Fs:(segmentSize_nSamp - 1)/proc.Fs);
dopplerShifts = gpuArray(exp((-floor(nDopplerBinds/2):floor(nDopplerBinds/2))' * 1j * 2 * pi * sampleTimes_s));

% Move data to GPU
RefData_gpu = gpuArray(RefData);
SurvData_gpu = gpuArray(SurvData);

% Preallocate A matrix and ZeroDopplerA to reuse memory across segments
A_gpu = gpuArray.zeros(segmentSize_nSamp, nRangeBins * nDopplerBinds, 'single');
ZeroDopplerA_gpu = gpuArray.zeros(segmentSize_nSamp, nRangeBins, 'single');

% Loop over segments with memory reuse
for segmentNo = 0:proc.nSegments - 1
    segStartSampleNo = segmentNo * segmentSize_nSamp + 1;
    segStopSampleNo = segStartSampleNo + segmentSize_nSamp - 1;
    
    % Refill ZeroDopplerA_gpu matrix directly for the current segment
    for i = 1:nRangeBins
        ZeroDopplerA_gpu(i:segmentSize_nSamp, i) = RefData_gpu(segStartSampleNo:segStopSampleNo - (i - 1));
    end
    
    % Fill A_gpu matrix by applying Doppler shifts
    Pos = 0;
    for i = 1:nDopplerBinds
        A_gpu(:, Pos + (1:nRangeBins)) = ZeroDopplerA_gpu .* dopplerShifts(i, :)';
        Pos = Pos + nRangeBins;
    end
    
    % Perform the least-squares calculation and subtraction
    SurvData_seg_gpu = SurvData_gpu(segStartSampleNo:segStopSampleNo);
    alpha_gpu = (A_gpu' * A_gpu) \ (A_gpu' * SurvData_seg_gpu);
    SurvData_gpu(segStartSampleNo:segStopSampleNo) = SurvData_seg_gpu - (A_gpu * alpha_gpu);
end

% Retrieve data back to CPU
SurvData = gather(SurvData_gpu);

end
