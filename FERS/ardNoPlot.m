function [y,ard_] = ardNoPlot(s1,s2,fs,fd_max,td_max,index,ard)


N=length(s1);                                                %number of points
Ndelay = floor(td_max*fs);                                   %number of points corresponding to td_max
Ndop = ceil(N*fd_max/fs);                                    %number of points corresponding to fd_max


%initialisation of temporary variables
s2_pad = [zeros(Ndelay,1);s2];
y1=zeros(Ndelay+1,2*fd_max+1);

%{
tic
for k = 1:Ndelay+1
    temp = s1.*conj(s2_pad(Ndelay+2-k:N+Ndelay+1-k));       %dot-product of the reference and scattered signals
    temp = temp.*hanning(N);                                %windowing the result |Using chebyshev window
    temp2 = fftshift(fft(temp,N));                          %FFT of the above dot-product
    y1(k,:) = temp2(floor(N/2)+1-Ndop:floor(N/2)+1+Ndop);   %Discarding frequency bins not of interest
end
%display('Range-Doppler computation')
toc
%}
%Faster execution time
chebWin = chebwin(N, 100);                       % Precompute Chebyshev window
s1_gpu = gpuArray(s1);                          % Move data to GPU
s2_pad_gpu = gpuArray(s2_pad);                  % Move data to GPU
chebWin_gpu = gpuArray(chebWin);                % Move Chebyshev window to GPU
y1_gpu = gpuArray.zeros(Ndelay+1, 2*Ndop+1);    % Preallocate y1 on GPU

% Loop over delays, use GPU and FFT
for k = 1:Ndelay+1
    % Dot-product of reference and shifted scattered signals (on GPU)
    temp = s1_gpu .* conj(s2_pad_gpu(Ndelay+2-k:N+Ndelay+1-k));
    
    % Apply Chebyshev window
    temp = temp .* chebWin_gpu;
    
    % Perform FFT and select frequency bins of interest
    temp2 = fftshift(fft(temp, N));
    y1_gpu(k,:) = temp2(floor(N/2)+1-Ndop:floor(N/2)+1+Ndop);
end

% Bring result back from GPU
y1 = gather(y1_gpu);

y = abs(y1).^2;                                             %Power conversion
y =y./max(max(abs(y)));                                     %Normalizing max to 1


%figure('Name','2D image');

if index==1
    ard =y ;
end

if index>1
    ard = y ;
end

ard_ = ard;

