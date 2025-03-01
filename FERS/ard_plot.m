%Generates an ARD plot for the reference and scattered signals using the
%frequency domain implementation.

function y = ard_plot(s1,s2,fs,fd_max,td_max)

  %Corresponds to:
  %1.  The dot-product of the reference and scattered signals
  %2.  Windowing the results to minimize the sidelobes in Doppler
  %3.  The FFT of the above dot-product
  %4.  Discarding frequency bins not of interest
  %5.  Delaying the reference signal by one sample and repeating the above process
  %6.  An ARD plot is then generated and displayed

  %Author:  Sebastiaan Heunis, Yoann paichard

  %function ard(s1,s2,fs,fd_max,Td)
  %s1:  Scattered signal, N samples in column [s1(1,1;...;s(N,1)]
  %s2:  Reference signal, N samples in column [s1(1,1;...;s(N,1)]
  %fs:  Sampling frequency
  %fd_max:  Maximum Doppler shift: fd_max < fs
  %td_max:  Maximum time delay: td_max < N/fs

  % Limits and resolution for computation :
  % Doppler  fd : range  = [-fd_max : fd_max], resolution = fs/N
  % Tau : range  = [0 : td_max], resolution = 1/fs  

  c = 3e8;                    %speed of the light
  N=length(s1);               %number of points
  Ndelay = floor(td_max*fs)   %number of points corresponding to td_max
  Ndop = ceil(N*fd_max/fs)    %number of points corresponding to fd_max

  %initialisation of temporary variables
  temp = zeros(N,1);
  temp2 = zeros(N,1);
  s2_pad = [zeros(Ndelay,1);s2];
  y1=zeros(Ndelay+1,2*fd_max+1);

  tic
  for k = 1:Ndelay+1
      temp = s1.*conj(s2_pad(Ndelay+2-k:N+Ndelay+1-k));       %dot-product of the reference and scattered signals
      temp = temp.*hanning(N);                                %windowing the result 
      temp2 = fftshift(fft(temp,N));                          %FFT of the above dot-product
      y1(k,:) = temp2(floor(N/2)+1-Ndop:floor(N/2)+1+Ndop);   %Discarding frequency bins not of interest
  end
  display('Range-Doppler computation')
  toc

  y = abs(y1).^2;             %Power conversion
  y =y./max(max(abs(y)));     %Normalizing max to 1

  %Time and frequency axis
  time = 0:1/fs:Ndelay/fs;
  range = time*c;
  frequency = -fd_max:1:fd_max;
  Dyn_dB = 40;                %Dynamic range (dB)
  max_dB = 10*log10(max(max(abs(y))));
    %{
  tic
  figure('Name','Contour plot');
  contourf(time,frequency,10*log10(y.'));
  colormap jet;
  colorbar;
  grid on;
  xlabel('Bistatic delay [s]','Fontsize',10);
  ylabel('Doppler frequency [Hz]','Fontsize',10);
  title('Range-Doppler response')
  display('contour plot computation')
  toc
    %}
  tic
  %figure('Name','2D image');
  imagesc(range,frequency,10*log10(y.'),[max_dB-100 max_dB]);
  axis xy;
  colorbar;
  %xlabel('Bistatic delay [s]','Fontsize',10);
  xlabel('Bistatic range [m]','Fontsize',10);

  ylabel('Doppler frequency [Hz]','Fontsize',10);
  grid on;
  title('Range-Doppler response')
  display('imagesc plot computation')
  drawnow
  toc