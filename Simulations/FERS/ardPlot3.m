function [y,ard_] = ardPlot3(s1,s2,fs,fd_max,td_max,index,ard,f)

%Parameters to allow zoom in
xlim_upper = 2e-4;
ylim_upper = 200;
ylim_lower = -200;


c = 299792458;                                               %speed of the light
N=length(s1);                                                %number of points
Ndelay = floor(td_max*fs);                                   %number of points corresponding to td_max
Ndop = ceil(N*fd_max/fs);                                    %number of points corresponding to fd_max


%initialisation of temporary variables
s2_pad = [zeros(Ndelay,1);s2];
y1=zeros(Ndelay+1,2*fd_max+1);


tic
for k = 1:Ndelay+1
    temp = s1.*conj(s2_pad(Ndelay+2-k:N+Ndelay+1-k));       %dot-product of the reference and scattered signals
    temp = temp.*chebwin(N,60);                                %windowing the result |Using chebyshev window
    temp2 = fftshift(fft(temp,N));                          %FFT of the above dot-product
    y1(k,:) = temp2(floor(N/2)+1-Ndop:floor(N/2)+1+Ndop);   %Discarding frequency bins not of interest
end
%display('Range-Doppler computation')
toc


y = abs(y1).^2;                                             %Power conversion
y =y./max(max(abs(y)));                                     %Normalizing max to 1

%Time and frequency axis
time = 0:1/fs:Ndelay/fs;
range = time*c;
range = range/1000;
frequency = -fd_max:1:fd_max;
Dyn_dB = 40;                                                %Dynamic range (dB)
max_dB = 10*log10(max(max(abs(y))));
tic


%figure('Name','2D image');

if index==1
    ard =y ;
end

if index>1
    ard = y ;
end

ard_ = ard;

figure(f);

[R, F] = meshgrid(range, frequency);
contourf(R, F, 10 * log10(ard.'), 20, 'LineColor', 'none');
colormap(jet);
colorbar; % To add a colorbar for reference

axis xy;
ax = gca; 
ax.FontSize = 14;

colorbar;
xlabel('Bistatic  Range [km]','Fontsize',20);
ylabel('Bistatic Doppler frequency [Hz]','Fontsize',20);
grid on;
c = colorbar;
c.Label.String = 'Level [dB]';
c.Label.FontSize = 18; 

title('Range-Doppler response','FontSize', 20);
%display('imagesc plot computation')
%xlim([0 xlim_upper]) 
ylim([ylim_lower ylim_upper])
%text(0,10,"Time:" + index+ "s");
drawnow