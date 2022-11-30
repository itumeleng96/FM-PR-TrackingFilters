%function [y,EKF_object_final,X_predict_arr_,X_estimate_arr_,Centroids_arr_,ard_,cfar_] = ardPlotEKF(s1,s2,fs,fd_max,td_max,EKF_object,X_predict_arr,X_estimate_arr,Centroids_arr,index,ard,cfar)
function [y,ard_,cfar_,tracks_,coefficients_,X_predicted_] = ardPlotGN(s1,s2,fs,fd_max,td_max,index,ard,cfar,tracks,X_predicted,coefficients)

%Parameters to allow zoom in
xlim_upper = 0.45e-4;
ylim_upper = 200;
ylim_lower = 80;

%Parameters for filter
NumberOfTargets=1;
initialValues = [[16;0;370;0],[14;0;310;0]];


c = 3e8;                     %speed of the light
N=length(s1);                %number of points
Ndelay = floor(td_max*fs);   %number of points corresponding to td_max
Ndop = ceil(N*fd_max/fs);    %number of points corresponding to fd_max


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
dt = 1;       
tic

f=figure(1);
%figure('Name','2D image');

if index==1
    ard =y ;
end

if index>1
    ard = ard+y ;
end

ard_ = ard;

imagesc(time,frequency,10*log10(ard.'),[max_dB-100 max_dB]);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Range-Doppler response')
display('imagesc plot computation')
movegui(f,'northwest');
xlim([0 1e-4]) 
ylim([0 200])
text(0,10,"Time:" + index+ "s");
drawnow


%GET CFAR and PLOT
f2=figure(2);
[RDM] = ca_cfar(10*log10(y.'),0.5);
if index==1
    cfar =RDM ;
end

if index>1
    cfar = cfar+RDM ;
end

cfar_ = cfar;

imagesc(time,frequency,cfar);
text(0,10,"Time:" + index+ "s");
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('CFAR');
disp('CFAR Plot');
movegui(f2,'northeast');
xlim([0 1e-4]) 
ylim([0 200])
drawnow

%GET Centroids using Kmeans Algorithm
[row,column] = find(RDM>0);
row= row.';
column = column.';
points = [column; row];

%PLOT Centroids and Kalman Estimate and Prediction

f3=figure(3);
imagesc(time,frequency,RDM*0);
text(0,ylim_lower+10,"Time:" + index+ "s");
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Target Centroids and EKF Prediction');
display('Target Centroids and EKF Prediction');

[cluster,centr] = kMeans(NumberOfTargets,points);

%Assign centroid to track
[tracks_]=trackAssign(tracks,centr);

%Based on the Number of Targets get the predictions
for j=1:NumberOfTargets
     GN_Object = GaussNewton(100,10^-5,10^-9,coefficients);
     coefficients_= GN_Object.fit(tracks_(1,:,j),tracks_(2,:,j),coefficients);
     %Get Predictions

end
X_predicted_=X_predicted;

for k=1:NumberOfTargets
    hold on;    
    plot(time((round(tracks_(1,:,k)))),frequency(round(tracks_(2,:,k))),'^-','MarkerFaceColor','black', 'MarkerSize', 5);
    %Plot kalman estimates
    %hold on;
    %plot(time((round(X_predicted_(1,:,k)))),frequency(round(X_predicted_(2,:,k))), 'y-o ', 'MarkerSize', 6); 
end
legend('Quadratic fit','Target Centroids');
movegui(f3,'southwest');
xlim([0 xlim_upper]) 
ylim([ylim_lower ylim_upper])
drawnow


