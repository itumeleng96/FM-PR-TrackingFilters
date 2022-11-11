%function [y,EKF_object_final,X_predict_arr_,X_estimate_arr_,Centroids_arr_,ard_,cfar_] = ardPlotEKF(s1,s2,fs,fd_max,td_max,EKF_object,X_predict_arr,X_estimate_arr,Centroids_arr,index,ard,cfar)
function [y,ard_,cfar_,tracks_,EKF_objects_,X_predicted_,X_estimated_] = ardPlotEKF(s1,s2,fs,fd_max,td_max,index,ard,cfar,tracks,EKF_objects,X_predicted,X_estimated)

%Parameters to allow zoom in
xlim_upper = 2e-4;
ylim_upper = 200;
ylim_lower = -200;

%Parameters for filter
NumberOfTargets=1;
initialValues = [[17;0;0;330;0;0],[14;0;0;310;0;0]];

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
%Based on the Number of Targets get the predictions
for j=1:NumberOfTargets
    if index ==1
        EKF_object = EKF(dt, 0.1, 0.1, 1, 0.01,0.01,initialValues(:,j));
        [X,EKF_object_]= predict(EKF_object);
        EKF_objects_(j)=EKF_object_;
        X_predicted(1,index,j) =X(1,1);
        X_predicted(2,index,j) =X(4,1);
    else
        [X,EKF_object_]= predict(EKF_objects(j));
        EKF_objects_(j)=EKF_object_;
        X_predicted(1,index,j) =X(1,1);
        X_predicted(2,index,j) =X(4,1);
    end
end
X_predicted_=X_predicted;

f=figure(1);
f.Position = [4000 10 1000 800]; 
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
xlim([0 xlim_upper]) 
ylim([ylim_lower ylim_upper])
text(0,10,"Time:" + index+ "s");
drawnow


%GET CFAR and PLOT
f2=figure(2);
f2.Position = [4000 10 1000 800]; 

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
display('CFAR Plot');
movegui(f2,'northeast');
xlim([0 xlim_upper]) 
ylim([ylim_lower ylim_upper])
drawnow

%GET Centroids using Kmeans Algorithm
[row,column] = find(RDM>0);
row= row.';
column = column.';
points = [column; row];

%PLOT Centroids and Kalman Estimate and Prediction

f3=figure(3);
f3.Position = [4000 10 1000 800]; 
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

for k=1:NumberOfTargets
    hold on;    
    plot(time((round(tracks_(1,:,k)))),frequency(round(tracks_(2,:,k))),'^-','MarkerFaceColor',	[0 0 0], 'MarkerSize', 7);
    %Plot kalman estimates
    %hold on;
    %plot(time((round(X_predicted_(1,:,k)))),frequency(round(X_predicted_(2,:,k))), 'o- ','MarkerFaceColor',[1 0 0], 'MarkerSize', 8); 
end
legend('Target Centroids','Kalman Prediction');
movegui(f3,'southwest');
xlim([0 xlim_upper]) 
ylim([ylim_lower ylim_upper])
drawnow
%update Kalman filter

f4=figure(4);
f4.Position = [4000 10 1000 800]; 

imagesc(time,frequency,RDM*0);
text(0,ylim_lower+10,"Time:" + index+ "s");
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Target Centroids and EKF Estimation');
display('Target Centroids and EKF Estimation');

for n=1:NumberOfTargets
    [X1,EKF_objects_(n)] = update(EKF_objects_(n),[tracks_(1,index,n);tracks_(2,index,n)]); 
    %Save previous Kalman estimates
    X_estimated(1,index,n) = X1(1,1);
    X_estimated(2,index,n) = X1(4,1);    
end

X_estimated_ = X_estimated; 

for l=1:NumberOfTargets
    hold on;
    plot(time((round(tracks_(1,:,l)))),frequency(round(tracks_(2,:,l))),'^-','MarkerFaceColor','black', 'MarkerSize', 5);

    hold on;
    plot(time((round(X_estimated_(1,:,l)))),frequency(round(X_estimated_(2,:,l))), 'y-o ', 'MarkerSize', 6); 
end

legend('Target Centroid','Kalman Estimate');
xlim([0 xlim_upper]) 
ylim([ylim_lower ylim_upper])
drawnow
toc
