function [y,KF_object_final,X_predict_arr_,X_estimate_arr_,Centroids_arr_] = ardPlotKalman(s1,s2,fs,fd_max,td_max,KF_object,X_predict_arr,X_estimate_arr,Centroids_arr,index)

c = 3e8;                    %speed of the light
N=length(s1);               %number of points
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

tic

%Predict using Kalman Filter
[X,KF_object_]= predict(KF_object);
KF_object=KF_object_;
%Save previous Kalman Estimates
X_predict_arr(index,1) = X(1,1) ;
X_predict_arr(index,2) = X(3,1) ;
X_predict_arr_ = X_predict_arr;


figure(1);
%figure('Name','2D image');
imagesc(time,frequency,10*log10(y.'),[max_dB-100 max_dB]);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Range-Doppler response')
display('imagesc plot computation')
drawnow


%GET CFAR and PLOT
figure(2);
[RDM] = ca_cfar(10*log10(y.'),0.5);
imagesc(time,frequency,RDM);
text(0,0,"Time:" + index+ "s");
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('CFAR');
display('CFAR Plot');
drawnow

%GET Centroids using Kmeans Algorithm
[row,column] = find(RDM>0);
row= row.';
column = column.';
points = [column; row];

%PLOT Centroids and Kalman Estimate and Prediction
figure(3);
imagesc(time,frequency,RDM*0);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Target Centroids and Kalman Prediction');
display('Target Centroids and Kalman Prediction');

%Plot kalman estimates
hold on;
plot(time((round(X_predict_arr(:,1)))),frequency(round(X_predict_arr(:,2))), 'y-o ', 'MarkerSize', 8);

[cluster,centr] = kMeans(1,points);

%Save previous Centroids
Centroids_arr(1,index) = centr(1,:) ;
Centroids_arr(2,index) = centr(2,:) ;
Centroids_arr_ = Centroids_arr;

hold on;
plot(time((round(Centroids_arr(1,:)))),frequency(round(Centroids_arr(2,:))),'^-','MarkerFaceColor','black', 'MarkerSize', 5);
text(0,0,"Time:" + index+ "s");
legend('Kalman Precition','Target Centroids');

%update Kalman filter

figure(4);
[X1,KF_object1] = update(KF_object,[centr(1,:);centr(2,:)]);
imagesc(time,frequency,RDM*0);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Target Centroids and Kalman Estimation');
display('Target Centroids and Kalman Estimation');

%Save previous Kalman estimates
X_estimate_arr(1,index) = X1(1,1) ;
X_estimate_arr(2,index) = X1(3,1) ;
X_estimate_arr_ = X_estimate_arr;

hold on;
plot(time((round(Centroids_arr(1,:)))),frequency(round(Centroids_arr(2,:))),'^-','MarkerFaceColor','black', 'MarkerSize', 5);
text(0,0,"Time:" + index+ "s");

hold on;
plot(time((round(X_estimate_arr(1,:)))),frequency(round(X_estimate_arr(2,:))), 'r-o', 'MarkerSize', 5);
KF_object_final = KF_object1;
legend('Target Centroid','Kalman Estimate');


drawnow
toc


