function [particles,weights,X_estimate_arr_,Centroids_arr_,ard_,cfar_] = ardPlotParticle(s1,s2,fs,fd_max,td_max,particles,weights,index,X_estimate_arr,Centroids_arr,ard,cfar)

  
c = 3e8;                     %speed of the light
N=length(s1);                %number of points
Ndelay = floor(td_max*fs);   %number of points corresponding to td_max
Ndop = ceil(N*fd_max/fs);    %number of points corresponding to fd_max


%initialisation of temporary variables
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

figure(1);
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
title('Range-Doppler response');
display('imagesc plot computation');
drawnow


%GET CFAR and Plot
figure(2);
[RDM] = ca_cfar(10*log10(y.'),0.5);
if index==1
    cfar =RDM ;
end

if index>1
    cfar = cfar+RDM ;
end

cfar_ = cfar;

imagesc(time,frequency,cfar);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('CFAR');
display('imagesc plot computation');


%GET Centroids using Kmeans Algorithm and Plot
[row,column] = find(RDM>0);
row= row.';
column = column.';
points = [column; row];

[cluster,centr] = kMeans(1,points);
%Save previous Centroids
Centroids_arr(1,index) = centr(1,:) ;
Centroids_arr(2,index) = centr(2,:) ;
Centroids_arr_ = Centroids_arr;


%Particle filter iteration
U = [1,1]; %input
std= [2,5];
dt = 0.000005;
meas_err = 1;
%Predict and Plot
particles =predict(particles,U,std,dt);


figure(3);
imagesc(time,frequency,RDM*0);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Target Centroids and Particle Filter Estimation');
display('Target Centroids and Particle Filter  Estimation');

hold on
plot(time((round(Centroids_arr(1,:)))),frequency(round(Centroids_arr(2,:))),'+','MarkerFaceColor','black', 'MarkerSize', 8);
text(0,0,"Time:" + index+ "s");


%hold on;
%plot(time((round(particles(:,1)))),frequency(round(particles(:,2))), 'r-o', 'MarkerSize', 5);
%legend('Target Centroid','Particle Filter Prediction');

drawnow
%incorporate Measurement and update
[particles,weights] = update(particles,weights,[centr(1,:),centr(2,:)],meas_err);

%resample if too few effective particles
neff = NEFF(weights);
if neff< N/2 
    indexes = resampleSystematic(weights);
    [particles,weights]= resampleFromIndex(particles,indexes);
end

%Mean and variance from Particles
[mean,var] = estimate(particles,weights);

display([centr(1,:),centr(2,:)]);
display(mean(:));

%Save previous Particle Filter estimates
X_estimate_arr(1,index) = mean(1) ;
X_estimate_arr(2,index) = mean(2) ;
X_estimate_arr_ = X_estimate_arr;

figure(4);
imagesc(time,frequency,RDM*0);
axis xy;
colorbar;
xlabel('Bistatic delay [s]','Fontsize',10);
ylabel('Doppler frequency [Hz]','Fontsize',10);
grid on;
title('Target Centroids and Particle Filter Estimation');
display('Target Centroids and Particle Filter  Estimation');

hold on
plot(time((round(Centroids_arr(1,:)))),frequency(round(Centroids_arr(2,:))),'^-','MarkerFaceColor','black', 'MarkerSize', 5);
text(0,0,"Time:" + index+ "s");

hold on;
plot(time((round(X_estimate_arr(1,:)))),frequency(round(X_estimate_arr(2,:))), 'r-o', 'MarkerSize', 5);
legend('Target Centroid','Particle Filter Estimate');

drawnow
toc


