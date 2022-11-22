clc; clear all; close all;
addpath('../FERS/','../CFAR/','../KmeansCentroids');

system("fers ../FERS/Simulation_60_direct.fersxml");
system("fers ../FERS/Simulation_60_echo_2.fersxml");


% h5Import from FERS simulation
[Ino Qno scale_no] = loadfersHDF5('direct.h5');
[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');


I_Qmov = Imov + j*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + j*Qno;
I_Qno = I_Qno.*scale_no;

I_Qmov=I_Qmov-I_Qno;

% run_ard
fs = 200000;
dopp_bins = 200;
delay = 233e-6;

s1 = I_Qmov;
s2 = I_Qno;

initial=1;
current=fs;  %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs ; %Simulation time: number of data points/sampling frequency

%initialize tracking filter and run every second
EKF_objects =[];
X_predicted =[];  %Arrary to store kalman predicted values
X_estimated =[];  %Arrary to store kalman estimated values
Centroids = [];
ard = [];
cfar = [];
tracks =[];

figure('Name','2D image');
%ARD
f=figure(1);
f.Position = [4000 10 1000 800]; 

%CFAR
f2=figure(2);
f2.Position = [4000 10 1000 800]; 

%Prediction
f3=figure(3);
f3.Position = [4000 10 1000 800];

%Estimation
f4=figure(4);
f4.Position = [4000 10 1000 800]; 

for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);

    [y,ard_,cfar_,tracks_,EKF_objects_,X_predicted_,X_estimated_] = ardPlotEKF(s1,s2,fs,dopp_bins,delay,i,ard,cfar,tracks,EKF_objects,X_predicted,X_estimated,f,f2,f3,f4);
     
    X_predicted = X_predicted_;
    X_estimated = X_estimated_;
    EKF_objects = EKF_objects_;    
    ard = ard_;
    cfar = cfar_;
    tracks = tracks_;
    
   
    initial = current+1;
    current = current + fs;
end