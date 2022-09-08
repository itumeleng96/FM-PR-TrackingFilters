clc; clear all; close all;
addpath('FERS/','CFAR/','KmeansCentroids');

% h5Import from FERS simulation
[Ino Qno scale_no] = loadfersHDF5('direct.h5');
[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');


I_Qmov = Imov + j*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + j*Qno;
I_Qno = I_Qno.*scale_no;

I_Qmov=I_Qmov-I_Qno;

% run_ard
fs = 500000;
dopp_bins = 200;
delay = 233e-6;

s1 = I_Qmov;
s2 = I_Qno;

initial=1;
current=500000;  %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs ; %Simulation time: number of data points/sampling frequency

%initialize tracking filter and run every second
EKF_objects =[];
X_predicted =[];  %Arrary to store kalman predicted values
X_estimated =[];  %Arrary to store kalman estimated values
Centroids = [];
ard = [];
cfar = [];
tracks =[];

for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    %[y,EKF_object_,X_predicted_,X_estimated_,Centroids_,ard_,cfar_] = ardPlotEKF(s1,s2,fs,dopp_bins,delay,EKF_object,X_predicted,X_estimated,Centroids,i,ard,cfar);
    [y,ard_,cfar_,tracks_,EKF_objects_,X_predicted_,X_estimated_] = ardPlotGN(s1,s2,fs,dopp_bins,delay,i,ard,cfar,tracks,EKF_objects,X_predicted,X_estimated);
    
    X_predicted = X_predicted_;
    X_estimated = X_estimated_;
    EKF_objects = EKF_objects_;    
    ard = ard_;
    cfar = cfar_;
    tracks = tracks_;
    
    initial = current+1;
    current = current + fs;
end