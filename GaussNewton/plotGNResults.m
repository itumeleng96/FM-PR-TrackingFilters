clc; clear all; close all;
addpath('FERS/','CFAR/','KmeansCentroids');

system("fers ../FERS/Simulation_60_direct.fersxml");
system("fers ../FERS/Simulation_60_echo.fersxml");

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
X_predicted =[];  %Arrary to store kalman predicted values
Centroids = [];
ard = [];
cfar = [];
tracks =[];

initialValues = [16;0;0;370;0;0];
max_iterations=100;
tolerance = 1e-6;

GN_Object = GaussNewton(0.1, 0.1, 0.1, 1, 0.01,0.01,initialValues,max_iterations,tolerance);


for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    %[y,EKF_object_,X_predicted_,X_estimated_,Centroids_,ard_,cfar_] = ardPlotEKF(s1,s2,fs,dopp_bins,delay,EKF_object,X_predicted,X_estimated,Centroids,i,ard,cfar);
    [y,ard_,cfar_,tracks_,GN_Object_,X_predicted_] = ardPlotGN(s1,s2,fs,dopp_bins,delay,i,ard,cfar,tracks,X_predicted,GN_Object);
    
    X_predicted = X_predicted_;
    ard = ard_;
    cfar = cfar_;
    tracks = tracks_;
    GN_Object = GN_Object_;
    initial = current+1;
    current = current + fs;
end