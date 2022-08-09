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
lambda = ?;
%initialize tracking filter and run every second
EKF_object = EKF(.1,lambda, 1, 1, 1, 0.01,0.01);
X_predicted =[];  %Arrary to store kalman predicted values
X_estimated =[];  %Arrary to store kalman estimated values
Centroids = [];
ard = [];
cfar = [];

for i = 1:10
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    [y,EKF_object_,X_predicted_,X_estimated_,Centroids_,ard_,cfar_] = ardPlotEKF(s1,s2,fs,dopp_bins,delay,EKF_object,X_predicted,X_estimated,Centroids,i,ard,cfar);
    
    X_predicted = X_predicted_;
    X_estimated = X_estimated_;
    EKF_object = EKF_object_;
    Centroids = Centroids_;
    
    ard = ard_;
    cfar = cfar_;
    
    initial = current+1;
    current = current + fs;
end