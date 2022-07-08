clc; clear all; close all;

% h5Import
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
current=500000;

%initialize tracking filter and run every second
KF_object = kalmanFilter(1/fs, 1, 1, 1, 0.1,0.1);

for i = 1:9
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    [y,KF_object_] = ard_plot(s1,s2,fs,dopp_bins,delay,KF_object);
    
    KF_object = KF_object_;
    
    initial = current+1;
    current = current + fs;
end