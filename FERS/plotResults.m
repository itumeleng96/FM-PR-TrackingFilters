clc; clear all; close all;
addpath('../FERS/','../CFAR/','../KmeansCentroids');

system("fers ../FERS/Simulation_60_direct.fersxml");
system("fers ../FERS/Simulation_60_echo_2.fersxml");
%system("fers ../FERS/Simulation_60_Bistatic.fersxml");
%system("fers ../FERS/singleFile.fersxml");


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
delay = 133e-6;

s1 = I_Qmov;
s2 = I_Qno;

initial=1;
current=fs;  %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs ; %Simulation time: number of data points/sampling frequency
disp(["SimTime:",simulation_time])
%disp(["I_Qmov",size(I_Qmov,1)]);
%disp(["fs",fs]);
%disp(["s1,s2:",size(s1),size(s2)]);

figure('Name','2D image');

for i = 1:simulation_time
    s1 = I_Qmov(initial:current);
    s2 = I_Qno(initial:current);
    ard_plot(s1,s2,fs,dopp_bins,delay);

    initial = current+1;
    current = current + fs;

    disp(["Time(s):",i]);
    %disp(["Current:",current]);
end