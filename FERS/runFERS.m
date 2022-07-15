%Remove previous data
clear; clc; close all;

% Ask user for files
%prompt_direct = 'Please enter the direct signal simulation file: ';
%prompt_echo = '\nPlease enter the echo signal simulation file: ';
%prompt_waveform = '\nPlease enter waveform for file simulation: ';

%direct = input(prompt_direct, 's');
%fprintf('Done')
%echo = input(prompt_echo, 's');
%fprintf('Done')
%waveform = input(prompt_waveform, 's');
%fprintf('Done')

direct = 'Simulation_direct.fersxml'
fprintf('Done')
echo = 'Simulation_echo.fersxml'
fprintf('Done')

% Run the simulation for direct signal
fprintf('\nSimulating Reference Signal\n\n')
strCommand = sprintf('export LD_LIBRARY_PATH="" && fers %s',direct);
system(strCommand, '-echo');
fprintf('Done')

% Run the simulation for echo signal
fprintf('\nSimulating Surveilance Signal\n\n')
strCommand = sprintf('export LD_LIBRARY_PATH="" && fers %s',echo);
system(strCommand, '-echo');
fprintf('Done')

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
delay = 133e-6;

s1 = I_Qmov(1:500000);
s2 = I_Qno(1:500000);
% 
ard_plot(s1,s2,fs,dopp_bins,delay);
% 
s3 = I_Qmov(500001:1000000);
s4 = I_Qno(500001:1000000);
% 
ard_plot(s3,s4,fs,dopp_bins,delay);

s5 = I_Qmov(1000001:1500000);
s6 = I_Qno(1000001:1500000);
% 
ard_plot(s5,s6,fs,dopp_bins,delay);