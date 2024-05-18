clc; clear all; close all;

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

% Load and save fers output as csv
direct = loadfersHDF5('direct.h5');
echo = loadfersHDF5('echo.h5');

csvwrite('directRx.csv',direct);
csvwrite('echoRx.csv',echo);

fprintf('/nCSV write complete...')

% h5Import
%[Ino Qno scale_no] = loadfersHDF5('direct.h5');
%[Imov Qmov scale_mov] = loadfersHDF5('echo.h5');
%I_Qmov = Imov + j*Qmov;
%I_Qmov = I_Qmov.*scale_mov;
%I_Qno = Ino + j*Qno;
%I_Qno = I_Qno.*scale_no;
%I_Qmov=I_Qmov-I_Qno;
