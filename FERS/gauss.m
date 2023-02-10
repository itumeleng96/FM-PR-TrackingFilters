
sig = randn(1,200000*60)+1i*randn(1,200000*60);
I = real(sig);
Q = imag(sig);
fprintf('Writing HDF5 data...\n');
hdf5write('waveform.h5', '/I/value', real(I), '/Q/value', imag(Q));
fprintf('Complete.\n');

%NB:Using Gauss for Transmitted Signal, waiting for FM waverform from Dr Paine
