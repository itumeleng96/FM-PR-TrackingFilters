sig = randn(1,5000000)+j*randn(1,5000000);
I = real(sig);
Q = imag(sig);
fprintf('Writing HDF5 data...\n');
hdf5write('waveform.h5', '/I/value', real(I), '/Q/value', imag(Q));
fprintf('Complete.\n');