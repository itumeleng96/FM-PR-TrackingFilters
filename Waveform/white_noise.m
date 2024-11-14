timeDuration_s = 60;
Fs_Hz = 200000;

vSamples = randn(timeDuration_s * Fs_Hz, 1) + 1j * randn(timeDuration_s * Fs_Hz, 1);

%Normalise data
power = var(vSamples);
vSamples = vSamples .* (1 / sqrt(power));

plot(abs(fftshift(fft(vSamples))))

fprintf('Writing HDF5 data...\n');
hdf5write('white_noise.h5', '/I/value', real(vSamples), '/Q/value', imag(vSamples));
fprintf('Complete.\n');
