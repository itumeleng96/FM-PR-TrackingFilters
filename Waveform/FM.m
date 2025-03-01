fprintf('\n----RECORDED DATA----\n');
   

oRCF = cRCF;
filename = 'Malmesbury_1.rcf';
    
% Select where to start from in the recorded data (This data is about 10 min long so anything from 0 to 500s+ works)     

t_start = 60;
t_length_Tx = 120;
t_delay = 0;
Fs = 200000;

refORsur = input('Reference [0] or surveillance channel [1]: ');
name = input('Output file name [RecordedNormalised.h5]: ','s');

if isempty(name) == true
	name = 'RecordedNormalised.h5';
end

if isempty(refORsur) == true
    refORsur = 0;
end

if isempty(t_delay) == true
    t_delay = 60;
end

t_length = t_length_Tx - t_delay;
t_end = t_start+t_length_Tx;

oRCF = oRCF.readFromFile(filename, (t_start+1)*Fs, t_end*Fs); %180 seconds at a sample rate of 204800 Hz (180*204800=36864000)

if refORsur == 0
	%Normalise ref data
	ref = oRCF.getReferenceData;

	refPower = var(ref);
	ref = ref .* (1 / sqrt(refPower));

	% Pad with zeros if requested by user    
	if t_delay > 0
	    ref_new = [zeros(t_delay*Fs,1); ref];
	else
	    ref_new = ref;
	end  

	fprintf('\nWriting HDF5 data...\n');
	hdf5write(name, '/I/value', real(ref_new), '/Q/value', imag(ref_new));
	fprintf('Complete.\n');
    
	signalDuration = t_length; % To allow plotting
    L = length(ref_new);
    Y = fftshift(fft(ref_new));
    f = Fs*(-L/2:(L/2-1))/L;
    plot(f,abs(Y));
    title('Frequency Domain of the Signal');
    xlabel('Frequency (Hz)');
    ylabel('|Y(f)|');

elseif refORsur == 1
	%Normalise sur data
	ref = oRCF.getSurveillanceData;

	refPower = var(ref);
	ref = ref .* (1 / sqrt(refPower));

	% Pad with zeros if requested by user    
	if t_delay > 0
	    ref_new = [zeros(t_delay*Fs,1); ref];
	else
	    ref_new = ref;
	end  

	fprintf('\nWriting HDF5 data...\n');
	hdf5write(name, '/I/value', real(ref_new), '/Q/value', imag(ref_new));
	fprintf('Complete.\n');
    
	signalDuration = t_length; % To allow plotting
end
