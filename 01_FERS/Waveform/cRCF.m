%********************************************
%## Commensal Radar Project ##
%
%Filename:
%	cRCF.m
%Description
%	Matlab implementantion of the RCF class
%
%Author:
%	Craig Tong
%	Radar Remote Sensing Group
%	Department of Electrical Engineering
%	University of Cape Town
%	craig.tong@uct.ac.za
%	Copyright (C) 2014 Craig Tong
%********************************************

classdef cRCF < matlab.mixin.Copyable
    properties %Should be private but can't because it breaks mex file access.
        m_FileType = '';
        m_TimeStamp_us  = 0;
        m_Fc_Hz = 0;
        m_Fs_Hz = 0;
        m_Bw_Hz = 0;
        m_NSamples = 0; %range dimension
        m_CommentOffset_B = 0;
        m_CommentLength = 0;
        m_FileSize_B = 0;
        m_strComment = '';
        m_strFilename = '';
        m_fvSurveillanceData = complex(zeros(0, 0, 'single')); %data is stored in single precision here
        m_fvReferenceData = complex(zeros(0, 0, 'single'));
    end %end protected variables
    
    properties(SetAccess = public, Constant = true)
        DEFAULT_TCP_PORT = 5002;
        HEADER_SIZE = 52;
    end % const variables
    
    methods
        function oNewRCF = cRCF(oRCF)
            %Default constructor
            if nargin == 0
                oNewRCF.m_FileType = 'rcf';
                oNewRCF.m_FileType(end + 1) = 0; % null terminate string
                
            %Copy contructor
            else
                if(isa(oRCF,'cRCF') && nargin == 1)
                    allProperties = properties(oRCF);
                    for i=1:length(allProperties)
                        try
                            oNewRCF.(allProperties{i}) = oRCF.(allProperties{i});
                        catch
                        end
                    end
                end
            end
        end
        
        %% Mutators
        function oRCF = setFileType(oRCF, FileType)
            oRCF.m_FileType = FileType;
        end %end setFileType
        
        function oRCF = setTimeStamp_us(oRCF, TimeStamp_us)
            oRCF.m_TimeStamp_us  = TimeStamp_us;
        end %end setZeroDelayTime_us
        
        function oRCF = setFc_Hz(oRCF, Fc_Hz)
            oRCF.m_Fc_Hz = Fc_Hz;
        end %end setFc_Hz
        
        function oRCF = setFs_Hz(oRCF, Fs_Hz)
            oRCF.m_Fs_Hz = Fs_Hz;
        end %end setFs_Hz
        
        function oRCF = setBw_Hz(oRCF, Bw_Hz)
            oRCF.m_Bw_Hz = Bw_Hz;
        end %end setBw_Hz
        
        function oRCF = setNSamples(oRCF, NSamples)
            oRCF.m_NSamples = NSamples;
        end %end setNSamples
        
        function oRCF = setCommentOffset_B(oRCF, CommentOffset_B)
            oRCF.m_CommentOffset_B = CommentOffset_B;
        end %end setCommentOffset_B
        
        function oRCF = setCommentLength(oRCF, CommentLength)
            oRCF.m_CommentLength = CommentLength;
        end %end setCommentLength
        
        function oRCF = setFileSize_B(oRCF, FileSize_B)
            oRCF.m_FileSize_B = FileSize_B;
        end %end setFileSize_B
        
        function oRCF = setComment(oRCF, strComment)
            oRCF.m_strComment = strComment;
        end %end setComment
        
        function oRCF = setFilename(oRCF, strFilename)
            oRCF.m_strFilename = strFilename;
        end %end setFilename
        
        function oRCF = setReferenceData(oRCF, referenceData)
            oRCF.m_fvReferenceData = referenceData;
        end %end setReferenceData
        
        function oRCF = setSurveillanceData(oRCF, surveillanceData)
            oRCF.m_fvSurveillanceData = surveillanceData;
        end %end setSurveillanceData
        
        %% Accessors
        function FileType = getFileType(oRCF)
            FileType = oRCF.m_FileType;
        end %end getFileType
        
        function TimeStamp_us = getTimeStamp_us(oRCF)
            TimeStamp_us = oRCF.m_TimeStamp_us;
        end %end getZeroDelayTime_us
        
        function Fc_Hz = getFc_Hz(oRCF)
            Fc_Hz = oRCF.m_Fc_Hz;
        end %end getFc_Hz
        
        function Fs_Hz = getFs_Hz(oRCF)
            Fs_Hz = oRCF.m_Fs_Hz;
        end %end getFs_Hz
        
        function Bw_Hz = getBw_Hz(oRCF)
            Bw_Hz = oRCF.m_Bw_Hz;
        end %end getBw_Hz
        
        function NSamples = getNSamples(oRCF)
            NSamples = oRCF.m_NSamples;
        end %end getNSamples
        
        function CommentOffset_B = getCommentOffset_B(oRCF)
            CommentOffset_B = oRCF.m_CommentOffset_B;
        end %end getCommentOffset_B
        
        function CommentLength = getCommentLength(oRCF)
            CommentLength = oRCF.m_CommentLength;
        end %end getCommentLength
        
        function FileSize_B = getFileSize_B(oRCF)
            FileSize_B = oRCF.m_FileSize_B;
        end %end getFileSize_B
        
        function strComment = getComment(RCF)
            strComment = RCF.m_strComment;
        end %end getComment
        
        function strFilename = getFilename(oRCF)
            strFilename = oRCF.m_strFilename;
        end %end getFilename
        
        function referenceData = getReferenceData(oRCF)
            referenceData = oRCF.m_fvReferenceData;
        end %end getreferenceData
        
        function surveillanceData = getSurveillanceData(oRCF)
            surveillanceData = oRCF.m_fvSurveillanceData;
        end %end geSurveillanceData
        
        %% Other function:
        
        function strTimeStamp = timeStampToString(oRCF)
            %Work out time since epoc in Matlab format.
            %Matlab works in decimal days so divide but 1000000 to from
            %microseconds to seconds, then by 86400 to get from seconds to
            %days. Note we use 2 hours after epic for ZA timezone.
            
            matlabTime = datenum(1970,1,1,2,0,0) + oRCF.m_TimeStamp_us / 86400 / 1000000;
            strTimeStamp = datestr(matlabTime, 'yyyy-mm-ddTHH.MM.SS');
            strTimeStamp = sprintf('%s.%.6i', strTimeStamp, rem(oRCF.m_TimeStamp_us, 1000000));
        end
        
        %% readHeaderFromFile
        % Reads only RCF header information from file.
        function oRCF = readHeaderFromFile(oRCF, strFilename)
            f = fopen (strFilename, 'rb');
            
            %Check that file is open
            assert(f ~= -1, 'Error:Unable to open file: %s', strFilename);
            
            oRCF.m_strFilename = strFilename;
            
            oRCF.m_FileType = fscanf(f, '%c', 4);
            oRCF.m_TimeStamp_us  = fread(f, 1, 'uint64');
            oRCF.m_Fc_Hz = fread(f, 1, 'uint32');
            oRCF.m_Fs_Hz = fread(f, 1, 'uint32');
            oRCF.m_Bw_Hz = fread(f, 1, 'uint32');
            oRCF.m_NSamples = fread(f, 1, 'uint64');
            oRCF.m_CommentOffset_B = fread(f, 1, 'uint64');
            oRCF.m_CommentLength = fread(f, 1, 'uint32');
            oRCF.m_FileSize_B = fread(f, 1, 'uint64');
            
            fseek(f, oRCF.m_CommentOffset_B, 0);
            oRCF.m_strComment = fscanf(f, '%c', oRCF.m_CommentLength);
            fclose(f);
            
            %fprintf('Finished reading header\n');
        end %end readHeaderFromFile
        
        %% readSamplesFromFile
        % Reads samples from RCF file. Can only be used after calling
        % readHeaderFrom file.
        function oRCF = readSamplesFromFile(oRCF, StartSample, NSamples)
            f = fopen (oRCF.m_strFilename, 'rb');
            
            %Check that file is open
            assert(f ~= -1, 'Error: Unable to open file: %s', oRCF.m_strFilename);
            
            %Check that file is open
            assert(StartSample >= 1, 'Error: Invalid start sample. The first sample in the file is 1');
            
            %Allocate space (will overwrite any previous data)
            oRCF.m_fvReferenceData = complex(zeros(NSamples, 1, 'single'));
            oRCF.m_fvSurveillanceData = complex(zeros(NSamples, 1, 'single'));
            
            fseek(f, cRCF.HEADER_SIZE + (StartSample - 1) * 16, 'bof'); % go to request starting sample in file
            
            %fprintf('Reading sample data...\n');
            
            %Read NSamples
            for ii = 1:NSamples
                oRCF.m_fvReferenceData(ii) = fread(f, 1, 'float32') + 1i * fread(f, 1 ,'float32');
                oRCF.m_fvSurveillanceData(ii) = fread(f, 1, 'float32') + 1i * fread(f, 1, 'float32');
            end
            
            fclose(f);
            
            %Set the number of samples that is now held in this object
            oRCF.m_NSamples = NSamples;
            
            %fprintf('Finished reading sample data.\n');
            
        end %end readSamples
        
        %% readFromFile
        %Opens a file and reads the header information and samples for the
        %specified range:
        function oRCF = readFromFile(oRCF, strFilename, StartSample, NSamples)
            oRCF.readHeaderFromFile(strFilename);
            oRCF.readSamplesFromFile(StartSample, NSamples);
        end %end readFromFile
        
        %% writeToFile
        %Opens a file and writes the header information and all samples:
        function writeToFile(oRCF, strFilename)
            outFile = fopen(strFilename, 'wb');
            
            %Write header conventionally (Empties existing file on open)
            oRCF.m_FileType = 'rcf';
            oRCF.m_FileType = [oRCF.m_FileType 0];
            fwrite(outFile, oRCF.m_FileType, 'char*1');
            fwrite(outFile, oRCF.m_TimeStamp_us, 'int64');
            fwrite(outFile, oRCF.m_Fc_Hz, 'int32');
            fwrite(outFile, oRCF.m_Fs_Hz, 'int32');
            fwrite(outFile, oRCF.m_Bw_Hz, 'int32');
            oRCF.m_NSamples = length(oRCF.m_fvReferenceData); %Update NSamples
            fwrite(outFile, oRCF.m_NSamples, 'uint64');
            oRCF.m_CommentOffset_B = cRCF.HEADER_SIZE + oRCF.m_NSamples * 16; %Update comment offset
            fwrite(outFile, oRCF.m_CommentOffset_B, 'uint64');
            oRCF.m_CommentLength = length(oRCF.m_strComment); %Update comment length
            fwrite(outFile, oRCF.m_CommentLength, 'uint32');
            oRCF.m_FileSize_B = oRCF.m_CommentLength + oRCF.m_CommentOffset_B; %Update filesize
            fwrite(outFile, oRCF.m_FileSize_B, 'uint64');
            
            for ii = 1:oRCF.m_NSamples
                fwrite(outFile, real(oRCF.m_fvReferenceData(ii)), 'float32');
                fwrite(outFile, imag(oRCF.m_fvReferenceData(ii)), 'float32');
                fwrite(outFile, real(oRCF.m_fvSurveillanceData(ii)), 'float32');
                fwrite(outFile, imag(oRCF.m_fvSurveillanceData(ii)), 'float32');
            end
            
            fwrite(outFile, oRCF.m_strComment, 'char*1');
            
            fclose(outFile);
        end %end writeToFile
        
        %% writeHeaderToFile
        %Opens a file and writes the header information and optionally comment:
        %Note this function does not check that values are correct as per
        %sample data as it can be used to edit the header of the file in
        %place. It therefore also doesn't clear the file before writing.
        function writeHeaderToFile(oRCF, strFilename, bWriteComment)
             if nargin == 2
                 bWriteComment = false;
             elseif nargin ~= 3
                 assert('Error incorrect number of arguments')
             end
                 
            outFile = fopen(strFilename, 'rb+');
            fseek(outFile, 0, 'bof');
            
            %Write header conventionally (Empties existing file on open)
            oRCF.m_FileType = 'rcf';
            oRCF.m_FileType = [oRCF.m_FileType 0];
            fwrite(outFile, oRCF.m_FileType, 'char*1');
            fwrite(outFile, oRCF.m_TimeStamp_us, 'int64');
            fwrite(outFile, oRCF.m_Fc_Hz, 'int32');
            fwrite(outFile, oRCF.m_Fs_Hz, 'int32');
            fwrite(outFile, oRCF.m_Bw_Hz, 'int32');
            fwrite(outFile, oRCF.m_NSamples, 'uint64');
            fwrite(outFile, oRCF.m_CommentOffset_B, 'uint64');
            fwrite(outFile, oRCF.m_CommentLength, 'uint32');
            fwrite(outFile, oRCF.m_FileSize_B, 'uint64');
            
            if(bWriteComment == true)
                fseek(outfile, oRCF.m_NSamples * 16, 'cof');
                fwrite(outFile, oRCF.m_strComment, 'char*1');
            end
            
            fclose(outFile);
        end %end writeHeaderToFile
        
        function plotFFT(oRCF, StartSample, NSamples, plotAroundFc)
            %Check that we are not trying to FFT more samples that what are
            %available. And also that the starting sample is not beyond the
            %endof the available range.
            FFTStart = 1;
            FFTSize = 0;
            if(nargin >= 2)
                assert(oRCF.getNSamples() > NSamples, 'Error nSamples exceeds number of samples in RCF');
                assert(nargin >= 3, 'Error starting FFT sample must be paired with NSamples to go in the FFT');
                assert(oRCF.getNSamples() > StartSample, 'Error StartSamp is out of range');
                FFTStart = StartSample;
                FFTSize = NSamples;
            else
                FFTSize = oRCF.getNSamples();
            end
            
            frequencyTicks = (-FFTSize/2:(FFTSize-1)/2)*oRCF.m_Fs_Hz/1e6/FFTSize; %frequency ticks for fft plot (X axis)
            
            %plot about centre frequency else plot at base band
            if nargin ~= 4 || plotAroundFc == 1
                frequencyTicks = frequencyTicks + oRCF.m_Fc_Hz/1e6;
            end
            
            Fref=fftshift(fft(oRCF.m_fvReferenceData(FFTStart:FFTStart + FFTSize - 1)));
            Fsurv=fftshift(fft(oRCF.m_fvSurveillanceData(FFTStart:FFTStart + FFTSize - 1)));
            
            zeroRef = max(abs(Fref)); % 0 dB reference for both plots
            
            fprintf('Finished calculating FFT\n')
            fprintf('Plotting figures...\n')
            
            figure;
            subplot(2,1,1)
            plot(frequencyTicks,20*log10(abs(Fref)/zeroRef),'b')
            xlabel({'Frequency [MHz]'});
            ylabel({'Normalised Power [dB]'});
            title([oRCF.m_strFilename ' - Reference Channel'], 'Interpreter','none');
            
            subplot(2,1,2)
            plot(frequencyTicks,20*log10(abs(Fsurv)/zeroRef),'r')
            xlabel({'Frequency [MHz]'});
            ylabel({'Normalised Power [dB]'});
            title([oRCF.m_strFilename ' - Surveillance Channel'], 'Interpreter','none');
        end %end plotFFT
        
    end %end normal methods
    
    
end %end classdef

