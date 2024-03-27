classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,
        gatingThreshold,       %Radius around the predicted measurement to eliminate other measurements
        filterType,            %KalmanFilter:1 , GaussNewton:2
        newtracksCreated,
        M,                      %M-out-of N logic
        N,                      %M-out-of N logic

    end
    
    methods
        function obj = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType)
            %MTT Construct an instance of this class
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;
            obj.gatingThreshold = gatingThreshold;
            obj.filterType = filterType;
            obj.newtracksCreated = 0;
            obj.M =3;
            obj.N =5;
        end
        
        function obj = createNewTracks(obj,detections)
            %This method assigns Detections to the nearest Track, else
            numberOfDetections=size(detections,2);
            if isempty(obj.tracks)
                %disp("Assigning Detections");
                %create tracks = number of detections for the first time

                for i=1:numberOfDetections
                    if i ==1
                        obj.tracks = [track([detections(1,i);detections(2,i)],[;],i,0,0,0,obj.filterType)];
                        obj.newtracksCreated = 1;
                    
                    else
                        obj.tracks(end+1)=track([detections(1,i);detections(2,i)],[;],i,0,0,0,obj.filterType);
                        obj.newtracksCreated = 1;
                    end
                end            
            end
        end

        function obj = updateStage(obj,detections)
            %disp("UpdateStage");
            %Assign Tracks to Detection using GNN and update filter with new measurements
            %Get qualifying detections within radius if not create new tracks

            if ~isempty(obj.tracks)
                numOfTracks = length(obj.tracks);
                for i=1:numOfTracks
                    predictedCoodinate = obj.tracks(i).predictedTrack(:,end);
                    detectionsInRadius = obj.pruneDetections(detections,predictedCoodinate,obj.gatingThreshold);
                    [detection,detections] = obj.globalNearestNeighbour(detectionsInRadius,predictedCoodinate,detections);
                    if(detection)
                        obj.tracks(i)=obj.tracks(i).updateTrueTrack(detection); 
                    end
                end   
                %If detection is unassigned,create new track 
                
                if ~isempty(detections)
                    numberOfUnassignedDetections = size(detections,2);
                    
                    for i=1:numberOfUnassignedDetections
                        obj.tracks(end+1) = track([detections(1,i);detections(2,i)],[;],0,0,0,0,obj.filterType);  %Still need a workaround TrackIds
                    end
                end
            end
            obj.newtracksCreated = 0;
        end

        function tracks = deleteTracks(obj)
            %Delete Tracks based on deletion Treshold
            idx_to_delete =[];
            for i=1:max(size(obj.tracks))
                %if obj.tracks(i).sampleSinceLastUpdate > obj.deletionThreshold
                if obj.tracks(i).seenCountDel > obj.M  %M-out-of-N logic
                   idx_to_delete = [idx_to_delete, i];
                end
            end
            obj.tracks(:,idx_to_delete) = [];
            tracks = obj.tracks;
        end
        
        function tracks = confirmTracks(obj)
            %Confirm Tracks based on confirmation Threshold
            for i=1:max(size(obj.tracks))                
                if (obj.tracks(i).confirmed == 0 && obj.tracks(i).seenCount >= obj.M)
                    obj.tracks(i).confirmed = 1;
                end
            end

            tracks =obj.tracks;
        end

        function obj = maintainTracks(obj)
            %disp("Maintain Tracks");
            obj.tracks =obj.deleteTracks();
            obj.tracks =obj.confirmTracks();

        end

        function obj = predictionStage(obj)
            %disp("Prediction Stage");
            numberOfTracks = max(size(obj.tracks));
            
            for i=1:numberOfTracks
                obj.tracks(i)=obj.tracks(i).predictTrack();
            end

            
        end

        function plotMultiTargetTracking(obj,fs,fd_max,td_max,index,f,RDM)
                figure(f);
                c=3e8;
                Ndelay = floor(td_max*fs);                                 
                time = 0:1/fs:Ndelay/fs;
                range = time *c;
                frequency = -fd_max:1:fd_max;
                imagesc(range,frequency,RDM*0);
                colormap(gca, 'white'); % Set the colormap to 'gray'
    
                text(0,10,"Time:" + index+ "s");
                axis xy;
                xlabel('Bistatic Range [m]','Fontsize',10);
                ylabel('Bistatic Doppler frequency [Hz]','Fontsize',10);
                grid on;
                title('Targets centroids and  Prediction');
                
               
                hold on;
                for i = 1:length(obj.tracks)
                    if obj.tracks(i).confirmed == 0
                        plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor','none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Predicted Track');                        % Plot true track as open circles joined by a line
                        plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color','green' ,'MarkerSize', 4, 'DisplayName', 'Confirmed Track');
                    else
                        plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor', 'none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Predicted Track');                        % Plot true track as open circles joined by a line
                        plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4, 'DisplayName', 'Tentative Track');
                    end

                    % Shift and box the text annotations
                    %textOffset = 15; 
                    %annotationX = obj.tracks(i).predictedTrack(1, end) + textOffset;
                    %annotationY = obj.tracks(i).predictedTrack(2, end) + textOffset;
                    
                    %rectangle('Position', [annotationX-0.5, annotationY-0.5, 1, 10], 'EdgeColor', 'black', 'LineWidth', 1); % Box
                    %text(annotationX, annotationY, ['T',num2str(i)], 'Color', 'black', 'FontSize', 8);

                end
                predicted_marker = plot(nan, nan, '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue', 'MarkerSize', 6);
                tentative_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'green', 'MarkerSize', 4);
                confirmed_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4);
                
                % Create a legend with custom markers and labels
                legend([predicted_marker, tentative_marker, confirmed_marker], 'Predicted Track', 'Tentative Track', 'Confirmed Track', 'Location', 'best');
                hold off;
        end
        function [doppler_error, range_error,doppler_meas,range_meas] =getErrors(obj, i,doppler_error,range_error)
            
            % Calculate the Range and Doppler RMS
            for j = 1:1
                predictedTrack = obj.tracks(j).predictedTrack;
                trueTrack = obj.tracks(j).trueTrack;

                
                range_error(j,1:i) = predictedTrack(1, 1:i);
                doppler_error(j,1:i)= predictedTrack(2, 1:i);

                range_meas(j,1:i) = trueTrack(1, 1:i);
                doppler_meas(j,1:i)= trueTrack(2, 1:i);
                
            end


        end

        function [doppler_ll, range_ll] = calculateLogLikelihood(obj, i,doppler_ll,range_ll)
            
            
            %for j = 1:length(obj.tracks)
            for j = 1:1
                predictedTrack = obj.tracks(j).predictedTrack;
                trueTrack = obj.tracks(j).trueTrack;
                s_matrix = obj.tracks(j).trackingFilterObject.S;
                %----------------------------------------------------------------%
                %--Log-likelihood for Bistatic Range
                %----------------------------------------------------------------%
                range_sample=trueTrack(1,i);
                range_mean =predictedTrack(1,i);
                range_ll(j, i) = logLikelihood(range_mean,s_matrix(1,1),range_sample);

                %----------------------------------------------------------------%
                %--Log-likelihood for Bistatic Doppler
                %----------------------------------------------------------------%
                doppler_sample=trueTrack(2,i);
                doppler_mean =predictedTrack(2,i);
                
                doppler_ll(j, i) = logLikelihood(doppler_mean,s_matrix(2,2),doppler_sample);

            end
    
        end
        function [doppler_ll, range_ll] = plotLogLikelihood(obj, f, f1, i,doppler_ll,range_ll,plotResults)
            
            if(plotResults)
                time = 1:1:i;            
            
                 for j = 1:1
                    predictedTrack = obj.tracks(j).predictedTrack;
                    trueTrack = obj.tracks(j).trueTrack;
                    s_matrix = obj.tracks(j).trackingFilterObject.S;
                    %disp(s_matrix);
                    %----------------------------------------------------------------%
                    %--Log-likelihood for Bistatic Range
                    %----------------------------------------------------------------%
                    range_sample=trueTrack(1,i);
                    range_mean =predictedTrack(1,i);
                    
                    %range_ll(j, i) = normlike([(range_mean),s_matrix(1,1)],range_sample);
                    %range_ll(j, i) = lognlike([range_mean,sqrt(R(1,1))],range_sample);
                    range_ll(j, i) = logLikelihood(range_mean,s_matrix(1,1),range_sample);
    
                    %----------------------------------------------------------------%
                    %--Log-likelihood for Bistatic Doppler
                    %----------------------------------------------------------------%
                    doppler_sample=trueTrack(2,i);
                    doppler_mean =predictedTrack(2,i);
                    
                    %doppler_ll(j, i) = normlike([doppler_mean,s_matrix(2,2)],doppler_sample);
                    %doppler_ll(j, i) = lognlike([doppler_mean,sqrt(R(2, 2))],doppler_sample);
                    doppler_ll(j, i) = logLikelihood(doppler_mean,s_matrix(2,2),doppler_sample);
    
                end
    
                figure(f);
                plot(time, doppler_ll);
                xlabel('Time(s)');
                ylabel('Bistatic Doppler Log-likelihood');
                title('Bistatic Doppler Log-likelihood vs Time');
                %legend('Doppler Log-likelihood', 'Doppler Error'); % Add legend for the two plotted lines
    
                figure(f1);
                plot(time,range_ll);
                xlabel('Time(s)');
                ylabel('Bistatic Range Log-likelihood');
                title('Bistatic Range  Log-likelihood vs Time');
                %legend('Range Log-likelihood', 'Range Error'); % Add legend for the two plotted lines

            end
        end
        
        function [range_ll,doppler_ll] = logLikelihoodMatrix(z,Hx, S )    % Inputs:
            %   - sample: Array of observed  values from the predicted track.
            %   - mean: Mean of the distribution under consideration from the true track.
            %   - R: Covariance value for the parameter (scalar).
            % 
            % Output:
            %   - ll: Log-likelihood for the parameter.
        
           
            % Calculate the log-likelihood for the  parameter
            
            % Log-likelihood for a NORMAL distribution is given by:
            %ll = log(normpdf(sample,mean,S));
            likelihood = log(mvnpdf(z, Hx, S));
            disp(likelihood);
            range_ll =0;
            doppler_ll =0;
        
        end

        function logLikelihood= plotLogLikelihoodMatrix(obj, f, i,logLikelihood,plotResults)
            
            if(plotResults)
                time = 1:1:i;            
            
                
                for j = 1:length(obj.tracks)
                    predictedTrack = obj.tracks(j).predictedTrack;
                    trueTrack = obj.tracks(j).trueTrack;
                    S_matrix = obj.tracks(j).trackingFilterObject.S;
                    S_matrix = S_matrix(1:2,1:2);
                    

                    x=[predictedTrack(1,i);predictedTrack(2,i)];
                    z=[trueTrack(1,i);trueTrack(2,i)];
                    H= [1,0;0,1;];
                    Hx = H*x;
                    
     
                   logLikelihood(j,i) = log(mvnpdf(z, Hx, S_matrix));
    
                end
                figure(f);
                plot(time, logLikelihood);
                xlabel('Time(s)');
                ylabel('Log-likelihood');
                title('Log-likelihood vs Time');
                %legend('Doppler Log-likelihood', 'Doppler Error'); % Add legend for the two plotted lines
    

            end
        end
        function [range_mse,doppler_mse] = calculateMSE(obj,i,range_mse,doppler_mse,true_doppler,true_range)

        
            for j = 1:length(obj.tracks)
                predictedTrack = obj.tracks(j).predictedTrack;
    
                % Calculate squared errors for each track at the current time step
                range_mse(i) = (true_range(i) - predictedTrack(1, i))^2;
                doppler_mse(i) = (true_doppler(i) - predictedTrack(2, i))^2;
    
            end
                

        
        end
    
        function [crlb_doppler,crlb_range] = calculateCRLB(obj, i,crlb_doppler,crlb_range,true_doppler,true_range)
            
            % Calculate the CRLB for each track
            %{
            for j = 1:length(obj.tracks)
                observedTrack = obj.tracks(j).trueTrack;
                sample_variance_doppler = sum((observedTrack(2,i) - true_doppler(i)  ).^2) / 1;
                sample_variance_range = sum((observedTrack(1,i) - true_range(i)  ).^2) / 1;

                variance_mu_doppler = sample_variance_doppler / 1; % Variance of the sample mean
                variance_mu_range = sample_variance_range / 1; % Variance of the sample mean


                Fisher_information_mu_doppler = 1 / variance_mu_doppler;
                Fisher_information_mu_range = 1 / variance_mu_range;

                crlb_doppler(i) = 1/Fisher_information_mu_doppler;
                crlb_range(i) = 1/Fisher_information_mu_range;
                
            end
            %}
            FIM = obj.tracks(1).trackingFilterObject.P;
            crlb = FIM^(-1) ;
            crlb_doppler(i) = crlb(1,1);
            crlb_range(i) = crlb(2,2);

        end
    end

    methods(Static)

        function detectionsInRadius = pruneDetections(detections, predictedCoordinate, gatingThreshold)
            % For every detection check that it falls within the predicted coordinates' radius
            % Rectangular (1-norm)  Gating
         
            numberOfDetections = size(detections, 2);
            detectionsInRadius = [];
            x_std = gatingThreshold(1);
            y_std = gatingThreshold(2);
            k = 3;  %Threshold

            for i = 1:numberOfDetections
                y_dist = abs(detections(2, i) - predictedCoordinate(2));
                x_dist = abs(detections(1, i) - predictedCoordinate(1));
                
                if x_dist<k*x_std && y_dist<k*y_std
                    detectionsInRadius = [detectionsInRadius, detections(:, i)];
                end
            end
            
        end

        function [detection,detections] = globalNearestNeighbour(detectionsInRadius,predictedCoordinate,allDetections)
            numberOfDetections=size(detectionsInRadius,2);
            distances=[];
            measurementCovariance = [100,0;0,0.1];
            %Use Mahalanobis Distance
            for i=1:numberOfDetections
                %distances(i)=norm(detectionsInRadius(:,i)-predictedCoordinate);
                innovation = detectionsInRadius(:,i)-predictedCoordinate;
                %disp("Innovation");
                %disp(innovation);
                mahalanobisDistance = sqrt(innovation' / measurementCovariance*innovation);
                %disp(mahalanobisDistance);

                distances(i)=mahalanobisDistance;
            end

            %Get The index of the detection with min distances and delete from detections
            [~,indexOfNeighbour] = min(distances);
            detection = detectionsInRadius(:,indexOfNeighbour); 
            %Find index of detection and delete
            if detection
                indexInAllDetections = round(allDetections(1,:))==round(detection(1,1)) & round(allDetections(2,:))==round(detection(2,1));
                allDetections(:,indexInAllDetections) = [];
            end

            detections = allDetections;
        end

    end
end