classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,
        gatingThreshold,       %Radius around the predicted measurement to eliminate other measurements
        filterType,            %KalmanFilter:1 , GaussNewton:2
        newtracksCreated,
        idCounter,
        %M,                      %M-out-of N logic
        %N,                      %M-out-of N logic

    end
    
    methods
        function obj = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType)
            %MTT Construct an instance of this class
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;
            obj.gatingThreshold = gatingThreshold;
            obj.filterType = filterType;
            obj.newtracksCreated = 0;
            obj.idCounter =0;
            %obj.M =3;
            %obj.N =5;
        end
        
        function obj = createNewTracks(obj,detections,index)
            %This method assigns Detections to the nearest Track, else
            numberOfDetections=size(detections,2);
            if isempty(obj.tracks)
                %disp("Assigning Detections");
                %create tracks = number of detections for the first time

                for i=1:numberOfDetections
                    if i ==1
                        obj.idCounter =obj.idCounter+1;
                        obj.tracks = [track([detections(1,i);detections(2,i)],[;],obj.idCounter,index,0,0,0,obj.filterType)];
                        obj.newtracksCreated = 1;
                    
                    else
                        obj.idCounter =obj.idCounter+1;
                        obj.tracks(end+1)=track([detections(1,i);detections(2,i)],[;],obj.idCounter,index,0,0,0,obj.filterType);
                        obj.newtracksCreated = 1;
                    end
                end            
            end
        end

        function obj = updateStage(obj,detections,index)
            %disp("UpdateStage");
            %Assign Tracks to Detection using GNN and update filter with new measurements
            %Get qualifying detections within radius if not create new tracks

            if ~isempty(obj.tracks)
                numOfTracks = length(obj.tracks);
                for i=1:numOfTracks
                    if(~obj.tracks(i).deleted)
                        predictedCoodinate = obj.tracks(i).predictedTrack(:,end);
                        detectionsInRadius = obj.pruneDetections(detections,predictedCoodinate,obj.gatingThreshold,obj.tracks(i).trackingFilterObject.S);
                        [detection,detections] = obj.globalNearestNeighbour(detectionsInRadius,predictedCoodinate,detections,obj.tracks(i).trackingFilterObject.S);
                        if(detection)
                            obj.tracks(i)=obj.tracks(i).updateTrueTrack(detection); 
                        end
                    end
                end   
                %If detection is unassigned,create new track 
                
                if ~isempty(detections)
                    numberOfUnassignedDetections = size(detections,2);
                    
                    for i=1:numberOfUnassignedDetections
                        obj.idCounter = obj.idCounter+1;
                        obj.tracks(end+1) = track([detections(1,i);detections(2,i)],[;],obj.idCounter,index,0,0,0,obj.filterType); 
                    end
                end
            end
            obj.newtracksCreated = 0;
        end

        function tracks = deleteTracks(obj)
            %Delete Tracks based on deletion Treshold
            %idx_to_delete =[];
            for i=1:max(size(obj.tracks))
                %if obj.tracks(i).seenCountDel > obj.M
                if obj.tracks(i).sampleSinceLastUpdate > obj.deletionThreshold
                   %idx_to_delete = [idx_to_delete, i];
                   obj.tracks(i).deleted =1;
                end
            end
            tracks = obj.tracks;
        end
        
        function tracks = confirmTracks(obj)
            %Confirm Tracks based on confirmation Threshold
            for i=1:max(size(obj.tracks))                
                %if (obj.tracks(i).confirmed == 0 && obj.tracks(i).seenCount >= obj.M)
                if (obj.tracks(i).confirmed == 0 && obj.tracks(i).numberOfUpdates > obj.confirmationThreshold)
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
                if(~obj.tracks(i).deleted)
                    obj.tracks(i)=obj.tracks(i).predictTrack();
                end
            end

            
        end

        function obj = plotMultiTargetTrackingGT(obj, fs, fd_max, td_max, index, f, RDM, rangeTrueData, dopplerTrueData)
            figure(f);
            c = 3e8;
            Ndelay = floor(td_max * fs);                                 
            time = 0:1/fs:Ndelay/fs;
            range = time * c;
            range = range / 1000;
            frequency = -fd_max:1:fd_max;
            imagesc(range, frequency, RDM * 0);
            colormap(gca, 'white');
            
            text(0, 10, "Time:" + index + "s");
            axis xy;
            xlabel('Bistatic Range [km]', 'Fontsize', 18);
            ylabel('Bistatic Doppler frequency [Hz]', 'Fontsize', 18);
            xlim([0 80]);
            ylim([-200 200]);
            grid on;
            title('Targets centroids and Prediction');
            hold on;
            % Plot ground truth data with a distinct star marker for visibility
            % Plot the ground truth track
            gt_marker = plot(rangeTrueData, dopplerTrueData, '-', 'Color', 'black', 'MarkerSize', 4, 'DisplayName', 'Ground Truth', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black');
            
            hold on;
            
            for i = 1:length(obj.tracks)
                if(obj.tracks(i).deleted == 0)
                    trackID = obj.tracks(i).trackId;
                    % Display the track ID as text near the predicted track position
                    text(obj.tracks(i).predictedTrack(1, end) + 1000, obj.tracks(i).predictedTrack(2, end) - 5, ...
                         num2str(trackID), 'Color', 'red', 'FontSize', 8);
                    
                    if obj.tracks(i).confirmed == 0
                        plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor','none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Tracking Filter Prediction ');                        % Plot true track as open circles joined by a line
                        plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color','green' ,'MarkerSize', 4, 'DisplayName', 'Measurement');
                    else
                        plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor', 'none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Tracking Filter Prediction');                        % Plot true track as open circles joined by a line
                        plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4, 'DisplayName', 'Tentative');
                    end
                end
            end
            
            % Create dummy plots for legend clarity
            predicted_marker = plot(nan, nan, '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue', 'MarkerSize', 6);
            tentative_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'green', 'MarkerSize', 4);
            confirmed_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4);
               
            hold off;
            legend([predicted_marker, tentative_marker, confirmed_marker, gt_marker], 'Prediction', 'Tentative', 'Measurement', 'Ground Truth', 'Location', 'best');
            drawnow
        end


        function obj = plotMultiTargetTracking(obj,fs,fd_max,td_max,index,f,RDM)
                figure(f);
                c=3e8;
                Ndelay = floor(td_max*fs);                                 
                time = 0:1/fs:Ndelay/fs;
                range = time *c;
                range = range/1000;
                frequency = -fd_max:1:fd_max;
                imagesc(range,frequency,RDM*0);
                colormap(gca, 'white'); % Set the colormap to 'gray'
    
                text(0,10,"Time:" + index+ "s");
                axis xy;
                xlabel('Bistatic Range [km]','Fontsize',18);
                ylabel('Bistatic Doppler frequency [Hz]','Fontsize',18);
                grid on;
                title('Targets centroids and  Prediction');
                
               
                hold on;
                for i = 1:length(obj.tracks)
                    if(obj.tracks(i).deleted ==0)
                        trackID = obj.tracks(i).trackId;
                        text(obj.tracks(i).predictedTrack(1,end)+1000, obj.tracks(i).predictedTrack(2,end)-5, num2str(trackID), 'Color', 'red', 'FontSize', 8);
                        if obj.tracks(i).confirmed == 0
                            plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor','none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Predicted Track');                        % Plot true track as open circles joined by a line
                            plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color','green' ,'MarkerSize', 4, 'DisplayName', 'Confirmed Track');
                        else
                            plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor', 'none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Predicted Track');                        % Plot true track as open circles joined by a line
                            plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4, 'DisplayName', 'Tentative Track');
                        end
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
                legend([predicted_marker, tentative_marker, confirmed_marker], 'Prediction', 'Tentative', 'Measurement Track', 'Location', 'best');
                %hold on;
        end

        function plotErrors(obj,f,f2,index,rangeGroundTruth,dopplerGroundTruth)
            % Plot Range
            figure(f);
            hold on;
            title('Error comparison Range(m)');
            xlabel('Bistatic Range [m]','Fontsize',10);
            ylabel('Bistatic Doppler [Hz]','Fontsize',10);
            plot(rangeGroundTruth(1:index),dopplerGroundTruth(1:index),'-*');
            
            hold on;
            for i = 1:length(obj.tracks)
                % Check if the track has sufficient length
                if  obj.tracks(i).confirmed==1
                    % Plot the track
                    plot(obj.tracks(i).predictedTrack(1, :),obj.tracks(i).predictedTrack(2, :),'-^');
                    plot(obj.tracks(i).trueTrack(1, :),obj.tracks(i).trueTrack(2, :),'-o');

                end
            end
                ground__marker = plot(nan, nan, '-*');
                predicted_marker = plot(nan, nan, '-^');
                measured_marker = plot(nan, nan, '-o');
                legend([ground__marker, predicted_marker, measured_marker], 'Ground Truth', 'Predicted Track', 'Measured Track', 'Location', 'best');
                grid on;
            hold off;
        end

        function plotErrorTime(obj, f, f2, index,rangeGroundTruth,dopplerGroundTruth)
            % Plot Range
            figure(f);
            hold on;
            title('Error comparison Range(m)');
            xlabel('Time(s)','Fontsize',10);
            ylabel('Bistatic Range [m]','Fontsize',10);
            
            plot(rangeGroundTruth(1:index),'-*');
            for i = 1:length(obj.tracks)
                % Check if the track has sufficient length
                if size(obj.tracks(i).predictedTrack, 2) >= index && obj.tracks(i).confirmed==1
                    % Plot the track
                    plot(obj.tracks(i).predictedTrack(1, :),'-^');
                    plot(obj.tracks(i).trueTrack(1, :),'-');

                else
                    if size(obj.tracks(i).predictedTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).predictedTrack, 2))+1;
                        endIndex = min(index, size(obj.tracks(i).predictedTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).predictedTrack(1, :),'-^');
                    end
                    if size(obj.tracks(i).trueTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).trueTrack, 2))+1;
                        endIndex = min(index, size(obj.tracks(i).trueTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).trueTrack(1, :),'-');
                    end

                end
            end
            hold off;
        
            % Plot Doppler
            figure(f2);
            hold on;
            title('Error comparison Doppler (Hz)');
            xlabel('Time(s)','Fontsize',10);
            ylabel('Bistatic Doppler [Hz]','Fontsize',10);
            plot(dopplerGroundTruth(1:index),'-*');
            for i = 1:length(obj.tracks)
                % Check if the track has sufficient length
                if size(obj.tracks(i).predictedTrack, 2) >= index && obj.tracks(i).confirmed==1
                    plot(obj.tracks(i).predictedTrack(2, :),'-^');
                    plot(obj.tracks(i).trueTrack(2, :),'-');
                else
                    if size(obj.tracks(i).predictedTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).predictedTrack, 2) + 1);
                        endIndex = min(index, size(obj.tracks(i).predictedTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).predictedTrack(2, :),'-^');
                    end
                   if size(obj.tracks(i).trueTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).trueTrack, 2) + 1);
                        endIndex = min(index, size(obj.tracks(i).trueTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).trueTrack(2, :),'-');
                    end 
                end
            end
            hold off;
        end

        function [doppler_error, range_error,doppler_meas,range_meas] =getErrors(obj, i,doppler_error,range_error)
            
            % Calculate the Range and Doppler RMS
            predictedTrack = obj.tracks(1).predictedTrack;
            trueTrack = obj.tracks(1).trueTrack;
            
            range_error(1,1:i) = predictedTrack(1, 1:i);
            doppler_error(1,1:i)= predictedTrack(2, 1:i);

            range_meas(1,1:i) = trueTrack(1, 1:i);
            doppler_meas(1,1:i)= trueTrack(2, 1:i);
                
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

        function [doppler_ll, range_ll,t] = calculateLogLikelihoodGroundTruth(obj,trackId,doppler_ll,range_ll,dopplerGroundTruth,rangeGroundTruth,simTime)
            
             indexOfTrack = -1;
             Tracks = obj.tracks;
             for j = 1:length(Tracks)
                %look for track with the specified Id
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end
             t=[];
             if(indexOfTrack>0)

                 sMatrix = Tracks(indexOfTrack).sMatrix;
                 %sMatrix(1,end+1) = Tracks(indexOfTrack).trackingFilterObject.S(1,1);
                 %sMatrix(2,end) = Tracks(indexOfTrack).trackingFilterObject.S(2,2);
                 

                 predictedTrack = Tracks(indexOfTrack).predictedTrack;

                 range_ll = obj.logLikelihoodMatrix(predictedTrack(1,:),rangeGroundTruth(1,:),sMatrix(1,:),Tracks(indexOfTrack).deleted);
                 doppler_ll = obj.logLikelihoodMatrix(predictedTrack(2,:),dopplerGroundTruth(1,:),sMatrix(2,:),Tracks(indexOfTrack).deleted);
                 %Plot against time 
                 startTime = Tracks(indexOfTrack).startTime;
                 endTime = length(doppler_ll)+startTime -1;
                 
                 t = startTime:1:endTime;
                 
             end
        end

        function [predictedTrack] = getTrack(obj,trackId)
            
             indexOfTrack = -1;
             Tracks = obj.tracks;
             predictedTrack =[];
             for j = 1:length(Tracks)
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end

             if(indexOfTrack>0)
                 predictedTrack = Tracks(indexOfTrack).predictedTrack;
             end
        end

        function [measuredTrack] = getMeasuredTrack(obj,trackId)
            
             indexOfTrack = -1;
             Tracks = obj.tracks;
             measuredTrack =[];
             for j = 1:length(Tracks)
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end

             if(indexOfTrack>0)
                 measuredTrack = Tracks(indexOfTrack).trueTrack;
             end
        end

        function [doppler_ll, range_ll] = plotLogLikelihood(obj, f, f1,trackId,dopplerGroundTruth,rangeGroundTruth,simTime)
                        
             indexOfTrack = -1;
             Tracks = obj.tracks;
             for j = 1:length(Tracks)
                %look for track with the specified Id
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end

             if(indexOfTrack>0)

                 sMatrix = Tracks(indexOfTrack).sMatrix;
                 sMatrix(1,end+1) = Tracks(indexOfTrack).trackingFilterObject.S(1,1);
                 sMatrix(2,end) = Tracks(indexOfTrack).trackingFilterObject.S(2,2);

                 predictedTrack = Tracks(indexOfTrack).predictedTrack;
                 range_ll = obj.logLikelihoodMatrix(predictedTrack(1,:),rangeGroundTruth(1,:),sMatrix(1,:),Tracks(indexOfTrack).deleted);
                 doppler_ll = obj.logLikelihoodMatrix(predictedTrack(2,:),dopplerGroundTruth(1,:),sMatrix(2,:),Tracks(indexOfTrack).deleted);

                 %Plot against time 
                 startTime = Tracks(indexOfTrack).startTime;
                 endTime = length(doppler_ll)+startTime-1;

                 figure(f);
                 hold on;
                 t = startTime:1:endTime;

                 plot(t,doppler_ll);
                
                 xlabel('Time(s)');
                 ylabel('Bistatic Doppler Log-likelihood');
                 title(['Bistatic Doppler Log-likelihood vs Time for Track ID: ', num2str(trackId)]);
                 %legend('Doppler Log-likelihood', 'Doppler Error'); % Add legend for the two plotted lines
    
                 figure(f1);
                 plot(t,range_ll);
                 xlabel('Time(s)');
                 ylabel('Bistatic Range Log-likelihood');
                 title(['Bistatic Range  Log-likelihood vs Time for Track ID:',num2str(trackId)]);
                 %legend('Range Log-likelihood', 'Range Error'); % Add legend for the two plotted lines
                 %}
             end
        end
        
        function [doppler_ll, range_ll] = plotLogLikelihoodSingle(obj, f, f1, i,doppler_ll,range_ll,dopplerTrueData,rangeTrueData,plotResults)
            
            if(plotResults)
                time = 1:1:i;            
            
                 for j = 1:1
                    predictedTrack = obj.tracks(j).predictedTrack;
                    %trueTrack = obj.tracks(j).trueTrack;
                    s_matrix = obj.tracks(j).trackingFilterObject.S;
                    %disp(s_matrix);
                    %----------------------------------------------------------------%
                    %--Log-likelihood for Bistatic Range
                    %----------------------------------------------------------------%
                    range_sample=rangeTrueData(i);
                    range_mean =predictedTrack(1,i);
                    
                    %range_ll(j, i) = normlike([(range_mean),s_matrix(1,1)],range_sample);
                    %range_ll(j, i) = lognlike([range_mean,sqrt(R(1,1))],range_sample);
                    range_ll(j, i) = logLikelihood(range_mean,s_matrix(1,1),range_sample);
    
                    %----------------------------------------------------------------%
                    %--Log-likelihood for Bistatic Doppler
                    %----------------------------------------------------------------%
                    doppler_sample=dopplerTrueData(i);
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

        function [range_mse,doppler_mse] = calculateMSE(obj,i,range_mse,doppler_mse,true_doppler,true_range)

        
            for j = 1:length(obj.tracks)
                predictedTrack = obj.tracks(j).predictedTrack;
    
                % Calculate squared errors for each track at the current time step
                range_mse(i) = (true_range(i) - predictedTrack(1, i))^2;
                doppler_mse(i) = (true_doppler(i) - predictedTrack(2, i))^2;
    
            end
                

        
        end
    
    end

    methods(Static)

        function detectionsInRadius = pruneDetections(detections, predictedCoordinate, gatingThreshold, Smatrix)
            % Inputs:
            % detections - 2D matrix of detection coordinates, each column is a detection
            % predictedCoordinate - 2D vector of the predicted measurement coordinate
            % gatingThreshold - scalar, the gating threshold G in Mahalanobis distance
            % Smatrix - Covariance matrix for the detection residual
        
            % Number of detections
            numberOfDetections = size(detections, 2);
            
            % Preallocate memory for the valid detections matrix
            detectionsInRadius = zeros(size(detections));
            count = 0;  % Counter for valid detections
            
            % Loop through each detection
            for i = 1:numberOfDetections
                % Compute the residual vector
                residual = detections(:, i) - predictedCoordinate;
                
                % Compute the Mahalanobis distance using A\b instead of inv(A)*b
                mahalanobisDistSquared = residual' * (Smatrix \ residual);
                % Check if the Mahalanobis distance is within the gating threshold
                if mahalanobisDistSquared <= gatingThreshold
                    count = count + 1;
                    detectionsInRadius(:, count) = detections(:, i);
                end
            end
            
            % Remove unused preallocated columns
            detectionsInRadius = detectionsInRadius(:, 1:count);
        end
        

        function [nearestDetection, remainingDetections] = globalNearestNeighbour(detectionsInRadius, predictedCoordinate, allDetections, measurementCovariance)
            % Inputs:
            % detectionsInRadius - Detected points within the gating threshold
            % predictedCoordinate - The predicted state of the track
            % allDetections - All available detections
            % measurementCovariance - Covariance matrix for Mahalanobis distance calculation (sMatrix)
        
            % Get the number of detections in the radius
            numDetections = size(detectionsInRadius, 2);
            
            % Preallocate distances array to avoid size changes during loop
            distances = zeros(1, numDetections);
            
            % Loop through detections to calculate Mahalanobis distance
            for i = 1:numDetections
                % Compute innovation (residual between detection and predicted coordinate)
                innovation = detectionsInRadius(:, i) - predictedCoordinate;
                
                % Calculate the Mahalanobis distance
                mahalanobisDistance = sqrt(innovation' / measurementCovariance * innovation);
                
                % Store the Mahalanobis distance
                distances(i) = mahalanobisDistance;
            end
        
            % Find the detection with the minimum Mahalanobis distance
            [~, indexOfNearestNeighbour] = min(distances);
            nearestDetection = detectionsInRadius(:, indexOfNearestNeighbour);
        
            if nearestDetection 
                % Remove the assigned detection from the list of all detections
                % Use a tolerance-based approach to avoid rounding errors
                tolerance = 1e-5;  % Tolerance for matching floating-point numbers
                indexInAllDetections = abs(allDetections(1, :) - nearestDetection(1)) < tolerance & ...
                                       abs(allDetections(2, :) - nearestDetection(2)) < tolerance;
                
                allDetections(:, indexInAllDetections) = [];
            end

            % Remove the matched detection from allDetections
            remainingDetections = allDetections;
        end

        

         function [logLikelihoodValues] = logLikelihoodMatrix(samples,groundTruthValues,s_matrix,deleted)
            %Loop through the length of the track 

            logLikelihoodValues =[];
            if(size(samples,2)== size(groundTruthValues,2))
                for j=1:size(s_matrix,2)
                    logLikelihoodValues = [logLikelihoodValues,logLikelihood(groundTruthValues(j),s_matrix(j),samples(j))];
                end
            end

            %If the predicted values of the track are not the same length as groundTruthValues
            if(size(samples,2) < size(groundTruthValues,2) && deleted)
                %if deleted it means the missing values are at the end
                for j=1:size(s_matrix,2)
                    logLikelihoodValues = [logLikelihoodValues,logLikelihood(groundTruthValues(j),s_matrix(j),samples(j))];
                end
            end

            if(size(samples,2)< size(groundTruthValues,2) && ~deleted)
                %if deleted it means the missing values are at the
                %beginning
                for j=1:size(s_matrix,2)
                    %%Fix time ,so the plot does not start from time=0
                    diff = abs(size(samples,2)-size(groundTruthValues,2));
                    logLikelihoodValues = [logLikelihoodValues,logLikelihood(groundTruthValues(diff+j),s_matrix(j),samples(j))];
                end
            end


         end

    end

end