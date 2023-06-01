classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,
        gatingThreshold,       %Radius around the predicted measurement to eliminate other measurements
        filterType,            %KalmanFilter:1 , GaussNewton:2
        newtracksCreated,

    end
    
    methods
        function obj = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType)
            %MTT Construct an instance of this class
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;
            obj.gatingThreshold = gatingThreshold;
            obj.filterType = filterType;
            obj.newtracksCreated = 0;

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

            if ~isempty(obj.tracks) && ~obj.newtracksCreated
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
                if obj.tracks(i).sampleSinceLastUpdate > obj.deletionThreshold
                   idx_to_delete = [idx_to_delete, i];
                end
            end
            obj.tracks(:,idx_to_delete) = [];
            tracks = obj.tracks;
        end
        
        function tracks = confirmTracks(obj)
            %Confirm Tracks based on confirmation Threshold
            for i=1:max(size(obj.tracks))                
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
            %call tracking filter on all tracks
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
            ylabel('Doppler frequency [Hz]','Fontsize',10);
            grid on;
            title('Targets centroids and  Prediction');
            
           
            hold on;
            for i = 1:length(obj.tracks)
                if obj.tracks(i).confirmed == 0
                    % Plot tentative track as a line connected by open triangles
                    plot(obj.tracks(i).predictedTrack(1,:), obj.tracks(i).predictedTrack(2,:), 'b--^', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'DisplayName', 'Predicted Track');
                    % Plot true track as open circles joined by a line
                    plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), 'o:', 'MarkerFaceColor', [1 0.5 0], 'MarkerSize', 4, 'DisplayName', 'Tentative Track');
                else
                    % Plot confirmed track as a line connected by filled triangles
                    plot(obj.tracks(i).predictedTrack(1,:), obj.tracks(i).predictedTrack(2,:), 'b--^', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'DisplayName', 'Predicted Track');
                    plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), 'o:', 'MarkerFaceColor',[0 0 0], 'MarkerSize', 4, 'DisplayName', 'Confirmed');
                end
            end
            predicted_marker = plot(nan, nan, 'b--^', 'MarkerFaceColor', 'none', 'MarkerSize', 8);
            tentative_marker = plot(nan, nan, 'o:', 'MarkerFaceColor', [1 0.5 0], 'MarkerSize', 4);
            confirmed_marker = plot(nan, nan, 'o:', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 8);
        
            % Create a legend with custom markers and labels
            legend([predicted_marker,tentative_marker, confirmed_marker], 'Predicted Track', 'Tentative Track', 'Confirmed Track', 'Location', 'best');
                    
            hold off;
        end
        function [doppler_error, range_error] = plotError(obj, f, f1, plotDoppler_Error, plotRange_Error, i)
            time = 1:1:i;
            doppler_error = zeros(length(obj.tracks), length(time));
            range_error = zeros(length(obj.tracks), length(time));
            
            % Calculate the Range and Doppler RMS
            for j = 1:length(obj.tracks)
                predictedTrack = obj.tracks(j).predictedTrack;
                trueTrack = obj.tracks(j).trueTrack;
                
                trackLength = size(trueTrack, 2);
                
                range_diff = predictedTrack(1, 1:i) - trueTrack(1, 1:i);
                doppler_diff = predictedTrack(2, 1:i) - trueTrack(2, 1:i);
                
                range_error(j, 1:i) = abs(range_diff);
                doppler_error(j, 1:i) = abs(doppler_diff);
            end
        
            if plotDoppler_Error
                figure(f);
                plot(time, doppler_error);
                xlabel('Time(s)');
                ylabel('Doppler Error (Hz)');
                title('Doppler Error vs Time');
                legend('Track 1', 'Track 2', 'Track 3', 'Track 4', 'Track 5', 'Track 6', 'Track 7', 'Track 8', 'Track 9', 'Track 10');
            end
        
            if plotRange_Error
                figure(f1);
                plot(time, range_error);
                xlabel('Time(s)');
                ylabel('Range Error (m)');
                title('Range Error vs Time');
                legend('Track 1', 'Track 2', 'Track 3', 'Track 4', 'Track 5', 'Track 6', 'Track 7', 'Track 8', 'Track 9', 'Track 10');
            end
        end

        function [doppler_ll, range_ll] = calculateLogLikelihood(obj, f, f1, i)
            time = 1:1:i;

            doppler_ll = zeros(length(obj.tracks), length(time));
            range_ll = zeros(length(obj.tracks), length(time));
        
            % Calculate the Range and Doppler log-likelihood
            for j = 1:length(obj.tracks)
                predictedTrack = obj.tracks(j).predictedTrack;
                trueTrack = obj.tracks(j).trueTrack;

                trackLength = size(trueTrack, 2);
                
                range_error = predictedTrack(1, 1:i) - trueTrack(1, 1:i);
                doppler_error = predictedTrack(2, 1:i) - trueTrack(2, 1:i);
        
                % Calculate the log-likelihood for Range
                value_1 = -0.5 * log(2 * pi) - 0.5 * log(var(range_error)) - 0.5 * (range_error.^2) / var(range_error);
                range_ll(j,:) = abs(value_1);

        
                % Calculate the log-likelihood for Doppler
                value_2 = -0.5 * log(2 * pi) - 0.5 * log(var(doppler_error)) - 0.5 * (doppler_error.^2) / var(doppler_error);
                doppler_ll(j,:) = abs(value_2);
            
                
            end
            figure(f);
            plot(time, doppler_ll);
            xlabel('Time(s)');
            ylabel('Bistatic Doppler Log-likelihood(Hz)');
            title('Bistatic Doppler Log-likelihood vs Time');
            
            figure(f1);
            plot(time, range_ll);
            xlabel('Time(s)');
            ylabel('Bistatic Range Log-likelihood (m)');
            title('Bistatic Range  Log-likelihood vs Time');
        end
    end
    methods(Static)

        function detectionsInRadius = pruneDetections(detections, predictedCoordinate, gatingThreshold)
            % For every detection check that it falls within the predicted coordinates' radius
            numberOfDetections = size(detections, 2);
            detectionsInRadius = [];
            for i = 1:numberOfDetections
                if abs(detections(1, i) - predictedCoordinate(1)) < gatingThreshold(1) && ...
                   abs(detections(2, i) - predictedCoordinate(2)) < gatingThreshold(2)
                    detectionsInRadius = [detectionsInRadius, detections(:, i)];
                end
            end
        end

        function [detection,detections] = globalNearestNeighbour(detectionsInRadius,predictedCoordinate,allDetections)
            numberOfDetections=size(detectionsInRadius,2);
            distances=[];

            for i=1:numberOfDetections
                distances(i)=norm(detectionsInRadius(:,i)-predictedCoordinate);
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


