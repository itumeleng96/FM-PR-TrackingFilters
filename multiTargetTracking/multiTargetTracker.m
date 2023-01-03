classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,
        gatingThreshold,       %Radius around the predicted measurement to eliminate other measurements

    end
    
    methods
        function obj = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold)
            %MTT Construct an instance of this class
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;
            obj.gatingThreshold = gatingThreshold;

        end
        
        function obj = assignDetectionToTrack(obj,detections)
            %This method assigns Detections to the nearest Track, else
            numberOfDetections=size(detections,2);
            
            if isempty(obj.tracks)
                disp("Assigning Detections");
                %create tracks = number of detections for the first time

                for i=1:numberOfDetections
                    if i ==1
                        obj.tracks = [track([detections(1,i);detections(2,i)],[;],0,i,0,0)];
                    
                    else
                        obj.tracks(end+1)=track([detections(1,i);detections(2,i)],[;],0,i,0,0);
                    end
                end
                        
            else
                disp("Adding Detections");
                %Assign Tracks to Detection using GNN and update filter with new measurements
                %Get qualifying detections within radius
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
                        obj.tracks(end+1) = track([detections(1,i);detections(2,i)],[;],0,0,0,0);  %Still need a workaround TrackIds
                    end
                end
            end
        end
        function tracks = deleteTracks(obj)
            %Delete Tracks based on deletion Treshold
            for i=1:max(size(obj.tracks))
                if obj.tracks(i).sampleSinceLastUpdate > obj.deletionThreshold
                    obj.tracks(i)=[];
                end
            end
            tracks = obj.tracks;
        end
        
        function tracks = confirmTracks(obj)
            %Confirm Tracks based on confirmation Threshold
            for i=1:max(size(obj.tracks))
                if obj.tracks(i).numberOfUpdates > obj.confirmationThreshold
                    obj.tracks(i).confirmed = 1;
                end
            end

            tracks =obj.tracks;
        end

        function obj = maintainTracks(obj)
            disp("Maintain Tracks");
            obj.deleteTracks();
            obj.confirmTracks();

        end

        function obj = trackingFilter(obj)
            %call tracking filter on all tracks
            disp("Tracking filter");
            numberOfTracks = max(size(obj.tracks));

            for i=1:numberOfTracks
                obj.tracks(i)=obj.tracks(i).predictTrack();
            end

            
        end

        function plotMultiTargetTracking(obj,fs,fd_max,td_max,index,f,RDM)

            figure(f);
            Ndelay = floor(td_max*fs);                                 
            time = 0:1/fs:Ndelay/fs;
            frequency = -fd_max:1:fd_max;
            imagesc(time,frequency,RDM*0);

            text(0,10,"Time:" + index+ "s");
            axis xy;
            colorbar;
            xlabel('Bistatic delay [s]','Fontsize',10);
            ylabel('Doppler frequency [Hz]','Fontsize',10);
            grid on;
            title('Targets centroids and EKF Prediction');
            for i=1:max(size(obj.tracks))
                hold on;
                plot(time((round(obj.tracks(i).predictedTrack(1,:)))),frequency(round(obj.tracks(i).predictedTrack(2,:))),'^-','MarkerFaceColor',	[0 0 0], 'MarkerSize', 7);
                text(time((round(obj.tracks(i).predictedTrack(1,end)))),frequency(round(obj.tracks(i).predictedTrack(2,end))),{sprintf('T%d',i)},'VerticalAlignment','top','HorizontalAlignment','left')
                hold on 
                plot(time((round(obj.tracks(i).trueTrack(1,:)))),frequency(round(obj.tracks(i).trueTrack(2,:))),'y-o','MarkerFaceColor',	[0 0 0], 'MarkerSize', 6);
                legend('Predicted Track','TrueTrack')
            end


        end
    end
    methods(Static)

        function [detectionsInRadius] = pruneDetections(detections,predictedCoordinate,gatingThreshold)
            %For every detection check that it falls within the predicted
            %Coordinates'radius 
           numberOfDetections=size(detections,2);
           detectionsInRadius=[];
           for i=1:numberOfDetections 
                if norm(detections(:,i)-(predictedCoordinate)) < gatingThreshold
                    detectionsInRadius(:,end+1)=detections(:,i);
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


