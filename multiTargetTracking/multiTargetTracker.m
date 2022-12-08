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
            
            disp(numel(obj.tracks));
            if isempty(obj.tracks)
                disp("Assigning Detections");
                %create tracks = number of detections for the first time

                for i=1:numberOfDetections
                    if i ==1
                        obj.tracks = [track([detections(1,i);detections(2,i)],[],0,i,0,0)];
                    
                    else
                        obj.tracks(end+1)=track([detections(1,i);detections(2,i)],[],0,i,0,0);
                    end
                    disp(obj.tracks(i).trueTrack);
                end
                        
            else
                disp("Adding Detections");
                %Assign Tracks to Detection using GNN and update filter with new measurements
                %Get qualifying detections within radius
                numOfTracks = length(obj.tracks);

                for i=1:numOfTracks
                    predictedCoodinate = obj.tracks(i).predictedTrack(end,:);
                    detectionsInRadius = obj.pruneDetections(detections,predictedCoodinate,obj.gatingThreshold);
                    [detection,detections] = obj.globalNearestNeighbour(detectionsInRadius,predictedCoodinate,detections);
                    disp(obj.tracks(i).trueTrack);
                    disp(detection);
                    obj.tracks(i).updateTrueTrack(detection);
                end   
                %If detection is unassigned,create new track 
                
                if ~isempty(detections)
                    numberOfUnassignedDetections = size(detections,2);
                    
                    for i=1:numberOfUnassignedDetections
                        obj.tracks(end+1) = track([detections(1,i);detections(2,i)],[],0,0,0,0);  %Still need a workaround TrackIds
                    end
                end
            end
            disp("trueTrack");
            disp(obj.tracks(1).trueTrack);
            disp(obj.tracks(2).trueTrack);
            
        end
        function tracks = deleteTracks(obj)
            %Delete Tracks based on deletion Treshold
            disp("Delete tracks");
            for i=1:max(size(obj.tracks))
                if obj.tracks(i).sampleSinceLastUpdate > obj.deletionThreshold
                    obj.tracks(i)=[];
                end
            end
            tracks = obj.tracks;
        end
        
        function tracks = confirmTracks(obj)
            %Confirm Tracks based on confirmation Threshold
            disp("Confirm Tracks");
            for i=1:max(size(obj.tracks))
                if obj.tracks(i).numberOfUpdates > obj.confirmationThreshold
                    obj.tracks(i).confirmation = 1;
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
    end
    methods(Static)

        function [detectionsInRadius] = pruneDetections(detections,predictedCoordinate,gatingThreshold)
            disp("Prune");
            %For every detection check that it falls within the predicted
            %Coordinates'radius 
           numberOfDetections=size(detections,2);
           detectionsInRadius=[];
           for i=1:numberOfDetections    
                if norm(detections(:,i)-(predictedCoordinate')) < gatingThreshold
                    detectionsInRadius(:,end+1)=detections(:,i);
                end
           end
           disp(predictedCoordinate');
           disp(detectionsInRadius);
        end

        function [detection,detections] = globalNearestNeighbour(detectionsInRadius,predictedCoordinate,allDetections)
            disp("GNN");
            numberOfDetections=size(detectionsInRadius,2);
            distances=[];

            for i=1:numberOfDetections
                distances(i)=norm(detectionsInRadius(:,i)-predictedCoordinate');
            end

            %Get The index of the detection with min distances and delete from detections
            disp(distances);
            [~,indexOfNeighbour] = min(distances);
            detection = detectionsInRadius(:,indexOfNeighbour); 
            disp(detection);
            %Find index of detection and delete
            if detection
                indexInAllDetections = round(allDetections(1,:))==round(detection(1,1)) & round(allDetections(2,:))==round(detection(2,1));
                allDetections(:,indexInAllDetections) = [];
            end

            detections = allDetections;
        end
    end
end


