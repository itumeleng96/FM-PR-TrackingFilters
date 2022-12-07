classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,
        gatingThreshold,       %Radius around the predicted measurement to eliminate other measurements

    end
    
    methods
        function obj = MTT(tracks,confirmationThreshold,deletionThreshold,gatingThreshold)
            %MTT Construct an instance of this class
            obj.tracks = tracks;
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;
            obj.gatingThreshold = gatingThreshold;

        end
        
        function tracks = assignDetectionToTrack(obj,detections)
            %This method assigns Detections to the nearest Track, else
            disp("Assigning Detections");
            numberOfDetections=size(detections,2);
            if isempty(obj.tracks)
                %create tracks = number of detections for the first time
                for i=1:numberOfDetections
                    obj.tracks(i) = track([],[],i,0); 
                end

            else
                %Assign Tracks to Detection using GNN and update filter with new measurements
                %Get qualifying detections within radius
                numOfTracks = size(tracks,1);

                for i=1:numOfTracks
                    predictedCoodinate = obj.tracks(i).predictedTrack(end,:);
                    detectionsInRadius = pruneDetections(detections,predictedCoodinate,obj.gatingThreshold);
                    [detection,detections] = globalNearestNeighbour(detectionsInRadius,predictedCoodinate,detections);
                    obj.tracks(i).updateTrueTrack(detection);
                end   
                %If detection is unassigned,create new track 

                if ~isempty(detections)
                    numberOfUnassignedDetections = size(detections,2);
                    
                    for i=1:numberOfUnassignedDetections
                        obj.tracks(end+1,:) = track([],[],0,0);  %Still need a workaround TrackIds
                    end
                end
            end
            
            tracks=obj.tracks;

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

        function tracks = maintainTracks(obj)
            disp("Maintain Tracks");
            tracks = obj.deleteTracks();
            tracks = obj.confirmTracks();

        end

        function tracks = trackingFilter(obj)
            %call tracking filter on all tracks
            disp("Tracking filter");
            numberOfTracks = max(size(obj.tracks));

            for i=1:numberOfTracks
                obj.tracks(i)=obj.tracks(i).predictTrack();
            end
            tracks = obj.tracks;
        end

        function detectionsInRadius = pruneDetections(detections,predictedCoordinate,gatingThreshold)
            %For every detection check that it falls within the predicted
            %Coordinates'radius 
           numberOfDetections=size(detections,2);
           for i=1:numberOfDetections
                if norm(detections(1,:)-(predictedCoordinate)) < gatingThreshold
                    detectionsInRadius(end+1,:)=detections(1,:);
                end
           end
        end

        function [detection,detections] = globalNearestNeighbour(detectionsInRadius,predictedCoordinate,allDetections)
            numberOfDetections=size(detectionsInRadius,2);
            distances;
            for i=1:numberOfDetections
                distances(i)=norm(detectionsInRadius(i,:)-predictedCoordinate);
            end

            %Get The index of the detection with min distances and delete from detections
            [indexOfNeighbour,~] = min(distances);
            detection = detectionsInRadius(indexOfNeighbour,:); 

            %Find index of detection and delete
            indexInAllDetections = allDetections(1,:)==detection(1,1) & allDetections(2,:)==detection(2,1);
            detections(indexInAllDetections,:) = [];
        end
    end
end


