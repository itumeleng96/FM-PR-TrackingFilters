classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,

    end
    
    methods
        function obj = MTT(tracks,confirmationThreshold,deletionThreshold)
            %MTT Construct an instance of this class
            obj.tracks = tracks;
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;

        end
        
        function tracks = assignDetectionToTrack(obj,detections)
            %This method assigns Detections to the nearest Track, else
            disp("Assigning Detections");
            numberOfDetections=size(detections,2);
            if isempty(tracks)
                %create tracks = number of detections
                for i=1:numberOfDetections
                    obj.tracks(i) = track([],[],i,0); 
                end

            end

            %Assign Tracks to Detection using GNN and update filter with
            %new measurements
            %Create a gating threhold for each track
            %Based on detections inside gate assign and 
            %If detection is unassigned,create new track 


            %tracks().updateTrueTrack
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
        
    end
end

