classdef track
    %TRACK 
    %   This class holds all informations about a Track
    
    properties
       sampleSinceLastUpdate,
       trueTrack,
       predictedTrack,
       trackingFilterObject,
       trackId,
       confirmed,
       numberOfUpdates,
    end
    
    methods
        function obj = track(trueTrack,predictedTrack,trackId,sampleSinceLastUpdate)
            obj.trueTrack = trueTrack;
            obj.predictedTrack = predictedTrack;
            obj.sampleSinceLastUpdate = sampleSinceLastUpdate;
            obj.trackId = trackId;
            obj.confirmed = confrimation;
            obj.numberOfUpdates = numberOfUpdates;

            %Initialize Tracker 
            EKF_object = EKF(dt, 0.1, 0.1, 1, 0.01,0.01,[14;0;0;310;0;0]);
            obj.trackingFilterObject = EKF_object;

        end
        
        function obj = updateTrueTrack(obj,newTargetObservation)
            %Insert Observations from Target Observations
            index = size(obj.trueTrack);
            obj.trueTrack(1,index+1) = newTargetObservation(1,1);
            obj.trueTrack(2,index+1) = newTargetObservation(2,1);
             
            %Update Tracking Filter 
            [~,obj.trackingFilterObject] = update(obj.trackingFilterObject,[newTargetObservation(1,1);newTargetObservation(2,1)]); 

            %Update sampleSinceLastUpdate and number Of Updates
            obj.sampleSinceLastUpdate = 0;
            obj.numberOfUpdates = obj.numberOfUpdates+1;

        end

        function obj = predictTrack(obj)
            %Predict using Tracking filter and update predictedTrack
            [X,obj.trackingFilterObject]= predict(obj.trackingFilterObject);

            %Update the predicted track
            index = size(obj.predictedTrack,1)+1;
            obj.predictedTrack(1,index+1)=X(1,1);
            obj.predictedTrack(2,index+1)=X(4,1);
          
        end

        function obj = incrementSampleSinceLastUpdate(obj)
            obj.sampleSinceLastUpdate = obj.sampleSinceLastUpdate+1;
        end

    end
end

