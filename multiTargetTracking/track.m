classdef track
    %TRACK 
    %   This class holds all informations about a Track
    %   still need to find a way to make type of filter configurable
    properties
       sampleSinceLastUpdate,
       seenCountDel,
       totalUpdatesDel,
       seenCount,
       totalUpdates,
       trueTrack,
       predictedTrack,
       sMatrix,
       trackingFilterObject,
       trackId,
       confirmed,            %confirmed=1: Confirmed ,confirmed=0: Tentative
       deleted,
       numberOfUpdates,
       newTrack,
       x_initial,
       startTime,
    end
    
    methods
        function obj = track(trueTrack,predictedTrack,trackId,startTime,confirmation,sampleSinceLastUpdate,numberOfUpdates,filterType)
            obj.trueTrack = trueTrack;
            obj.predictedTrack = predictedTrack;
            obj.sMatrix = [];
            obj.sampleSinceLastUpdate = sampleSinceLastUpdate;
            obj.totalUpdates=0;
            obj.totalUpdatesDel=0;
            obj.seenCount=0;
            obj.seenCountDel=0;
            obj.trackId = trackId;
            obj.startTime = startTime;
            obj.confirmed = confirmation;
            obj.deleted = 0;
            obj.numberOfUpdates = numberOfUpdates;
            %Initialize Tracker
            %Create random initial points within observation space
            obj.x_initial = [trueTrack(1,1),trueTrack(2,1)];
            dt=1;                              %Time step between samples(update time)
            U=[0,0];                           %Input values x(Delay) and y(Doppler shift) 
            obj.newTrack =1;
            switch filterType
                case 1
                    disp("Initializing Kalman Filter");
                    std_meas=[500,0.1];                                %Standard Deviation of the measurements in the x and y
                    std_acc=0.09;                                       %Standard Deviation of the process noise
                    KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = KF_object;
                
                case 2
                    disp("Initializing Huber Covariance Scaling Kalman Filter");
                    std_meas=[500,0.1];                                 %Standard Deviation of the measurements in the x and y
                    std_acc=0.09;                                       %Standard Deviation of the process noise
                    HSCKF_object = HCSKF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = HSCKF_object;
                                   
                case 3
                    disp("Initializing Particle Filter");
                    N=10000;                                             %Number of particles
                    std_acc=1;                                           %Standard Deviation of the process noise
                    std_meas=[500,2];                                  %Standard Deviation of the measurements in the x and y
                    PF_object = particleFilter(dt,std_acc,std_meas,[obj.x_initial(1);0;obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = PF_object;
                
                case 4
                    disp("Initializing Unscented Kalman Filter");
                    std_acc=1;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[500,0.1];                              %Standard Deviation of the measurements in the x and y
                    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;]);
                    obj.trackingFilterObject = UKF_object;

                case 5
                    disp("Initializing Huber Covariance Scaling Unscented Kalman Filter");
                    std_acc=0.09;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[500,0.1];                              %Standard Deviation of the measurements in the x and y
                    CSUKF_object = CSUKF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;]);
                    obj.trackingFilterObject = CSUKF_object;

                case 6
                    disp("Initializing Recursive Gauss Newton Filter");
                    std_acc=1e3;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[500,0.1];                               %Standard Deviation of the measurements in the x and y
                    RGNF_object = RGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100);
                    obj.trackingFilterObject = RGNF_object;
    
               case 7
                    disp("Initializing Covariance Scaling Recursive Gauss Newton Filter");
                    std_acc=1;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[500,0.5];                               %Standard Deviation of the measurements in the x and y
                    CSRGNF_object = CSRGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100);
                    obj.trackingFilterObject = CSRGNF_object;
    
                case 8
                    disp("Initializing Filter");
                    std_acc=0.01;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[100,0.1];                               %Standard Deviation of the measurements in the x and y
                    EMP_object = FMP(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;]);
                    obj.trackingFilterObject = EMP_object;
                
                case 9
                    disp("Initializing Filter");
                    std_acc=0.01;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[100,0.1];                               %Standard Deviation of the measurements in the x and y
                    EMP_object = EMP(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;]);
                    obj.trackingFilterObject = EMP_object;

                case 10
                    disp("Initializing Filter");
                    std_acc=0.01;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[100,0.1];                               %Standard Deviation of the measurements in the x and y
                    compositePF_object = compositePF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;]);
                    obj.trackingFilterObject = compositePF_object;
                    
                otherwise
                    dt=1;
                    std_meas=[25,0.001];                              %Standard Deviation of the measurements in the x and y
                    std_acc=[1e-3,1];                                 %Standard Deviation of the acceleration in ms^2
                    KF_object = kalmanFilter(dt,U(1),U(2),std_acc,std_meas(1),std_meas(2),[x_initial(1);0;0;x_initial(2);0;0]);
                    obj.trackingFilterObject = KF_object;
            end
        end
        
        function obj = updateTrueTrack(obj,newTargetObservation)

            %True Track is already updated for the first observation
            if(~obj.newTrack)
                %Insert Observations from Target Observations
                obj.trueTrack(1,end+1) = newTargetObservation(1,1);
                obj.trueTrack(2,end) = newTargetObservation(2,1);
                
            end

            %Update Tracking Filter 
            [~,obj.trackingFilterObject] = update(obj.trackingFilterObject,[newTargetObservation(1,1);newTargetObservation(2,1)]); 

            %Update sampleSinceLastUpdate and number Of Updates
            obj.sampleSinceLastUpdate = 0;
            obj.numberOfUpdates = obj.numberOfUpdates+1;

            %obj.seenCountDel = obj.seenCountDel-1;
            %obj.seenCount = obj.seenCount+1;
            obj.newTrack=0;

        end

        function obj = predictTrack(obj)
            %Predict using Tracking filter and update predictedTrack
            %N=5;
            [X,obj.trackingFilterObject]= predict(obj.trackingFilterObject);
            %Update the predicted track
            obj.predictedTrack(1,end+1)=X(1,1);
            obj.predictedTrack(2,end)=X(3,1);   
            obj.sampleSinceLastUpdate  = obj.sampleSinceLastUpdate+1;

            if(size(obj.trackingFilterObject.S,1)>1)
                obj.sMatrix(1,end+1) = obj.trackingFilterObject.S(1,1);
                obj.sMatrix(2,end) = obj.trackingFilterObject.S(2,2);
            end

            %M Out of N logic 
            %obj.totalUpdatesDel =obj.totalUpdatesDel+1;
            %if obj.totalUpdatesDel>N
            %    obj.seenCountDel=obj.seenCountDel+1;
            %end

            %obj.totalUpdates =obj.totalUpdates+1;
            %if obj.totalUpdates>N
            %    obj.seenCount=obj.seenCount-1;
            %end
        end
 
        function obj = incrementSampleSinceLastUpdate(obj) %call function to keep track of 
            obj.sampleSinceLastUpdate = obj.sampleSinceLastUpdate+1;
        end

    end
end

