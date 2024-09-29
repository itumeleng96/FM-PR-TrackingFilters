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
                    std_meas=[4.9038,0.9985];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    std_acc=[0.0048354,0.0991];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = KF_object;

                case 15
                    disp("Initializing Kalman Filter");
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    std_acc=[0.1,0.02];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = KF_object;
                    
                case 16
                    disp("Initializing Kalman Filter");
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    std_acc=[0.05,0.025];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = KF_object;
                
                case 17
                    disp("Initializing Kalman Filter");
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    std_acc=[0.05,0.05];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = KF_object;

                case 2
                    disp("Initializing Covariance Scaling Kalman Filter");
                    std_meas=[1,0.2];                                 %Standard Deviation of the measurements in the x and y
                    std_acc=[0.001,0.02];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    CSKF_object = CSKF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;]);
                    obj.trackingFilterObject = CSKF_object;
                                   
                case 3
                    disp("Initializing Particle Filter");
                    N=2500;                                             %Number of particles
                    std_acc=2;                                           %Standard Deviation of the process noise
                    std_meas=[2,1];                                  %Standard Deviation of the measurements in the x and y
                    PF_object = particleFilter(dt,std_acc,std_meas,[obj.x_initial(1);0;obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = PF_object;
                case 18
                    disp("Initializing Particle Filter");
                    N=5000;                                             %Number of particles
                    std_acc=2;                                           %Standard Deviation of the process noise
                    std_meas=[3,1];                                  %Standard Deviation of the measurements in the x and y
                    PF_object = particleFilter(dt,std_acc,std_meas,[obj.x_initial(1);0;obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = PF_object;
                case 19
                    disp("Initializing Particle Filter");
                    N=7500;                                             %Number of particles
                    std_acc=2;                                           %Standard Deviation of the process noise
                    std_meas=[4,1];                                  %Standard Deviation of the measurements in the x and y
                    PF_object = particleFilter(dt,std_acc,std_meas,[obj.x_initial(1);0;obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = PF_object;
                
                case 20
                    disp("Initializing Particle Filter");
                    N=10000;                                             %Number of particles
                    std_acc=2;                                           %Standard Deviation of the process noise
                    std_meas=[4,1];                                  %Standard Deviation of the measurements in the x and y
                    PF_object = particleFilter(dt,std_acc,std_meas,[obj.x_initial(1);0;obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = PF_object;
                case 4
                    disp("Initializing Covariance Scaling Particle Filter");
                    N=10000;                                             %Number of particles
                    std_acc=2;                                           %Standard Deviation of the process noise
                    std_meas=[2,1];                                  %Standard Deviation of the measurements in the x and y
                    CSPF_object = CSPF(dt,std_acc,std_meas,[obj.x_initial(1);0;obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = CSPF_object;

                case 5
                    disp("Initializing Unscented Kalman Filter");
                    std_acc=[0.09,0.2];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;],0.01,0,0.01); %alpha,Kappa,Beta 
                    obj.trackingFilterObject = UKF_object;

                 case 9
                    disp("Initializing Unscented Kalman Filter");
                    std_acc=[0.09,0.2];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;],0.01,0,2); %alpha,Kappa,Beta 
                    obj.trackingFilterObject = UKF_object;
                
                case 10
                    disp("Initializing Unscented Kalman Filter");
                    std_acc=[0.09,0.2];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;],0.01,0,4); %alpha,Kappa,Beta 
                    obj.trackingFilterObject = UKF_object;
                case 11
                    disp("Initializing Unscented Kalman Filter");
                    std_acc=[0.09,0.2];                             %Standard Deviation of the process noise x(Range) and y(Doppler)
                    std_meas=[2,0.2];                                 %Standard Deviation of the measurements in the x(Range) and y(Doppler)
                    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;],0.01,0,8); %alpha,Kappa,Beta 
                    obj.trackingFilterObject = UKF_object;
               

                case 6
                    disp("Initializing Covariance Scaling Unscented Kalman Filter");
                    std_acc=0.9;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[5,0.5];                              %Standard Deviation of the measurements in the x and y
                    CSUKF_object = CSUKF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),0,obj.x_initial(2),0;]);
                    obj.trackingFilterObject = CSUKF_object;

                case 7
                    disp("Initializing Recursive Gauss Newton Filter");
                    std_acc=[0.09,0.2];                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[2,0.1];                               %Standard Deviation of the measurements in the x and y
                    RGNF_object = RGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100,1);
                    obj.trackingFilterObject = RGNF_object;
    
               case 8
                    disp("Initializing Covariance Scaling Recursive Gauss Newton Filter");
                    std_acc=0.5;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[5,0.1];                               %Standard Deviation of the measurements in the x and y
                    CSRGNF_object = CSRGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100);
                    obj.trackingFilterObject = CSRGNF_object;
                    
                case 12
                    disp("Initializing   Recursive Gauss Newton Filter");
                    std_acc=[0.09,0.2];                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[2,0.1];                               %Standard Deviation of the measurements in the x and y
                    RGNF_object = RGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100,0.25);
                    obj.trackingFilterObject = RGNF_object;
                case 13
                    disp("Initializing   Recursive Gauss Newton Filter");
                    std_acc=[0.09,0.2];                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[2,0.1];                               %Standard Deviation of the measurements in the x and y
                    RGNF_object = RGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100,0.5);
                    obj.trackingFilterObject = RGNF_object;
                case 14
                    disp("Initializing   Recursive Gauss Newton Filter");
                    std_acc=[0.09,0.2];                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[2,0.1];                               %Standard Deviation of the measurements in the x and y
                    RGNF_object = RGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);0;obj.x_initial(2);0;],100,0.75);
                    obj.trackingFilterObject = RGNF_object;

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

            if(size(obj.trackingFilterObject.S,1)>0)
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

