classdef track
    %TRACK 
    %   This class holds all informations about a Track
    %   still need to find a way to make type of filter configurable
    properties
       sampleSinceLastUpdate,
       trueTrack,
       predictedTrack,
       trackingFilterObject,
       trackId,
       confirmed,            %confirmed=1: Confirmed ,confirmed=0: Tentative
       numberOfUpdates,
       newTrack,
       x_initial,
    end
    
    methods
        function obj = track(trueTrack,predictedTrack,trackId,confirmation,sampleSinceLastUpdate,numberOfUpdates,filterType)
            obj.trueTrack = trueTrack;
            obj.predictedTrack = predictedTrack;
            obj.sampleSinceLastUpdate = sampleSinceLastUpdate;
            obj.trackId = trackId;
            obj.confirmed = confirmation;
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
                    std_meas=[1000,10];                                %Standard Deviation of the measurements in the x and y
                    std_acc=1;                                       %Standard Deviation of the process noise
                    KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;]);
                    obj.trackingFilterObject = KF_object;
                                   
                case 2
                    disp("Initializing Particle Filter");
                    N=5000;                                          %Number of particles
                    std_acc=1;                                     %Standard Deviation of the process noise
                    std_meas=[1000,5];                                %Standard Deviation of the measurements in the x and y
                    PF_object = particleFilter(dt,std_acc,std_meas,[obj.x_initial(1);obj.x_initial(2);0;],N);
                    obj.trackingFilterObject = PF_object;
                
                case 3
                    disp("Initializing Unscented Kalman Filter");
                    std_acc=1;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[1000,0.1];                              %Standard Deviation of the measurements in the x and y
                    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1),obj.x_initial(2),0;]);
                    obj.trackingFilterObject = UKF_object;

                case 4
                    disp("Initializing Recursive Gauss Newton Filter");
                    std_acc=1;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[1000,0.1];                               %Standard Deviation of the measurements in the x and y
                    RGNF_object = RGNF(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;],100);
                    obj.trackingFilterObject = RGNF_object;
    
                case 5
                    disp("Initializing Fading Polynomial Filter");
                    std_acc=0.01;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[100,0.1];                               %Standard Deviation of the measurements in the x and y
                    EMP_object = FMP(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;]);
                    obj.trackingFilterObject = EMP_object;
                
                case 6
                    disp("Initializing Expanding Polynomial Filter");
                    std_acc=0.01;                                     %Standard Deviation of the acceleration in ms^2
                    std_meas=[100,0.1];                               %Standard Deviation of the measurements in the x and y
                    EMP_object = EMP(dt,std_acc,std_meas(1),std_meas(2),[obj.x_initial(1);obj.x_initial(2);0;]);
                    obj.trackingFilterObject = EMP_object;

                case 7
                    disp("Initializing Composite Polynomial Filter");
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
            %disp("true Track");
            %disp(obj.trueTrack);
            obj.newTrack=0;

        end

        function obj = predictTrack(obj)
            %Predict using Tracking filter and update predictedTrack

            [X,obj.trackingFilterObject]= predict(obj.trackingFilterObject);
            %Update the predicted track
            obj.predictedTrack(1,end+1)=X(1,1);
            obj.predictedTrack(2,end)=X(2,1);   
            %disp("predicted Track");
            %disp(obj.predictedTrack);


            
        end
 
        function obj = incrementSampleSinceLastUpdate(obj) %call function to keep track of 
            obj.sampleSinceLastUpdate = obj.sampleSinceLastUpdate+1;
        end

    end
end

