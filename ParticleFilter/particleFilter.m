classdef particleFilter

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,particles;
    end
    
    methods
        function obj = particleFilter(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,initialCentroid,N)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % u_x : acceleration in the x-direction
            % u_y : acceleration in the y-direction
            % x_std_meas : standard deviation of the measurement in the x-direction
            % y_std_meas : standard deviation of the measurement in the y-direction
            
            
            %Control Input Variables
            obj.U = [u_x;
                     u_y];
    
           %Initial State of the particles
            %Create Gaussian or uniformly distributed particles on Initialization
            obj.X= [0;0;0;0;0;0];
            obj.particles = createGaussianParticles(initialCentroid,[10,10],N);

            %State transition matrix
            obj.dt = dt;


            obj.A = [1,dt,(1/2)*dt^2, 0 ,0 ,0;
                     0, 1, dt,0,0,0;
                     0, 0, 1, 0,0,0;
                     0, 0  0, 1,dt,(1/2)*dt^2;
                     0, 0, 0, 0,1,dt ;
                     0, 0, 0, 0,0,1;];


            %The control input matrix B
            obj.B = [0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;
                     0,0;];

            %Measurement Mapping Matrix 
            obj.H = [1,0,0,0,0,0;
                     0,0,0,1,0,0];

            %Process Noise Covariance
            obj.Q = [(dt^4)/4, (dt^3)/2, (dt^2)/2, 0,0,0;
                     (dt^3)/2, dt^2, dt, 0,0,0;
                     (dt^2)/2, dt, 1, 0,0,0;
                     0, 0, 0, (dt^4)/4, (dt^3)/2, (dt^2)/2;
                     0, 0, 0, (dt^3)/2, dt^2, dt;
                     0, 0, 0, (dt^2)/2, dt,1].*std_acc^2;

            
            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.A,2));

        end
        
        function [X_pred,PF_obj] = predict(obj)
            %Calculate the predicted time state

            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            %obj.X= obj.A * obj.X + obj.B * obj.U;
            
            %particles : 6XN matrix where N is number of particles
            % Calculate the predicted time state for all particles
            X_pred = obj.A * obj.particles';
            obj.particles = X_pred';

            
            %calculate error covariance
            %P= A*P*A' + Q 
            obj.P = eye(size(obj.A,2));
            
            obj.P = (obj.A * obj.P) * obj.A.' + obj.Q;
            
            X_pred = obj.X;
            PF_obj  = obj;
        end


                
     end
 end
        
    


