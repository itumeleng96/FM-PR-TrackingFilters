classdef GN

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y;
    end
    
    methods
        function obj = GN(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial)
        
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
    
            %Initial State
            obj.X= X_initial;
         

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
            obj.Q = [0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0;
                     0, 0, 0, 0,0,0];

            
            %Initial Measurement Noise Covariance
            obj.R = [x_std_meas^2,0;
                     0,y_std_meas^2];

            %Initial covariance matrix -  Identity matrix same shape as A
            obj.P = eye(size(obj.A,2));

        end

        %Function to predict the next state
        function [X_pred,GN_obj] = predict(obj)
            %Calculate the predicted time state
            
            %Update time state
            %x_k = Ax_(k-1) + Bu_(k-1) 
            obj.X= obj.A * obj.X + obj.B * obj.U;
                        
            
            X_pred = obj.X;
            GN_obj  = obj;
        end

        %sum of squared differences
        function [residual] = objectiveFunction(X_pred,z,obj)
            residual = (z-obj.H*X_pred)' * (z-obj.H*X_pred);
        end

        function [J] =jacobian(z,obj)
            %h(p) = [fn,rn]
            % p = [fn,fdotn,fdotdotn;rn,rdotn,rdotdotn]

            %J_h(i, j) = ∂h(p)/∂p_j = ∂h_i/∂p_j
            

            %J = -2 * J_h^T * (z - h(p))

        end
        
        function [GN_obj] = update(obj)
            
        end

     end
 end
        
    


