classdef GN

    properties
        dt,U,X,A,B,H,Q,R,P,coeff,measured_x,measured_y,max_iter,tolerance;
    end
    
    methods
        function obj = GN(dt,u_x,u_y,std_acc,x_std_meas,y_std_meas,X_initial,max_iter,tolerance)
        
            %Init funtion
            %Inputs: 
            % dt : smapling time
            % u_x : acceleration in the x-direction
            % u_y : acceleration in the y-direction
            % x_std_meas : standard deviation of the measurement in the x-direction
            % y_std_meas : standard deviation of the measurement in the y-direction
            
            %variables for Stopping criterion
            obj.max_iter=max_iter;
            obj.tolerance = tolerance;
            
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

        function [J] =jacobian(obj,z)
            %h(p) = [fn,rn]
            % p = [fn,fdotn,fdotdotn;rn,rdotn,rdotdotn]

            %J_h(i, j) = ∂h(p)/∂p_j = ∂h_i/∂p_j
            

            %J = -2 * J_h^T * (z - h(p))

        end
        
        function [GN_obj] = update(obj,z)
            
            iteration = 0;
            while true
                %objective function
                residual = objectiveFunction(obj.X,z);

                %Compute the jacobian matrix
                J= jacobian(obj,z);
    
                %Compute delta ,dx =(J^T * J)^-1 * J^T * r
                delta = (J' * J)\ J' * residual;
                
                %Update state parameters
                obj.X = obj.X +delta;
                
                %stopping criterion
                if norm(delta) < obj.tolerance || iteration >=obj.max_iter
                    break;
                end
                
                %increment iterations
                iteration = iteration+1;
               
            end            
            
            GN_obj=obj;
        end

     end
 end
        
    


