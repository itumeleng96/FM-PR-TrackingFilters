classdef CSRGNF

    properties
        dt,U,X,F,A,B,H,Q,R,P,S,coeff,measured_x,measured_y,max_iter,wk,count,updater,update1,epsDoppler,epsRange;
    end
    
    methods
        function obj = CSRGNF(dt,std_acc,r_std,rdot_std,X_initial,max_iter)
            % Constructor function to initialize the filter
            
            % Variables for Stopping criterion
            obj.max_iter = max_iter;
            
    
            % Initial State
            obj.X = X_initial;
         
            % State transition matrix
            obj.dt = dt;
            %wave number k=-lambda=c/f
            k= -299792458/94e6; 

            obj.F = [1, 0, k*dt,0;
                     0, 0, k, k*dt;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
                    
            
            obj.H = [1,0,0,0;
                     0,0,1,0;];                             % Measurement Function

            

            obj.Q = [(dt^4)/4,(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                     (dt^3)/2, dt^2, 0, 0;
                     0, 0, (dt^4)/4,(dt^3)/2;
                     0, 0, (dt^3)/2, dt^2];

            obj.R = [r_std,0;0,rdot_std];              % Measurement Uncertainty
            
            obj.P = [100,0,0,0;                              % Initial Error Covariance Matrix
                     0, 10, 0, 0;
                     0, 0, 1,0;
                     0, 0, 0, 0.1];                           

            obj.A = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];
            
            obj.wk = std_acc*[dt^2;dt;dt^2;dt];

            obj.epsDoppler =[];
            obj.epsRange = [];

        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
           
            obj.X = obj.F*obj.X;
        
            % Initial covariance matrix
            obj.P = (obj.A * obj.P * obj.A.') + obj.Q;
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            X_pred = obj.X;
            GN_Obj = obj;
        end


        function [X_est, RGNF_obj] = update(obj, Y_n)
            x_new = obj.X;

            %Convergence tolerance 
            tolerance = 1e-1; 

            %The forgetting  factor(Lambda) - between 0 and 1
            lambda = 1;

            %S = H*P*H'+ R - Innovation Covariance Matrix
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% ADAPTIVE FILTERING - Covariance scaling on Outliers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ek = Y_n - obj.H*obj.X;
            eps_doppler = (ek(2)^2)/ obj.S(2,2);
            eps_range = (ek(1)^2)/ obj.S(1,1);

            obj.epsDoppler = [obj.epsDoppler,eps_doppler];
            obj.epsRange = [obj.epsRange,eps_range];
            

            M=6; %Number of samples to average
            r_adapt=obj.R;
            if(size(obj.epsDoppler,2)>M && eps_doppler > 1 && eps_doppler>mean(obj.epsDoppler(end-M:end-1)))
                r_adapt(2,2) =(eps_doppler-1)*obj.P(3,3) + eps_doppler*r_adapt(2,2)*10;
                obj.epsDoppler =obj.epsDoppler(end-M:end-1);
            end

            if(size(obj.epsRange,2)>M && eps_range > 1 && eps_range>20*mean(obj.epsRange(end-M:end-1)))
                r_adapt(1,1) =(eps_range-1)*obj.P(1,1) + eps_range*r_adapt(1,1);
                obj.epsRange =obj.epsRange(end-M:end-1);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Recursive Gauss Newton Update Step
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           

            S1 = obj.H * obj.P * obj.H.' + r_adapt;
            I = eye(size(obj.H, 2));

            % Recursive Update
            for i = 1:obj.max_iter                
                % Observer gain Kn
                % K_n = Wn^-1_n-1*H^T[R_n + H W^-1_n-1 H^T]^-1
                K_n = obj.P * obj.H.' *S1^(-1); 
                dx = K_n*(Y_n - obj.H * x_new - obj.H * (obj.X - x_new));
                %weights = obj.huberFunction(obj.R,residual, 1.5); % Compute weights using Huber loss function
            
                % Compute the update
                x_temp = obj.X + dx;

                % Check for convergence
                if norm(x_temp-x_new) < tolerance
                    x_new = x_temp;
                    % Converged, exit the loop
                    break; 
                else
                    x_new = x_temp;
                end
            end

            % W^−1_n = λ^−1[I −Kn H]W^−1_n-1
            obj.P = (I - K_n * obj.H) * obj.P / lambda;

                            
            % Update Covariance Matrix
            obj.R = [500,0;0,0.5];
        
            obj.X = x_new;
            X_est = obj.X;
            RGNF_obj = obj;
        end
    end
    methods(Static)
        function weights = huberFunction(residual, delta)                    
            weights = zeros(size(residual));
            quadratic_indices = abs(norm(res)) <= delta;
            
            weights(quadratic_indices) = 1; % Quadratic loss
            weights(~quadratic_indices) = delta ./ abs(residual(~quadratic_indices)); % Linear loss
        end
    end
    
 end
