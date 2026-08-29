classdef CSRGNF

    properties
        dt, X, F, H, Q, Q_nominal, R, R_nominal, P, S, max_iter, wk, lambda, ...
        cs_alpha, cs_adapt_q;
    end
    
    methods
        function obj = CSRGNF(dt,std_acc,r_std,rdot_std,X_initial,max_iter,lambda)
            % Constructor function to initialize the filter
            
            % Variables for Stopping criterion
            obj.max_iter = max_iter;
            
    
            % Initial State
            obj.X = X_initial;
         
            obj.dt = dt;

            obj.F = [1, dt, 0, 0;
                     0, 1, 0, 0;
                     0, 0, 1, dt;
                     0, 0, 0, 1;];

            obj.H = [1,0,0,0;
                     0,0,1,0;];                            % Measurement Function

            

            obj.Q = [std_acc(1)*(dt^4)/4,std_acc(1)*(dt^3)/2,0,0;                % Process Noise Covariance Matrix
                 std_acc(1)*(dt^3)/2, std_acc(1)*dt^2, 0, 0;
                 0, 0, std_acc(2)*(dt^4)/4,std_acc(2)*(dt^3)/2;
                 0, 0, std_acc(2)*(dt^3)/2, std_acc(2)*dt^2];


            %Measurement Error covariance matrix
            obj.R = [r_std,0;
                     0,rdot_std;];
            obj.R_nominal = obj.R;                 % Nominal R for fading-memory decay

            obj.P = [5,0,0,0;                                      % Initial Error Covariance Matrix
                     0, 2.5, 0, 0;
                     0, 0, 2,0;
                     0, 0, 0, 1];  
            
            obj.wk = [std_acc(1)*dt^2;std_acc(1)*dt;std_acc(2)*dt^2;std_acc(2)*dt];
            obj.lambda = lambda;

            % Akhlaghi (2018) covariance-scaling defaults (settable post-construction).
            obj.cs_alpha   = 0.98;   % forgetting factor (both R and Q)
            obj.cs_adapt_q = true;   % master switch for Q adaptation
            obj.Q_nominal  = obj.Q;  % nominal Q retained as PSD floor

        end

        % Function to predict the next state
        function [X_pred, GN_Obj] = predict(obj)
           
            obj.X = obj.F*obj.X;
        
            % Initial covariance matrix
            obj.P = (obj.F * obj.P * obj.F.') + obj.Q;
            obj.S = obj.H * obj.P * obj.H.' + obj.R;

            X_pred = obj.X;
            GN_Obj = obj;
        end


        function [X_est, RGNF_obj] = update(obj, Y_n)
            % Akhlaghi (2018) verbatim continuous R + Q inside RGNF iteration.
            % Outlier rejection handled upstream at the association step.
            %   R <- alpha R + (1-alpha)( eps eps' + H P^- H' )       Eq. 11
            %   Q <- alpha Q + (1-alpha)( K d d' K' )                 Eq. 15

            alphaFactor = obj.cs_alpha;
            x_new       = obj.X;
            tolerance   = 1e-1;

            % Innovation (pre-update) using current R.
            obj.S = obj.H * obj.P * obj.H.' + obj.R;
            ek    = Y_n - obj.H * obj.X;
            P_pred = obj.P;

            % Recursive Gauss-Newton update.
            S1 = obj.H * obj.P * obj.H.' + obj.R;
            I  = eye(size(obj.H, 2));
            for i = 1:obj.max_iter
                K_n = obj.P * obj.H.' * S1^(-1);
                dx = K_n * (Y_n - obj.H * x_new - obj.H * (obj.X - x_new));
                x_temp = obj.X + dx;
                if norm(x_temp - x_new) < tolerance
                    x_new = x_temp;
                    break;
                else
                    x_new = x_temp;
                end
            end
            obj.P = (I - K_n * obj.H) * obj.P / obj.lambda;
            obj.X = x_new;

            % Akhlaghi Eq. 11 — continuous R update.
            eps_post = Y_n - obj.H * obj.X;
            HPH      = obj.H * P_pred * obj.H.';
            obj.R    = alphaFactor * obj.R + (1 - alphaFactor) * ...
                       (eps_post * eps_post.' + HPH);

            % Akhlaghi Eq. 15 — continuous Q update (per block).
            if obj.cs_adapt_q
                KddK = K_n * (ek * ek.') * K_n.';
                obj.Q(1:2,1:2) = alphaFactor * obj.Q(1:2,1:2) + ...
                                 (1 - alphaFactor) * KddK(1:2,1:2);
                obj.Q(3:4,3:4) = alphaFactor * obj.Q(3:4,3:4) + ...
                                 (1 - alphaFactor) * KddK(3:4,3:4);
                obj.Q(1:2,1:2) = max(obj.Q(1:2,1:2), obj.Q_nominal(1:2,1:2));
                obj.Q(3:4,3:4) = max(obj.Q(3:4,3:4), obj.Q_nominal(3:4,3:4));
            end

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
