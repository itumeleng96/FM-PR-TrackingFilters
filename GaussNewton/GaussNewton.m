classdef GaussNewton

    properties
        max_iterations,
        tol_difference,
        tolerance,
        initial_guess,
        coefficients,
        x,y;
    end
    
    methods
        function obj = GaussNewton(max_iteration,tol_difference,tolerance,initial_guess)
            obj.max_iterations = max_iteration;
            obj.tol_difference = tol_difference;
            obj.tolerance = tolerance;
            obj.initial_guess = initial_guess;
        end
        
        function [coefficients]= fit(obj,x,y,initial_guess)
            
            % param x: independant variable (Bistatic range)
            % param y: Response vector      (Bistatic Doppler shift)
          
            obj.x = x;
            obj.y = y;

            obj.coefficients = initial_guess;
            
            rmse_prev = inf;

            for i=0:obj.max_iterations
                residual = obj.getResidual();
                jacobian = obj.calculateJacobian(obj.coefficients,10^(-6));

                obj.coefficients = obj.coefficients - obj.calculateInverse(jacobian);
                rmse = sqrt(sum(residual.^2));
                if(~isnan(obj.tol_difference)) 
                    diff = abs(rmse_prev-rmse);
                    if (diff < obj.tol_difference)
                        coefficients= obj.coefficients;
                        return
                    end

                end
                if rmse < obj.tolerance
                    coefficients =  obj.coefficients;
                    return
                end

                rmse_prev = rmse;
            end

            coefficients = obj.coefficients;

        end

        function [y]= predict(obj,x)
            
            %Predict Bistatic Doppler shift(y) given Bistatic range(x)
            y=obj.fit_function(x,obj.coefficients);
        end

        function [residual] = getResidual(obj)
            residual = obj.calculateResidual(obj.coefficients);
        end

        function [y_estimate] = getEstimate(obj)
            y_estimate = obj.fit_function(obj.x,obj.coefficients);
        end

        function [residual] = calculateResidual(obj,coefficients)
            y_fit = obj.fit_function(obj.x,coefficients);
            residual = y_fit- obj.y;
        end

        function [jacobian] =calculateJacobian(obj,x0,step)
            %Calculate Jacobian Numerically 
            %J_ij = d(r_i)/d(x_j)
            
            y0 = obj.calculateResidual(x0);
            jacobian = [];
    
            for i=size(x0,1)
                x_  = x0;
                x_(i) = x_(i) +step;

                y_ = obj.calculateResidual(x_);
                derivative = (y_-y0)/step;
                jacobian=[jacobian,derivative];
            end

            jacobian = jacobian.';
        end

        function[inverse] = calculateInverse(obj,x)
            inverse_ = pinv((x.')*x) *(x.');
            inverse = inverse_.';
            disp(inverse);
        end

        function [y]  = fit_function(obj,x,coeff)
            y = coeff(1)*x.^2 + coeff(2)*x.^2 + coeff(3);

        end

     end
 end
        
    


