classdef HuberScore
    properties
        delta % Delta parameter for Huber score
    end
    
    methods
        function obj = HuberScore(delta)
            obj.delta = delta;
        end
        
        function score = evaluate(obj, residual)
            % Evaluate the Huber score for a given residual
            score = zeros(size(residual));
            for i = 1:numel(residual)
                if abs(residual(i)) <= obj.delta
                    score(i) = 0.5 * residual(i)^2;
                else
                    score(i) = obj.delta * (abs(residual(i)) - 0.5 * obj.delta);
                end
            end
        end
    end
end
