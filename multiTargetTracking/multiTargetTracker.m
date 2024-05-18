classdef multiTargetTracker
    %MTT Multi-Tracker
    
    properties
        tracks,
        confirmationThreshold,
        deletionThreshold,
        gatingThreshold,       %Radius around the predicted measurement to eliminate other measurements
        filterType,            %KalmanFilter:1 , GaussNewton:2
        newtracksCreated,
        idCounter,
        %M,                      %M-out-of N logic
        %N,                      %M-out-of N logic

    end
    
    methods
        function obj = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType)
            %MTT Construct an instance of this class
            obj.confirmationThreshold = confirmationThreshold;
            obj.deletionThreshold = deletionThreshold;
            obj.gatingThreshold = gatingThreshold;
            obj.filterType = filterType;
            obj.newtracksCreated = 0;
            obj.idCounter =0;
            %obj.M =3;
            %obj.N =5;
        end
        
        function obj = createNewTracks(obj,detections,index)
            %This method assigns Detections to the nearest Track, else
            numberOfDetections=size(detections,2);
            if isempty(obj.tracks)
                %disp("Assigning Detections");
                %create tracks = number of detections for the first time

                for i=1:numberOfDetections
                    if i ==1
                        obj.idCounter =obj.idCounter+1;
                        obj.tracks = [track([detections(1,i);detections(2,i)],[;],obj.idCounter,index,0,0,0,obj.filterType)];
                        obj.newtracksCreated = 1;
                    
                    else
                        obj.idCounter =obj.idCounter+1;
                        obj.tracks(end+1)=track([detections(1,i);detections(2,i)],[;],obj.idCounter,index,0,0,0,obj.filterType);
                        obj.newtracksCreated = 1;
                    end
                end            
            end
        end

        function obj = updateStage(obj,detections,index)
            %disp("UpdateStage");
            %Assign Tracks to Detection using GNN and update filter with new measurements
            %Get qualifying detections within radius if not create new tracks

            if ~isempty(obj.tracks)
                numOfTracks = length(obj.tracks);
                for i=1:numOfTracks
                    if(~obj.tracks(i).deleted)
                        predictedCoodinate = obj.tracks(i).predictedTrack(:,end);
                        detectionsInRadius = obj.pruneDetections(detections,predictedCoodinate,obj.gatingThreshold);
                        [detection,detections] = obj.globalNearestNeighbour(detectionsInRadius,predictedCoodinate,detections);
                        if(detection)
                            obj.tracks(i)=obj.tracks(i).updateTrueTrack(detection); 
                        end
                    end
                end   
                %If detection is unassigned,create new track 
                
                if ~isempty(detections)
                    numberOfUnassignedDetections = size(detections,2);
                    
                    for i=1:numberOfUnassignedDetections
                        obj.idCounter = obj.idCounter+1;
                        obj.tracks(end+1) = track([detections(1,i);detections(2,i)],[;],obj.idCounter,index,0,0,0,obj.filterType); 
                    end
                end
            end
            obj.newtracksCreated = 0;
        end

        function tracks = deleteTracks(obj)
            %Delete Tracks based on deletion Treshold
            %idx_to_delete =[];
            for i=1:max(size(obj.tracks))
                %if obj.tracks(i).seenCountDel > obj.M
                if obj.tracks(i).sampleSinceLastUpdate > obj.deletionThreshold
                   %idx_to_delete = [idx_to_delete, i];
                   obj.tracks(i).deleted =1;
                end
            end
            tracks = obj.tracks;
        end
        
        function tracks = confirmTracks(obj)
            %Confirm Tracks based on confirmation Threshold
            for i=1:max(size(obj.tracks))                
                %if (obj.tracks(i).confirmed == 0 && obj.tracks(i).seenCount >= obj.M)
                if (obj.tracks(i).confirmed == 0 && obj.tracks(i).numberOfUpdates > obj.confirmationThreshold)
                    obj.tracks(i).confirmed = 1;
                end
            end

            tracks =obj.tracks;
        end

        function obj = maintainTracks(obj)
            %disp("Maintain Tracks");
            obj.tracks =obj.deleteTracks();
            obj.tracks =obj.confirmTracks();

        end

        function obj = predictionStage(obj)
            %disp("Prediction Stage");
            numberOfTracks = max(size(obj.tracks));
            
            for i=1:numberOfTracks
                if(~obj.tracks(i).deleted)
                    obj.tracks(i)=obj.tracks(i).predictTrack();
                end
            end

            
        end

        function obj = plotMultiTargetTracking(obj,fs,fd_max,td_max,index,f,RDM)
                figure(f);
                c=3e8;
                Ndelay = floor(td_max*fs);                                 
                time = 0:1/fs:Ndelay/fs;
                range = time *c;
                frequency = -fd_max:1:fd_max;
                imagesc(range,frequency,RDM*0);
                colormap(gca, 'white'); % Set the colormap to 'gray'
    
                text(0,10,"Time:" + index+ "s");
                axis xy;
                xlabel('Bistatic Range [m]','Fontsize',10);
                ylabel('Bistatic Doppler frequency [Hz]','Fontsize',10);
                grid on;
                title('Targets centroids and  Prediction');
                
               
                hold on;
                for i = 1:length(obj.tracks)
                    if(obj.tracks(i).deleted ==0)
                        trackID = obj.tracks(i).trackId;
                        text(obj.tracks(i).predictedTrack(1,end)+1000, obj.tracks(i).predictedTrack(2,end)-5, num2str(trackID), 'Color', 'red', 'FontSize', 8);
                        if obj.tracks(i).confirmed == 0
                            plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor','none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Predicted Track');                        % Plot true track as open circles joined by a line
                            plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color','green' ,'MarkerSize', 4, 'DisplayName', 'Confirmed Track');
                        else
                            plot(obj.tracks(i).predictedTrack(1, :), obj.tracks(i).predictedTrack(2, :), '^', 'MarkerFaceColor', 'none','MarkerEdgeColor', 'blue', 'MarkerSize', 6, 'DisplayName', 'Predicted Track');                        % Plot true track as open circles joined by a line
                            plot(obj.tracks(i).trueTrack(1,:), obj.tracks(i).trueTrack(2,:), '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4, 'DisplayName', 'Tentative Track');
                        end
                    end
                    % Shift and box the text annotations
                    %textOffset = 15; 
                    %annotationX = obj.tracks(i).predictedTrack(1, end) + textOffset;
                    %annotationY = obj.tracks(i).predictedTrack(2, end) + textOffset;
                    
                    %rectangle('Position', [annotationX-0.5, annotationY-0.5, 1, 10], 'EdgeColor', 'black', 'LineWidth', 1); % Box
                    %text(annotationX, annotationY, ['T',num2str(i)], 'Color', 'black', 'FontSize', 8);

                end
                predicted_marker = plot(nan, nan, '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue', 'MarkerSize', 6);
                tentative_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'green', 'MarkerSize', 4);
                confirmed_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4);
                
                % Create a legend with custom markers and labels
                legend([predicted_marker, tentative_marker, confirmed_marker], 'Predicted Track', 'Tentative Track', 'Measurement Track', 'Location', 'best');
                %hold on;
        end

        function plotErrors(obj,f,f2,index,rangeGroundTruth,dopplerGroundTruth)
            % Plot Range
            figure(f);
            hold on;
            title('Error comparison Range(m)');
            xlabel('Bistatic Range [m]','Fontsize',10);
            ylabel('Bistatic Doppler [Hz]','Fontsize',10);
            plot(rangeGroundTruth(1:index),dopplerGroundTruth(1:index),'-*');
            
            hold on;
            for i = 1:length(obj.tracks)
                % Check if the track has sufficient length
                if  obj.tracks(i).confirmed==1
                    % Plot the track
                    plot(obj.tracks(i).predictedTrack(1, :),obj.tracks(i).predictedTrack(2, :),'-^');
                    plot(obj.tracks(i).trueTrack(1, :),obj.tracks(i).trueTrack(2, :),'-o');

                end
            end
                ground__marker = plot(nan, nan, '-*');
                predicted_marker = plot(nan, nan, '-^');
                measured_marker = plot(nan, nan, '-o');
                legend([ground__marker, predicted_marker, measured_marker], 'Ground Truth', 'Predicted Track', 'Measured Track', 'Location', 'best');
                grid on;
            hold off;
        end

        function plotErrorTime(obj, f, f2, index,rangeGroundTruth,dopplerGroundTruth)
            % Plot Range
            figure(f);
            hold on;
            title('Error comparison Range(m)');
            xlabel('Time(s)','Fontsize',10);
            ylabel('Bistatic Range [m]','Fontsize',10);
            
            plot(rangeGroundTruth(1:index),'-*');
            for i = 1:length(obj.tracks)
                % Check if the track has sufficient length
                if size(obj.tracks(i).predictedTrack, 2) >= index && obj.tracks(i).confirmed==1
                    % Plot the track
                    plot(obj.tracks(i).predictedTrack(1, :),'-^');
                    plot(obj.tracks(i).trueTrack(1, :),'-');

                else
                    if size(obj.tracks(i).predictedTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).predictedTrack, 2))+1;
                        endIndex = min(index, size(obj.tracks(i).predictedTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).predictedTrack(1, :),'-^');
                    end
                    if size(obj.tracks(i).trueTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).trueTrack, 2))+1;
                        endIndex = min(index, size(obj.tracks(i).trueTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).trueTrack(1, :),'-');
                    end

                end
            end
            hold off;
        
            % Plot Doppler
            figure(f2);
            hold on;
            title('Error comparison Doppler (Hz)');
            xlabel('Time(s)','Fontsize',10);
            ylabel('Bistatic Doppler [Hz]','Fontsize',10);
            plot(dopplerGroundTruth(1:index),'-*');
            for i = 1:length(obj.tracks)
                % Check if the track has sufficient length
                if size(obj.tracks(i).predictedTrack, 2) >= index && obj.tracks(i).confirmed==1
                    plot(obj.tracks(i).predictedTrack(2, :),'-^');
                    plot(obj.tracks(i).trueTrack(2, :),'-');
                else
                    if size(obj.tracks(i).predictedTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).predictedTrack, 2) + 1);
                        endIndex = min(index, size(obj.tracks(i).predictedTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).predictedTrack(2, :),'-^');
                    end
                   if size(obj.tracks(i).trueTrack,2) >1 && obj.tracks(i).confirmed==1
                        startIndex = max(1, index - size(obj.tracks(i).trueTrack, 2) + 1);
                        endIndex = min(index, size(obj.tracks(i).trueTrack, 2))-1;
                        plot(startIndex:startIndex+endIndex,obj.tracks(i).trueTrack(2, :),'-');
                    end 
                end
            end
            hold off;
        end

        function [doppler_error, range_error,doppler_meas,range_meas] =getErrors(obj, i,doppler_error,range_error)
            
            % Calculate the Range and Doppler RMS
            predictedTrack = obj.tracks(1).predictedTrack;
            trueTrack = obj.tracks(1).trueTrack;
            
            range_error(1,1:i) = predictedTrack(1, 1:i);
            doppler_error(1,1:i)= predictedTrack(2, 1:i);

            range_meas(1,1:i) = trueTrack(1, 1:i);
            doppler_meas(1,1:i)= trueTrack(2, 1:i);
                
        end


        function [doppler_ll, range_ll] = calculateLogLikelihood(obj, i,doppler_ll,range_ll)

            %for j = 1:length(obj.tracks)
            for j = 1:1
                predictedTrack = obj.tracks(j).predictedTrack;
                trueTrack = obj.tracks(j).trueTrack;
                s_matrix = obj.tracks(j).trackingFilterObject.S;
                %----------------------------------------------------------------%
                %--Log-likelihood for Bistatic Range
                %----------------------------------------------------------------%
                range_sample=trueTrack(1,i);
                range_mean =predictedTrack(1,i);
                range_ll(j, i) = logLikelihood(range_mean,s_matrix(1,1),range_sample);

                %----------------------------------------------------------------%
                %--Log-likelihood for Bistatic Doppler
                %----------------------------------------------------------------%
                doppler_sample=trueTrack(2,i);
                doppler_mean =predictedTrack(2,i);
                
                doppler_ll(j, i) = logLikelihood(doppler_mean,s_matrix(2,2),doppler_sample);
            end
        end

        function [doppler_ll, range_ll,t] = calculateLogLikelihoodGroundTruth(obj,trackId,doppler_ll,range_ll,dopplerGroundTruth,rangeGroundTruth,simTime)
            
             indexOfTrack = -1;
             Tracks = obj.tracks;
             for j = 1:length(Tracks)
                %look for track with the specified Id
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end
             t=[];
             if(indexOfTrack>0)

                 sMatrix = Tracks(indexOfTrack).sMatrix;
                 %sMatrix(1,end+1) = Tracks(indexOfTrack).trackingFilterObject.S(1,1);
                 %sMatrix(2,end) = Tracks(indexOfTrack).trackingFilterObject.S(2,2);
                 
                 disp(sMatrix);
                 disp(dopplerGroundTruth);

                 predictedTrack = Tracks(indexOfTrack).predictedTrack;
                 disp(predictedTrack(2,:));

                 range_ll = obj.logLikelihoodMatrix(predictedTrack(1,:),rangeGroundTruth(1,:),sMatrix(1,:),Tracks(indexOfTrack).deleted);
                 doppler_ll = obj.logLikelihoodMatrix(predictedTrack(2,:),dopplerGroundTruth(1,:),sMatrix(2,:),Tracks(indexOfTrack).deleted);
                 %Plot against time 
                 startTime = Tracks(indexOfTrack).startTime;
                 endTime = length(doppler_ll)+startTime -1;
                 
                 t = startTime:1:endTime;
                 
             end
        end

        function [predictedTrack] = getTrack(obj,trackId)
            
             indexOfTrack = -1;
             Tracks = obj.tracks;
             predictedTrack =[];
             for j = 1:length(Tracks)
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end

             if(indexOfTrack>0)
                 predictedTrack = Tracks(indexOfTrack).predictedTrack;
             end
        end

        function [measuredTrack] = getMeasuredTrack(obj,trackId)
            
             indexOfTrack = -1;
             Tracks = obj.tracks;
             measuredTrack =[];
             for j = 1:length(Tracks)
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end

             if(indexOfTrack>0)
                 measuredTrack = Tracks(indexOfTrack).trueTrack;
             end
        end

        function [doppler_ll, range_ll] = plotLogLikelihood(obj, f, f1,trackId,dopplerGroundTruth,rangeGroundTruth,simTime)
                        
             indexOfTrack = -1;
             Tracks = obj.tracks;
             for j = 1:length(Tracks)
                %look for track with the specified Id
                if(obj.tracks(j).trackId == trackId)
                    indexOfTrack = j;
                end
             end

             if(indexOfTrack>0)

                 sMatrix = Tracks(indexOfTrack).sMatrix;
                 sMatrix(1,end+1) = Tracks(indexOfTrack).trackingFilterObject.S(1,1);
                 sMatrix(2,end) = Tracks(indexOfTrack).trackingFilterObject.S(2,2);

                 predictedTrack = Tracks(indexOfTrack).predictedTrack;
                 disp(Tracks(indexOfTrack).trueTrack);
                 range_ll = obj.logLikelihoodMatrix(predictedTrack(1,:),rangeGroundTruth(1,:),sMatrix(1,:),Tracks(indexOfTrack).deleted);
                 doppler_ll = obj.logLikelihoodMatrix(predictedTrack(2,:),dopplerGroundTruth(1,:),sMatrix(2,:),Tracks(indexOfTrack).deleted);

                 %Plot against time 
                 startTime = Tracks(indexOfTrack).startTime;
                 endTime = length(doppler_ll)+startTime-1;

                 figure(f);
                 hold on;
                 t = startTime:1:endTime;

                 plot(t,doppler_ll);
                
                 xlabel('Time(s)');
                 ylabel('Bistatic Doppler Log-likelihood');
                 title(['Bistatic Doppler Log-likelihood vs Time for Track ID: ', num2str(trackId)]);
                 %legend('Doppler Log-likelihood', 'Doppler Error'); % Add legend for the two plotted lines
    
                 figure(f1);
                 plot(t,range_ll);
                 xlabel('Time(s)');
                 ylabel('Bistatic Range Log-likelihood');
                 title(['Bistatic Range  Log-likelihood vs Time for Track ID:',num2str(trackId)]);
                 %legend('Range Log-likelihood', 'Range Error'); % Add legend for the two plotted lines
                 %}
             end
        end
        
        function [doppler_ll, range_ll] = plotLogLikelihoodSingle(obj, f, f1, i,doppler_ll,range_ll,dopplerTrueData,rangeTrueData,plotResults)
            
            if(plotResults)
                time = 1:1:i;            
            
                 for j = 1:1
                    predictedTrack = obj.tracks(j).predictedTrack;
                    %trueTrack = obj.tracks(j).trueTrack;
                    s_matrix = obj.tracks(j).trackingFilterObject.S;
                    %disp(s_matrix);
                    %----------------------------------------------------------------%
                    %--Log-likelihood for Bistatic Range
                    %----------------------------------------------------------------%
                    range_sample=rangeTrueData(i);
                    range_mean =predictedTrack(1,i);
                    
                    %range_ll(j, i) = normlike([(range_mean),s_matrix(1,1)],range_sample);
                    %range_ll(j, i) = lognlike([range_mean,sqrt(R(1,1))],range_sample);
                    range_ll(j, i) = logLikelihood(range_mean,s_matrix(1,1),range_sample);
    
                    %----------------------------------------------------------------%
                    %--Log-likelihood for Bistatic Doppler
                    %----------------------------------------------------------------%
                    doppler_sample=dopplerTrueData(i);
                    doppler_mean =predictedTrack(2,i);
                    
                    %doppler_ll(j, i) = normlike([doppler_mean,s_matrix(2,2)],doppler_sample);
                    %doppler_ll(j, i) = lognlike([doppler_mean,sqrt(R(2, 2))],doppler_sample);
                    doppler_ll(j, i) = logLikelihood(doppler_mean,s_matrix(2,2),doppler_sample);
    
                end
    
                figure(f);
                plot(time, doppler_ll);
                xlabel('Time(s)');
                ylabel('Bistatic Doppler Log-likelihood');
                title('Bistatic Doppler Log-likelihood vs Time');
                %legend('Doppler Log-likelihood', 'Doppler Error'); % Add legend for the two plotted lines
    
                figure(f1);
                plot(time,range_ll);
                xlabel('Time(s)');
                ylabel('Bistatic Range Log-likelihood');
                title('Bistatic Range  Log-likelihood vs Time');
                %legend('Range Log-likelihood', 'Range Error'); % Add legend for the two plotted lines

            end
        end

        function [range_mse,doppler_mse] = calculateMSE(obj,i,range_mse,doppler_mse,true_doppler,true_range)

        
            for j = 1:length(obj.tracks)
                predictedTrack = obj.tracks(j).predictedTrack;
    
                % Calculate squared errors for each track at the current time step
                range_mse(i) = (true_range(i) - predictedTrack(1, i))^2;
                doppler_mse(i) = (true_doppler(i) - predictedTrack(2, i))^2;
    
            end
                

        
        end
    
        function [crlb_doppler,crlb_range] = calculateCRLB(obj, i,crlb_doppler,crlb_range,true_doppler,true_range)
            
            % Calculate the CRLB for each track
            %{
            for j = 1:length(obj.tracks)
                observedTrack = obj.tracks(j).trueTrack;
                sample_variance_doppler = sum((observedTrack(2,i) - true_doppler(i)  ).^2) / 1;
                sample_variance_range = sum((observedTrack(1,i) - true_range(i)  ).^2) / 1;

                variance_mu_doppler = sample_variance_doppler / 1; % Variance of the sample mean
                variance_mu_range = sample_variance_range / 1; % Variance of the sample mean


                Fisher_information_mu_doppler = 1 / variance_mu_doppler;
                Fisher_information_mu_range = 1 / variance_mu_range;

                crlb_doppler(i) = 1/Fisher_information_mu_doppler;
                crlb_range(i) = 1/Fisher_information_mu_range;
                
            end
            %}
            FIM = obj.tracks(1).trackingFilterObject.P;
            crlb = FIM^(-1) ;
            crlb_doppler(i) = crlb(1,1);
            crlb_range(i) = crlb(2,2);

        end
    end

    methods(Static)

        function detectionsInRadius = pruneDetections(detections, predictedCoordinate, gatingThreshold)
            % For every detection check that it falls within the predicted coordinates' radius
            % Rectangular (1-norm)  Gating
         
            numberOfDetections = size(detections, 2);
            detectionsInRadius = [];
            x_std = gatingThreshold(1);
            y_std = gatingThreshold(2);
            k = 3;  %Threshold

            for i = 1:numberOfDetections
                y_dist = abs(detections(2, i) - predictedCoordinate(2));
                x_dist = abs(detections(1, i) - predictedCoordinate(1));
                
                if x_dist<k*x_std && y_dist<k*y_std
                    detectionsInRadius = [detectionsInRadius, detections(:, i)];
                end
            end
        end

        function [detection,detections] = globalNearestNeighbour(detectionsInRadius,predictedCoordinate,allDetections)
            numberOfDetections=size(detectionsInRadius,2);
            distances=[];
            measurementCovariance = [100,0;0,0.1];
            %Use Mahalanobis Distance
            for i=1:numberOfDetections
                %distances(i)=norm(detectionsInRadius(:,i)-predictedCoordinate);
                innovation = detectionsInRadius(:,i)-predictedCoordinate;
                %disp("Innovation");
                %disp(innovation);
                mahalanobisDistance = sqrt(innovation' / measurementCovariance*innovation);
                %disp(mahalanobisDistance);
                distances(i)=mahalanobisDistance;
            end

            %Get The index of the detection with min distances and delete from detections
            [~,indexOfNeighbour] = min(distances);
            detection = detectionsInRadius(:,indexOfNeighbour); 
            %Find index of detection and delete
            if detection
                indexInAllDetections = round(allDetections(1,:))==round(detection(1,1)) & round(allDetections(2,:))==round(detection(2,1));
                allDetections(:,indexInAllDetections) = [];
            end

            detections = allDetections;
         end

         function [logLikelihoodValues] = logLikelihoodMatrix(samples,groundTruthValues,s_matrix,deleted)
            %Loop through the length of the track 

            logLikelihoodValues =[];
            if(size(samples,2)== size(groundTruthValues,2))
                for j=1:size(s_matrix,2)
                    logLikelihoodValues = [logLikelihoodValues,logLikelihood(groundTruthValues(j),s_matrix(j),samples(j))];
                end
            end

            %If the predicted values of the track are not the same length as groundTruthValues
            if(size(samples,2) < size(groundTruthValues,2) && deleted)
                %if deleted it means the missing values are at the end
                for j=1:size(s_matrix,2)
                    logLikelihoodValues = [logLikelihoodValues,logLikelihood(groundTruthValues(j),s_matrix(j),samples(j))];
                end
            end

            if(size(samples,2)< size(groundTruthValues,2) && ~deleted)
                %if deleted it means the missing values are at the
                %beginning
                for j=1:size(s_matrix,2)
                    %%Fix time ,so the plot does not start from time=0
                    diff = abs(size(samples,2)-size(groundTruthValues,2));
                    logLikelihoodValues = [logLikelihoodValues,logLikelihood(groundTruthValues(diff+j),s_matrix(j),samples(j))];
                end
            end


         end

    end

end