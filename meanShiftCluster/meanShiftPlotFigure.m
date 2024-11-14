function [clustCent, prevCent, xVariances, yVariances, numPoints] = meanShiftPlotFigure(dataPoints, bandWidthX, bandWidthY, prevCentroids, f)
    %MEANSHIFT Summary of this function goes here
    % ---INPUT---
    % dataPts           - input data, (numDim x numPts)
    % bandWidthX        - bandwidth parameter for x-axis (scalar)
    % bandWidthY        - bandwidth parameter for y-axis (scalar)
    % prevCentroids     - previous centroids for plotting historical positions
    % f                 - figure handle for plotting
    
    % ---OUTPUT---
    % clustCent         - is locations of cluster centers (numDim x numClust)
    
    [~, numPts] = size(dataPoints);
    numClust = 0;
    bandSqX = bandWidthX^2;
    bandSqY = bandWidthY^2;
    initPtInds = 1:numPts;
    stopThresh = 1e-3 * [bandWidthX, bandWidthY]; % when mean has converged
    clustCent = [];
    beenVisitedFlag = zeros(1, numPts, 'uint8'); % track if a points been seen already
    numInitPts = numPts;
    clusterVotes = zeros(1, numPts, 'uint16'); % used to resolve conflicts on cluster membership

    while numInitPts
        tempInd = ceil((numInitPts - 1e-6) * rand); % pick a random seed point
        stInd = initPtInds(tempInd);
        myMean = dataPoints(:, stInd); % initialize mean to this points location
        myMembers = [];
        thisClusterVotes = zeros(1, numPts, 'uint16');

        while 1
            sqDistToAll = ((myMean(1) - dataPoints(1, :)).^2) / bandSqX + ((myMean(2) - dataPoints(2, :)).^2) / bandSqY;
            inInds = find(sqDistToAll < 1); % points within bandWidth
            thisClusterVotes(inInds) = thisClusterVotes(inInds) + 1;

            myOldMean = myMean;
            myMean = mean(dataPoints(:, inInds), 2);
            myMembers = [myMembers, inInds];
            beenVisitedFlag(myMembers) = 1;

            if norm(myMean - myOldMean) < stopThresh
                mergeWith = 0;
                for cN = 1:numClust
                    distToOther = norm(myMean - clustCent(:, cN));
                    if distToOther < bandWidthX / 2 && distToOther < bandWidthY / 2
                        mergeWith = cN;
                        break;
                    end
                end

                if mergeWith > 0
                    clustCent(:, mergeWith) = 0.5 * (myMean + clustCent(:, mergeWith));
                    clusterVotes(mergeWith, :) = clusterVotes(mergeWith, :) + thisClusterVotes;
                else
                    numClust = numClust + 1;
                    clustCent(:, numClust) = myMean;
                    clusterVotes(numClust, :) = thisClusterVotes;
                end
                break;
            end
        end

        initPtInds = find(beenVisitedFlag == 0);
        numInitPts = length(initPtInds);
    end

    xVariances = zeros(1, numClust);
    yVariances = zeros(1, numClust);
    numPoints = zeros(1, numClust);

    for cN = 1:numClust
        clusterMembers = clusterVotes(cN, :) > 0;
        clusterPoints = dataPoints(:, clusterMembers);
        xVariances(cN) = var(clusterPoints(1, :));
        yVariances(cN) = var(clusterPoints(2, :));
        numPoints(cN) = size(clusterPoints, 2);
    end

    % Plot centroids on CFAR Plot using the same scales and axis
    prevCent = prevCentroids;
    figure(f);
    hold on;
    if ~isempty(prevCent)
        plot(prevCent(1, :), prevCent(2, :), '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.3, 0.3, 0.3], 'MarkerSize', 4);
    end
    plot(clustCent(1, :), clustCent(2, :), 's', 'MarkerEdgeColor', 'white', 'MarkerFaceColor', 'red', 'MarkerSize', 8);

    
    xlabel('Bistatic Range [km]', 'FontSize', 18, 'Color', 'black');
    ylabel('Bistatic Doppler frequency [Hz]', 'FontSize', 18, 'Color', 'black');
    grid on;
    set(gca, 'YDir', 'normal'); % Ensure y-axis is not inverted

    prevCentroids = [prevCentroids, clustCent];
    prevCent = prevCentroids;

    hold off;
end
