
function [clustCent,xVariances,yVariances,numPoints] = meanShift(dataPoints,bandWidthX,bandWidthY)
%MEANSHIFT Summary of this function goes here
% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% bandWidthX        - bandwidth parameter for x-axis (scalar)
% bandWidthY        - bandwidth parameter for y-axis (scalar)

% ---OUTPUT---
% clustCent         - is locations of cluster centers (numDim x numClust)

[~,numPts] = size(dataPoints);
numClust        = 0;
bandSqX          = bandWidthX^2;
bandSqY          = bandWidthY^2;
initPtInds      = 1:numPts;
stopThresh      = 1e-3*[bandWidthX, bandWidthY];              %when mean has converged
clustCent       = [];                                         %center of clust
beenVisitedFlag = zeros(1,numPts,'uint8');                    %track if a points been seen already
numInitPts      = numPts;                                     %number of points to posibaly use as initilization points
clusterVotes    = zeros(1,numPts,'uint16');                   %used to resolve conflicts on cluster membership

while numInitPts
    tempInd         = ceil( (numInitPts-1e-6)*rand);           %pick a random seed point
    stInd           = initPtInds(tempInd);                     %use this point as start of mean
    myMean          = dataPoints(:,stInd);                     % intilize mean to this points location
    myMembers       = [];                                      % points that will get added to this cluster                          
    thisClusterVotes    = zeros(1,numPts,'uint16');            %used to resolve conflicts on cluster membership
    while 1     %loop untill convergence
        
        sqDistToAll = ((myMean(1)-dataPoints(1,:)).^2)/bandSqX + ((myMean(2)-dataPoints(2,:)).^2)/bandSqY;    %dist squared from mean to all points still active
        inInds      = find(sqDistToAll < 1);                    %points within bandWidth
        thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;   %add a vote for all the in points belonging to this cluster
        
        myOldMean   = myMean;                                   %save the old mean
        myMean      = mean(dataPoints(:,inInds),2);             %compute the new mean
        myMembers   = [myMembers inInds];                       %add any point within bandWidth to the cluster
        beenVisitedFlag(myMembers) = 1;                         %mark that these points have been visited
        
        %**** if mean doesn't move much stop this cluster ***
        if norm(myMean-myOldMean) < stopThresh
            %check for merge posibilities
            mergeWith = 0;
            for cN = 1:numClust
                distToOther = norm(myMean-clustCent(:,cN));     %distance from posible new clust max to old clust max
                if distToOther < bandWidthX/2 && distToOther < bandWidthY/2   %if its within bandwidth/2 merge new and old
                    mergeWith = cN;
                    break;
                end
            end
            
            
            if mergeWith > 0    % something to merge
                clustCent(:,mergeWith)       = 0.5*(myMean+clustCent(:,mergeWith));             %record the max as the mean of the two merged
                %clustMembsCell{mergeWith}    = unique([clustMembsCell{mergeWith} myMembers]);   %record which points inside 
                clusterVotes(mergeWith,:)    = clusterVotes(mergeWith,:) + thisClusterVotes;    %add these votes to the merged cluster
            else    %its a new cluster
                numClust                    = numClust+1;                                           %increment clusters
                clustCent(:,numClust)       = myMean;                                               %record the mean  
                %clustMembsCell{numClust}    = myMembers;                                           %store my members
                clusterVotes(numClust,:)    = thisClusterVotes;
            end
            break;
        end
    end
    
    
    initPtInds      = find(beenVisitedFlag == 0);           %we can initialize with any of the points not yet visited
    numInitPts      = length(initPtInds);                   %number of active points in set
end

xVariances = zeros(1, numClust);
yVariances = zeros(1, numClust);
numPoints = zeros(1,numClust);

% Calculate and print x and y variance for each cluster
for cN = 1:numClust
    clusterMembers = clusterVotes(cN, :) > 0;
    clusterPoints = dataPoints(:, clusterMembers);

    xVariances(cN) = var(clusterPoints(1, :));
    yVariances(cN) = var(clusterPoints(2, :));
    numPoints(cN) = size(clusterPoints,2);
end


end

