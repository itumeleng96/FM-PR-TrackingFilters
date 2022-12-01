function [clustCent] = meanShiftPlot(dataPoints,bandWidth,fs,fd_max,td_max)
%MEANSHIFT Summary of this function goes here

% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% bandWidth         - is bandwidth parameter (scalar)

% ---OUTPUT---
% clustCent         - is locations of cluster centers (numDim x numClust)

%**** Initialize stuff ***
[numDim,numPts] = size(dataPoints);
numClust        = 0;
bandSq          = bandWidth^2;
initPtInds      = 1:numPts;
stopThresh      = 1e-3*bandWidth;                                %when mean has converged
clustCent       = [];                                            %center of clust
beenVisitedFlag = zeros(1,numPts,'uint8');                       %track if a points been seen already
numInitPts      = numPts;                                        %number of points to posibaly use as initilization points
clusterVotes    = zeros(1,numPts,'uint16');                      %used to resolve conflicts on cluster membership

Ndelay = floor(td_max*fs);   %number of points corresponding to td_max
time = 0:1/fs:Ndelay/fs;
frequency = -fd_max:1:fd_max;

while numInitPts
    tempInd         = ceil( (numInitPts-1e-6)*rand);            %pick a random seed point
    stInd           = initPtInds(tempInd);                      %use this point as start of mean
    myMean          = dataPoints(:,stInd);                      % intilize mean to this points location
    myMembers       = [];                                       % points that will get added to this cluster                          
    thisClusterVotes    = zeros(1,numPts,'uint16');             %used to resolve conflicts on cluster membership
    while 1     %loop untill convergence
        
        sqDistToAll = sum((repmat(myMean,1,numPts) - dataPoints).^2);    %dist squared from mean to all points still active
        inInds      = find(sqDistToAll < bandSq);               %points within bandWidth
        thisClusterVotes(inInds) = thisClusterVotes(inInds)+1;  %add a vote for all the in points belonging to this cluster
        
        
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
                if distToOther < bandWidth/2                    %if its within bandwidth/2 merge new and old
                    mergeWith = cN;
                    break;
                end
            end
            
            
            if mergeWith > 0    % something to merge
                clustCent(:,mergeWith)       = 0.5*(myMean+clustCent(:,mergeWith));             %record the max as the mean of the two merged (I know biased twoards new ones)
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

%Plot centroids on CFAR Plot
hold on;
plot(time((round(clustCent(1,:)))),frequency(round(clustCent(2,:))),'^','MarkerFaceColor','black', 'MarkerSize', 5);

end

