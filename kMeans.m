function [ cluster, centr ] = kMeans( k, P )


numP = size(P,2); % number of points
dimP = size(P,1); % dimension of points


%% Choose k data points as initial centroids
% choose k unique random indices between 1 and size(P,2) (number of points)
randIdx = randperm(numP,k);
% initial centroids
centr = P(:,randIdx);


%% Repeat until stopping criterion is met

% init cluster array
cluster = zeros(1,numP);

% init previous cluster array clusterPrev (for stopping criterion)
clusterPrev = cluster;

% for reference: count the iterations
iterations = 0;

% init stopping criterion
stop = false; % if stopping criterion met, it changes to true

while stop == false
    
    % for each data point 
    for idxP = 1:numP
        % init distance array dist
        dist = zeros(1,k);
        % compute distance to each centroid
        for idxC=1:k
            dist(idxC) = norm(P(:,idxP)-centr(:,idxC));
        end
        % find index of closest centroid (= find the cluster)
        [~, clusterP] = min(dist);
        cluster(idxP) = clusterP;
    end
    
    % Recompute centroids using current cluster memberships:
    
    centr = zeros(dimP,k);
    % for every cluster compute new centroid
    for idxC = 1:k
        % find the points in cluster number idxC and compute row-wise mean
        centr(:,idxC) = mean(P(:,cluster==idxC),2);
    end
    
    if clusterPrev==cluster
        stop = true;
    end
    
    clusterPrev = cluster;
    
    iterations = iterations + 1;
    
end

end