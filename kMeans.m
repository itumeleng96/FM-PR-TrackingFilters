function [ cluster, centroids ] = kMeans( k, dataset )

numberOfPoints = size(dataset,2); 
dimensionOfPoints = size(dataset,1); 

randomIndices = randperm(numberOfPoints,k);
centroids = dataset(:,randomIndices);

cluster = zeros(1,numberOfPoints);

clusterPrev = cluster;

iterations = 0;

stop = false;

while stop == false    
    % for each data point 
    for indexPoint = 1:numberOfPoints
        
        dist = zeros(1,k);
        % compute distance to each centroid
        for indexCluster=1:k
            dist(indexCluster) = norm(dataset(:,indexPoint)-centroids(:,indexCluster));
        end
        % find index of closest centroid (= find the cluster)
        [~, clusterP] = min(dist);
        cluster(indexPoint) = clusterP;
    end
    
    % Recompute centroids using current cluster memberships:
    
    centroids = zeros(dimensionOfPoints,k);
    for indexC = 1:k
        centroids(:,indexC) = mean(dataset(:,cluster==indexC),2);
    end
    
    if clusterPrev==cluster
        stop = true;
    end
    
    clusterPrev = cluster;
    iterations = iterations + 1;
    
end
end