function [track_] = trackAssign(track,centroids)

    N = size(centroids,2);
    %Check if Track is empty
    if isempty(track)
        track_=[];
        for i =1:N
            centroids_ =[];
            centroids_(1,1)= centroids(1,i);
            centroids_(2,1)= centroids(2,i);
            
            %Add The centroids to the 3D track matrix 
            track_(:,:,i)=centroids_;
        end
    else
        %for each centroid get the minimum distance
        track_=[];
        for i=1:N
            %Get the centroids from 
            euclidTrack= track;
            
            Distances = zeros();
            
            for j=1:N
                %Get euclidean distance
                euclidCentroids_j = euclidTrack(:,:,j);
                Distances(j)=((centroids(1,i) - euclidCentroids_j(1,end))^2+(centroids(2,i) - euclidCentroids_j(2,end))^2)^0.5;
            end
            %Get the min value and assign to track
            [~,I] = min(Distances);
            %assign centroid to Track
            centroids_ = track(:,:,I);
            
            temp = [centroids(1,i);centroids(2,i)];
            centroids_ = [centroids_,temp]
            
            track_(:,:,I) = centroids_;
        end 
    end
    
    
           
end
