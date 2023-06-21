function [ClassID,CentroidSpectrum,SAD,iter,loss] = SAKMeans(S,K,Mask,r,alpha,MaxIter)
% ClassID: class id for each sample; CentroidPos: k centroid positions
% S: Data matrix; sz: size of image; k: Number of clusters; r: SA range; alpha: SA coeff
[N1,N2,M] = size(S); % Number of samples (N1*N2) and features (M)
ClassID = randi(K,N1,N2); % Initial class id
ClassID(~Mask) = 0;
CentroidOld = zeros(K,M); % Centroid positions before update
SAD = zeros(N1,N2,K);
iter = 1;
while iter<MaxIter
    iter = iter+1;
    % Step A: Calculate new centroid
    CentroidSpectrum = zeros(K,M);
    for k = 1:K
        [rows, cols] = find(ClassID==k);
        nk = length(rows);
        for i = 1:nk
            CentroidSpectrum(k,:) = CentroidSpectrum(k,:) + reshape(S(rows(i),cols(i),:)/nk,1,[]);
        end
    end
    loss(iter-1) =  norm(CentroidSpectrum-CentroidOld,"fro")/norm(CentroidSpectrum,"fro");
    if loss(iter-1)<0.00001
        iter = iter-1;
        break
    else
        CentroidOld = CentroidSpectrum;
    end
    % Step B: Assign pixels to clusters
    for sample_i = 1:N1
        for sample_j = 1:N2 % For each sample
            if ~Mask(sample_i,sample_j)
                continue;
            end
            % Step 1: Find adjacent sample pixels
            X1 = zeros(2*r+1,2*r+1,M); % Adjacent sample data
            for di = -r:r
                for dj = -r:r
                    adjacent_i = sample_i+di;
                    adjacent_j = sample_j+dj;
                    if adjacent_i>0 && adjacent_i<=N1 && adjacent_j>0 && adjacent_j<=N2
                        X1(di+r+1,dj+r+1,:) = S(adjacent_i,adjacent_j,:);
                    end
                end
            end
            % Step 2: Calculate SA distance
            SADistance = zeros(K,1);
            for k = 1:K % For each centroid
                X2 = zeros(2*r+1,2*r+1,M);
                for xi = 1:2*r+1
                    for xj = 1:2*r+1
                        X2(xi,xj,:) = CentroidSpectrum(k,:);
                    end
                end
                Distance = sum((X1-X2).^2,3);
                SADistance(k) = sum(Distance.*alpha,'all');
            end
            SAD(sample_i,sample_j,:) = SADistance;
            [~,ClusterID] = min(SADistance);
            ClassID(sample_i,sample_j) = ClusterID;
        end
    end
    

end


end