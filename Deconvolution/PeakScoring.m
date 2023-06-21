
function [DistanceScore, SimilarityScore] = PeakScoring(ParentMatrix)
[Nt,Np] = size(ParentMatrix);
DistanceScore = zeros(Np,Np); %unper: peak distance
SimilarityScore = ones(Np,Np); %uper: peak shape similarity
for i = 1:Np
    for j = 1:Np
        ParentCorr = xcorr(ParentMatrix(:,i),ParentMatrix(:,j))/(norm(ParentMatrix(:,i))*norm(ParentMatrix(:,j)));
        [Similarity,Distance] = max(ParentCorr);
        DistanceScore(i,j) = abs(Distance-Nt);
        SimilarityScore(i,j) = Similarity;
    end
end
end