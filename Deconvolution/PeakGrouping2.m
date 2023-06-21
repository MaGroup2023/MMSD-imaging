
function Groups = PeakGrouping2(DistanceScore,SimilarityScore,DistanceThreshold,SimilarityThreshold)
Np = size(DistanceScore,1);
Groups = cell(Np,1);
for i = 1:Np   
    DistanceList = DistanceScore(i,:);
    SimilarityList = SimilarityScore(i,:);
    Groups{i} = find(DistanceList<DistanceThreshold & SimilarityList>SimilarityThreshold);
end
end