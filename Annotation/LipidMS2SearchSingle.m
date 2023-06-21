

function [MatchID,MatchScore,MatchFragmentsRef,MatchFragmentsInput] = LipidMS2SearchSingle(ParentMass,FragmentMass,FragmentIntensity,Database,options)
    
arguments
    ParentMass (1,1) double {mustBePositive}
    FragmentMass (1,:) double {mustBePositive}
    FragmentIntensity (1,:) double {mustBeNonnegative}
    Database (1,:) cell
    options.ParentMassTolerance (1,1) double {mustBePositive} = 0.008 % Accurate parent ion mass match tolerance
    options.FragmentMassTolerance (1,1) double {mustBePositive} = 0.01 % Fragment ion mass match tolerance
    options.FragmentIntensityThreshold (1,1) double {mustBePositive} = 0.02 % Fragment intensity threshold
    options.MS2ScoreThreshold (1,1) double {mustBePositive} = 0.4 % MS/MS score threshold
    options.FragmentNumTolerance (1,1) {mustBeInteger,mustBePositive} = 3;
end

MatchFragmentsInput = [];
FragmentMassReduced = FragmentMass(FragmentIntensity>(options.FragmentIntensityThreshold));

Nc = length(Database); % Number of lipid classes
MatchID = []; % Lipid match class and id
MatchScore = []; % Lipid fragment match proportion and number
MatchFragmentsRef = {};

if isempty(FragmentMassReduced)
    return
end

for ClassId = 1:Nc
    DB = Database{ClassId}; % Lipid class database
    % Step 1: Search accurate mass
    DBParentMass = DB{:,2}; % Database parent ion m/z
    MS1MatchID = find(abs(DBParentMass-ParentMass)<options.ParentMassTolerance); % MS1 search
    MS1MatchNum = length(MS1MatchID); % Number of MS1 match
    if MS1MatchNum==0 % If no MS1 match
        continue;
    end
    % Step 2: Search MS/MS spectrum
    for i = 1:MS1MatchNum
        MS2Reference = DB{MS1MatchID(i),8:end}; % Reference MS/MS m/z
        MS2Reference(MS2Reference==0) = [];
        if length(FragmentMassReduced)==1 || length(MS2Reference)==1
            FragmentMatchMap = abs(FragmentMassReduced'-MS2Reference)>options.FragmentMassTolerance;
        else
            FragmentMatchMap = cumprod(abs(FragmentMassReduced'-MS2Reference)>options.FragmentMassTolerance);
        end
        FragmentMatchMap = FragmentMatchMap(end,:);
        FragmentMatchNum = length(MS2Reference)-sum(FragmentMatchMap); % Fragment match number
        FragmentMatchProportion = FragmentMatchNum/length(MS2Reference); % Fragment match proportion
        if FragmentMatchProportion>options.MS2ScoreThreshold || FragmentMatchNum >= options.FragmentNumTolerance
            MatchID = [MatchID;ClassId,MS1MatchID(i)];
            MatchScore = [MatchScore;FragmentMatchProportion,FragmentMatchNum];
            MatchFragmentsRef{end+1} = MS2Reference(~FragmentMatchMap);
        end
    end
end
% Sort
if size(MatchID,1)>0
    OveralScore = MatchScore(:,1).*MatchScore(:,2);
    [~,SortId] = sort(OveralScore,'descend');
    MatchID = MatchID(SortId,:);
    MatchScore = MatchScore(SortId,:);
    MatchFragmentsRef = MatchFragmentsRef(SortId);

    %
    for i = 1:size(MatchID,1)
        for j = 1:length(MatchFragmentsRef{i})
            F = MatchFragmentsRef{i};
            MatchFragmentsInput = [MatchFragmentsInput,FragmentMass(abs(FragmentMass-F(j))<options.FragmentMassTolerance)];
        end
    end
end




end