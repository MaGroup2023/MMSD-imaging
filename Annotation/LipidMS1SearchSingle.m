function [MatchID,MassError] = LipidMS1SearchSingle(ParentMass,Database,MassTolerance)

arguments
    ParentMass (1,1) double {mustBePositive}
    Database cell
    MassTolerance (1,1) double {mustBePositive} = 0.008
end

Nc = length(Database); % Number of lipid classes
MatchID = [];
MassError = [];
for ClassID = 1:Nc % For each lipid class
    DB = Database{ClassID}; % Lipid class database

    DBParentMass = DB{:,2}; % Database parent ion m/z
    MS1MatchID = find(abs(DBParentMass-ParentMass)<MassTolerance); % MS1 search
    MS1MatchNum = length(MS1MatchID); % Number of MS1 match
    for i = 1:MS1MatchNum
        MatchID(end+1,:) = [ClassID,MS1MatchID(i)];
        MassError(end+1) = DBParentMass(MS1MatchID(i))-ParentMass;
    end
end

[~,SID] = sort(abs(MassError),'ascend');
MatchID = MatchID(SID,:);
MassError = MassError(SID);

end