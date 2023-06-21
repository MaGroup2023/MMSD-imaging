clear;close all;clc;

%%
% ParentMass = [716.5228,728.5607,792.5484];
% FragmentMass = [255.2,281.2,331.2,446.3,452.2,464.2,478.2];
% FragmentInt = [0.2 0.3 0.0 0.0 0.15 0 0.1;
%                0.0 0.1 0.0 0.1 0 0.2 0;
%                0 0.5 0.1 0 0 0 0];
load ParentIonList.mat
load FragmentIonList.mat
load CombinedFragmentationMap.mat

%%
LipidDB{1} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PA-H",VariableNamingRule="preserve");
LipidDB{2} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PS-H",VariableNamingRule="preserve");
LipidDB{3} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PE-H",VariableNamingRule="preserve");
LipidDB{4} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PEO-H",VariableNamingRule="preserve");
LipidDB{5} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PG-H",VariableNamingRule="preserve");
LipidDB{6} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PI-H",VariableNamingRule="preserve");
%% Calibration
CaliParaParent = [0,1.00003438357182,-0.00434807785177298];
ParentIonList = CaliParaParent(1)*(ParentIonList.^2)+CaliParaParent(2)*ParentIonList+CaliParaParent(3);
CaliParaFragment = [0,1.00003156985421,-0.00205971211641812];
FragmentIonList = CaliParaFragment(1)*(FragmentIonList.^2)+CaliParaFragment(2)*FragmentIonList+CaliParaFragment(3);

%%
ID = cell(size(ParentIonList));
Score = cell(size(ParentIonList));
MatchFragments = cell(size(ParentIonList));
MatchFragmentsInput = [];
for i = 1:length(ParentIonList)
    [ID{i},Score{i},MatchFragments{i},M] = LipidMS2SearchSingle(ParentIonList(i),FragmentIonList,CombinedFragmentationMap(i,:),LipidDB,...
        "FragmentIntensityThreshold",0.01,"MS2ScoreThreshold",0.5,"FragmentNumTolerance",3,"FragmentMassTolerance",0.005);
    MatchFragmentsInput = unique([MatchFragmentsInput,M]);
end
MatchFragmentsOriginal = ((MatchFragmentsInput-CaliParaFragment(3))/CaliParaFragment(2))';
Anno = cell(size(ParentIonList));
LipidTYPE = {'PA ','PS ','PE ','PEO ','PG ','PI ','ST '};
for i = 1:length(ParentIonList)
    A = ID{i};
    K = size(A,1);
    for k=1:K
        LDB = LipidDB{A(k,1)};
        SN = LDB{A(k,2),4:7};
        Anno{i} = [Anno{i},' & ',LipidTYPE{A(k,1)},num2str(SN(1)),':',num2str(SN(2)),'_',num2str(SN(3)),':',num2str(SN(4))];
    end
end

MatchFragmentID = cell(size(MatchFragments));
MatchFragmentIDMerge = cell(size(MatchFragments));
C = [1.00003156985421,-0.00205971211641812];
for i = 1:length(MatchFragments)
    FID = cell(size(MatchFragments{i}));
    for j = 1:length(MatchFragments{i})
        F = MatchFragments{i}{j};
        F = (F-C(2))/C(1);
        I = zeros(size(F));
        for k = 1:length(F)
            I(k) = find(abs(MatchFragmentsOriginal-F(k))<0.005);
        end
        FID{j} = sort(I,'ascend');
        MatchFragmentIDMerge{i} = [MatchFragmentIDMerge{i},FID{j}];
    end
    MatchFragmentID{i} = FID;
    MatchFragmentIDMerge{i} = unique(MatchFragmentIDMerge{i});
end

