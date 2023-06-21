clear;clc;close all;
%%
load ParentIonList.mat;

%%
LipidDB{1} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PA-H",VariableNamingRule="preserve");
LipidDB{2} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PS-H",VariableNamingRule="preserve");
LipidDB{3} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PE-H",VariableNamingRule="preserve");
LipidDB{4} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PEO-H",VariableNamingRule="preserve");
LipidDB{5} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PG-H",VariableNamingRule="preserve");
LipidDB{6} = readtable("LipidFragmentDatabaseReduced.xlsx",Sheet="PI-H",VariableNamingRule="preserve");
%% Calibration
CaliPara = [0,1.00003438357182,-0.00434807785177298];
ParentIonList = CaliPara(1)*(ParentIonList.^2)+CaliPara(2)*ParentIonList+CaliPara(3);

%%
ID = cell(size(ParentIonList));
MassError = cell(size(ParentIonList));
for i = 1:length(ParentIonList)
    [ID{i},MassError{i}] = LipidMS1SearchSingle(ParentIonList(i),LipidDB,0.008);
end
Anno = cell(size(ParentIonList));
LipidTYPE = {'PA ','PS ','PE ','PEO ','PG ','PI '};
TepO = '';
TepN = '';
for i = 1:length(ParentIonList)
    MassError{i} = unique(MassError{i});
    A = ID{i};
    K = size(A,1);
%     Anno{i} = [num2str(i),' : ',num2str(ParentIonList(i))];
    for k=1:K
        LDB = LipidDB{A(k,1)};
        SN = LDB{A(k,2),4:7};
        TepN = [LipidTYPE{A(k,1)},num2str(SN(1)+SN(3)),':',num2str(SN(2)+SN(4))];
        if ~strcmp(TepN,TepO)
            Anno{i} = [Anno{i},' ',TepN];
        end
        TepO = TepN;
    end
end