%Calcultating pixel average
clear;clc;close all;
FilePathLE = "../DataConversion/Output_mat_DataFile/LEPixelData";
FilePathHE = "../DataConversion/Output_mat_DataFile/HEPixelData";
Rows = 263; %Number of rows
Lines = 114; %Number of lines
DTScanNumber = 200; 
dm = 0.04/2;
RawParent = importdata("LE_MS.txt");
RawFragment = importdata("HE_MS.txt");
ILG = IonListGen;
ParentIonList = ILG.Relative(RawParent,700,900,0.01);
FragmentIonList = ILG.Relative(RawFragment,100,900,0.005);
% ParentIonList = importdata("PixelDecFeature\ParentIonListMSI.txt");
% FragmentIonList = importdata("PixelDecFeature\FragmentIonListMSI.txt");
Np = length(ParentIonList);
Nf = length(FragmentIonList);
AverageParentMatrix = zeros(DTScanNumber,Np);
AverageFragmentMatrix = zeros(DTScanNumber,Nf);
%%
% load mask
load("ClassID.mat");

ROIPosition = find(ClassID>0);
sz = size(ClassID);
%%
tic
for p = 1:length(ROIPosition)
    k = ROIPosition(p);
    kk = k2kk(sz,k);
    FileId = ceil(kk/2);
    if mod(kk,2)==1 % LE data
        FileNameLE = FilePathLE+num2str(FileId,'%06d')+".mat";
        load(FileNameLE);
        if isempty(PixelData)
            continue;
        end
        P = DataMatrixGenSingle3(PixelData,200,ParentIonList,dm);      
        AverageParentMatrix = AverageParentMatrix+P;
    else % HE data
        FileNameHE = FilePathHE+num2str(FileId,'%06d')+".mat";
        load(FileNameHE);
        if isempty(PixelData)
            continue;
        end
        F = DataMatrixGenSingle3(PixelData,200,FragmentIonList,dm);      
        AverageFragmentMatrix = AverageFragmentMatrix+F;
    end
end
toc
AverageParentMatrix = AverageParentMatrix/(Rows*Lines);
AverageFragmentMatrix = AverageFragmentMatrix/(Rows*Lines);
%%
DTStart = 1;
DTEnd = 200;
ParentMatrix = AverageParentMatrix(DTStart:DTEnd,:);
FragmentMatrix = AverageFragmentMatrix(DTStart:DTEnd,:);
FragmentMatrix = circshift(FragmentMatrix,4,1);
%%

ParentSpectra = sum(ParentMatrix,1);

%%
[Nt,Np] = size(ParentMatrix);
Nf = size(FragmentMatrix,2);
%%
%Calculting fragmentation coefficient by nonnegtive least square method
tic
C = zeros(Np,Nf);
for i=1:Nf
    C(:,i) = lsqnonneg(ParentMatrix,FragmentMatrix(:,i));
end
toc

%% Scoring
[DistanceScore,SimilarityScore] = PeakScoring(ParentMatrix);
%% Grouping
ParentGroups = PeakGrouping2(DistanceScore,SimilarityScore,2,0.98);
%%
%Combining fragmentation matrix
CombinedFragmentationMap = zeros(Np,Nf);
for i = 1:Np
    PeakIndex = ParentGroups{i};
    CombinedFragmentationMap(i,:) = ParentSpectra(PeakIndex)*C(PeakIndex,:)/sum(ParentSpectra(PeakIndex));
end
save CombinedFragmentationMap.mat CombinedFragmentationMap
%%
%disp
mass = FragmentIonList;
figure(1)
id = 79; % REPLACE.
I = CombinedFragmentationMap(id,:);
I = 100*I/(max(I));
% idxmag = mass>200 & mass<600;
% I(idxmag) = 5*I(idxmag);
msms = [mass,I'];
bar(mass,I,50,'k');
xlim([100,900]);
ylim([0,100]);
xticks(200:200:900);
yticks(0:50:100);
% set(gcf,'Position',[200,200,600,200])
box off;
set(gcf,'Units','centimeters','Position',[20,15,15,5]);
set(gca,'Position',[0.13,0.17,0.8,0.74]);
set(gca,'TickDir','out');
set(gca,'FontName','Arial','FontSize',7);
%%
save ParentIonList.mat ParentIonList
save FragmentIonList.mat FragmentIonList
save ParentGroupsForMSI.mat ParentGroups
save CombinedFragmentationMap.mat CombinedFragmentationMap




