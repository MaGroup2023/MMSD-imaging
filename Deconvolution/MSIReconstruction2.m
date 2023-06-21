 
% Reconstruct IM-MS/MS imaging data from mzML

clear;clc;close all;
DTScanNumber = 200; %Number of scans in drift time dimension
%%
ParentIonList = importdata("PixelDecFeature\ParentIonListMSI.txt");
FragmentIonList = importdata("PixelDecFeature\FragmentIonListMSI.txt");
load ParentGroupsForMSI.mat;
Np = length(ParentIonList);
Nf = length(FragmentIonList);
PixelParentMatrixLE = zeros(DTScanNumber,Np); 
PixelParentMatrixHE = zeros(DTScanNumber,Np);
PixelFragmentMatrix = zeros(DTScanNumber,Nf);

%% Step3: Defining features to be shown
FragmentFeatures=importdata("PixelDecFeature/FragmentationFeature.txt",' ',Np);
FragmentFeatureList = [];
for i = 1:Np
    Fragments = str2num(FragmentFeatures{i});
    for j = 1:length(Fragments)
        FragmentFeatureList = [FragmentFeatureList;[i,Fragments(j)]];
    end
end
Nff = length(FragmentFeatureList);
Rows = 263;
Lines = 228;
NpixelTotal = Rows*Lines;
NpixelLE = ceil(NpixelTotal/2);
NpixelHE = floor(NpixelTotal/2);
ParentImageLE = zeros(NpixelTotal,Np);
TICLE = ones(NpixelTotal,1);
ParentImageHE = zeros(NpixelTotal,Np);
TICHE = ones(NpixelTotal,1);
FragmentSumImage = zeros(NpixelTotal,Nf);
TICHE2 = ones(NpixelTotal,1);
FragmentImage = zeros(NpixelTotal,Nff);
FragmentImage_k = zeros(NpixelTotal,Nff);

%%
load("ClassID.mat");
ROIPosition = find(ClassID>0);
nROI = length(ROIPosition);
dm = 0.05/2; %Peak detection tolerance
FilePathHE = "../DataConversion/Output_mat_DataFile/HEPixelData";
FilePathLE = "../DataConversion/Output_mat_DataFile/LEPixelData";
tic
for k = 1:nROI % Calculating parent features from low energy data and high energy data
    i = k2kk([Rows,Lines],ROIPosition(k));
    if mod(i,2)==1 %LE pixel
        PixelId = ceil(i/2);
        FileNameLE = FilePathLE+num2str(PixelId,'%06d')+".mat";
        PixelDataLE = load(FileNameLE);
        PixelDataLE = PixelDataLE.PixelData;
        if ~isempty(PixelDataLE)
                [ParentImageLE(i,:),TICLE(i)] = PixelSpectraGen(PixelDataLE,ParentIonList,dm);
        end
    else
        PixelId = i/2;
        FileNameHE = FilePathHE+num2str(PixelId,'%06d')+".mat";
        PixelDataHE = load(FileNameHE);
        PixelDataHE = PixelDataHE.PixelData;
        if ~isempty(PixelDataHE)
                [ParentImageHE(i,:),TICHE(i)] = PixelSpectraGen(PixelDataHE,ParentIonList,dm);
                [FragmentSumImage(i,:),TICHE2(i)] = PixelSpectraGen(PixelDataHE,FragmentIonList,dm);
        end  
    end
end
toc
%%

figure(99)
subplot(121)
TLE = sum(ParentImageLE,2);
TLE = 100*TLE/max(TLE);
TLE = reshape(TLE,Lines,[]);
TLE = TLE';
imshow(TLE,[],'InitialMagnification','fit');colormap hot;
THE = sum(ParentImageHE,2);
THE = 100*THE/max(THE);
THE = reshape(THE,Lines,[]);
THE = THE';
subplot(122)
imshow(THE,[],'InitialMagnification','fit');colormap hot;
%%
idLE = 10;
f1 = figure(1);
f1.Name = ['LE precursor ion:', num2str(ParentIonList(idLE))];
subplot(211)
PLE = ParentImageLE(:,idLE)./TICLE;
PLE = 100*PLE/max(PLE);
PhistLE = PLE(PLE>0);
histogram(PhistLE,100)
hold on;
histfit(PhistLE,100,'kernel')
hold off;
set(gca,'TickDir','out')
xlim([0,max(PLE)]);
xticks(0:10:100);
subplot(212)
PLE = reshape(PLE,Lines,[]);
PLE = PLE';
PLE = ImageRangeAdjust(PLE,"LowerLimit",0.5,"UpperLimit",99.5);
rgLE = [];
imshow(PLE,rgLE,'InitialMagnification','fit');colormap hot;colorbar;

%% Step 8: Calculating fragment features
DTStart = 1;
DTEnd = 200;
tic
for k = 1:nROI
    if k==100
        toc        
         k;
    end
    kk = k2kk([Rows,Lines],ROIPosition(k));
    if mod(kk,2)==1
        i = (kk+1)/2;
        FileNameLE = FilePathLE+num2str(i,'%06d')+".mat";
        FileNameHE = FilePathHE+num2str(i,'%06d')+".mat";
        PixelDataLE = load(FileNameLE);
        PixelDataLE = PixelDataLE.PixelData;
        PixelDataHE = load(FileNameHE);
        PixelDataHE = PixelDataHE.PixelData;
        if isempty(PixelDataLE) || isempty(PixelDataHE)
                continue;
        end
        ParentMatrix = DataMatrixGenSingle3(PixelDataLE,DTScanNumber,ParentIonList,dm);
        FragmentMatrix = DataMatrixGenSingle3(PixelDataHE,DTScanNumber,FragmentIonList,dm);
        FragmentMatrix = circshift(FragmentMatrix,4,1);
        ParentSpectra = sum(ParentMatrix,1);
        C = zeros(Np,Nf);
        for j=1:Nf
            C(:,j) = lsqnonneg(ParentMatrix,FragmentMatrix(:,j));
        end
        CC = zeros(Np,Nf);
        CC_k = zeros(Np,Nf);
        for j = 1:Np
            PeakIndex = ParentGroups{j};
            CC(j,:) = ParentSpectra(PeakIndex)*C(PeakIndex,:);
            CC_k(j,:) = ParentSpectra(PeakIndex)*C(PeakIndex,:)/sum(ParentSpectra(PeakIndex));
        end
        for j = 1:Nff
            FragmentImage(kk,j) = CC(FragmentFeatureList(j,1),FragmentFeatureList(j,2));
            FragmentImage_k(kk,j) = CC_k(FragmentFeatureList(j,1),FragmentFeatureList(j,2));
        end
    end
end
toc

%% Step 9: Show fragment features
fx = figure(777);
id = 38;
fx.Name = [num2str(ParentIonList(FragmentFeatureList(id,1))),'-->',num2str(FragmentIonList(FragmentFeatureList(id,2)))];
subplot(211)
FHEk = FragmentImage_k(:,id)./TICHE;
FHEk = 100*FHEk/max(FHEk);
Fhist = FHEk(FHEk>0.1 & FHEk<20);
histogram(Fhist)
set(gca,'TickDir','out');
xticks(0:10:100);
subplot(212)
FHEk = reshape(FHEk,Lines,[]);
FHEk = FHEk';
FHEk(:,2:2:end)=FHEk(:,1:2:end);
FHEk = ImageRangeAdjust(FHEk,"LowerLimit",0.5,"UpperLimit",99.5);
imshow(FHEk,[],'InitialMagnification','fit');colormap hot;
%%
f3 = figure(3);
id = 38;
f3.Name = [num2str(ParentIonList(FragmentFeatureList(id,1))),'-->',num2str(FragmentIonList(FragmentFeatureList(id,2)))];
subplot(211)
FHE = FragmentImage(:,id)./TICHE;
FHE = 100*FHE/max(FHE);
Fhist = FHE(FHE>1);
histogram(Fhist)
xlim([0,max(FHE)]);
xticks(0:10:100);
subplot(212)
FHE = reshape(FHE,Lines,[]);
FHE = FHE';
FHE(:,2:2:end)=FHE(:,1:2:end);
FHE = ImageRangeAdjust(FHE,"LowerLimit",.5,"UpperLimit",99.5);
rgF = [];
imshow(FHE,rgF,'InitialMagnification','fit');colormap hot;


idFSum = FragmentFeatureList(id,2);
f77 = figure(77);
f77.Name = ['HE fragment ion:', num2str(FragmentIonList(idFSum))];
subplot(211)
FHES = FragmentSumImage(:,idFSum)./TICHE2;
FHES = 100*FHES/max(FHES);
FhistHES = FHES(FHES>0);
histogram(FhistHES)
xlim([0,max(FHES)]);
xticks(0:10:100);
subplot(212)
FHES = reshape(FHES,Lines,[]);
FHES = FHES';
FHES = ImageRangeAdjust(FHES,"LowerLimit",.5,"UpperLimit",99.5);
rgFHES = [];
imshow(FHES,rgFHES,'InitialMagnification','fit');colormap hot;colorbar;
%% show
f4 = figure(4);
f4.Name = [num2str(idLE),': ', num2str(ParentIonList(idLE))];
PLE(:,2:2:end)=PLE(:,1:2:end);
imshow(PLE,[0,100],'InitialMagnification','fit');
colormap('hot');
set(gcf,'Position',[100,200,520,450]);
set(gca,'Position',[0,0,1,1]);
f5 = figure(5);
f5.Name = [num2str(ParentIonList(FragmentFeatureList(id,1))),'-->',num2str(FragmentIonList(FragmentFeatureList(id,2)))];
FHE = ImageRangeAdjust(FHE,"LowerLimit",.5,"UpperLimit",99.5);
imshow(FHE,[],'InitialMagnification','fit');
colormap('hot');
set(gcf,'Position',[500,300,520,450]);
set(gca,'Position',[0,0,1,1]);
f18 = figure(18);
f18.Name = [num2str(FragmentFeatureList(id,2)), ': ', num2str(FragmentIonList(FragmentFeatureList(id,2)))];
FHES(:,1:2:end)=FHES(:,2:2:end);
imshow(FHES,[0,100],'InitialMagnification','fit');
colormap('hot');
set(gcf,'Position',[950,400,520,450]);
set(gca,'Position',[0,0,1,1]);
