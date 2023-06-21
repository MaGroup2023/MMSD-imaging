% Given a m/z-intensity list of parent/fragment ions, this three functions
% generate a reduced, centroid m/z list, which is further used to extracte the parent
% and fragment matrixes from the original 2D IM-MS/MS raw data. The first
% and second functions generates ion list based on absolute and relative
% threshold, the third function finds the top N (intensity) ions, and the
% last function , based on relative intensity, takes in an exception list,
% which is special for fragment ion list

function ILG = IonListGen
    ILG.Absolute = @IonListAbsoluteGen;
    ILG.Relative = @IonListRelativeGen;
    ILG.Top = @IonListTopGen;
    ILG.Except = @IonListExcept;
end

%Function 1: IonListAbsoluteGen
%Input items:
%IonMassSpectra : Profile mode of the mass spectra of the targeted ions.
%LowMass/HighMass : Mass range of the targeted ions to be determined by the
%user.
%IntensityThreshold : The threshold of the intensity (Absolute) of targeted
%ions.
function IonList = IonListAbsoluteGen(IonMassSpectra, LowMass, HighMass, IntensityThreshold)
    [~,MassId] = findpeaks(IonMassSpectra(:,2));
    PeakList = IonMassSpectra(MassId,:);
    TotalMassRange = PeakList(:,1); %Mass range of input spectra
    TargetId = TotalMassRange>LowMass & TotalMassRange<HighMass; %Mass range limited by LowMass and HighMass
    TargetMass = PeakList(TargetId,1);
    TargetIntensity = PeakList(TargetId,2);
    ResultId = TargetIntensity > IntensityThreshold;
    IonList = TargetMass(ResultId);
end

%Function 2: IonListRelativeGen
%Input items:
%IonMassSpectra : Profile mode of the mass spectra of the targeted ions.
%LowMass/HighMass : Mass range of the targeted ions to be determined by the
%user.
%IntensityThreshold : The threshold of the intensity (Relative) of targeted
%ions.

function IonList = IonListRelativeGen(IonMassSpectra, LowMass, HighMass, IntensityThreshold)
    [~,MassId] = findpeaks(IonMassSpectra(:,2));
    PeakList = IonMassSpectra(MassId,:);
    TotalMassRange = PeakList(:,1); %Mass range of input spectra
    TargetId = TotalMassRange>LowMass & TotalMassRange<HighMass; %Mass range limited by LowMass and HighMass
    TargetMass = PeakList(TargetId,1);
    TargetIntensity = PeakList(TargetId,2);
    TargetIntensity = TargetIntensity/max(TargetIntensity);
    ResultId = TargetIntensity > IntensityThreshold;
    IonList = TargetMass(ResultId);
end

%Function 3: IonListTopGen
%Input items:
%IonMassSpectra : Profile mode of the mass spectra of the target ions.
%LowMass/HighMass : Mass range of the target ions to be determined by the
%user.
%Number : Number of target ions.

function IonList = IonListTopGen(IonMassSpectra, LowMass, HighMass, Number)
    [~,MassId] = findpeaks(IonMassSpectra(:,2));
    PeakList = IonMassSpectra(MassId,:);
    TotalMassRange = PeakList(:,1); %Mass range of input spectra
    TargetId = TotalMassRange>LowMass & TotalMassRange<HighMass; %Mass range limited by LowMass and HighMass
    TargetMass = PeakList(TargetId,1);
    TargetIntensity = PeakList(TargetId,2);
    [~, SortId] = sort(TargetIntensity,'descend');
    IonList = TargetMass(SortId(1:Number));
    IonList = sort(IonList);
end

%Function 4: IonListExcept
%Input items:
%IonMassSpectra : Profile mode of the mass spectra of the target ions.
%NegList : An exception list.
%Tolerance : Mass tolerance of the exception list.
%LowMass/HighMass : User defined mass range of the target ions.
%IntensityThreshold : Relative intensity threshold.

function IonList = IonListExcept(IonMassSpectra,NegList,Tolerance, LowMass,HighMass,IntensityThreshold)
    [~,MassId] = findpeaks(IonMassSpectra(:,2));
    PeakList = IonMassSpectra(MassId,:);
    TotalMassRange = PeakList(:,1); %Mass range of input spectra
    TargetId = TotalMassRange>LowMass & TotalMassRange<HighMass; %Mass range limited by LowMass and HighMass
    TargetMass = PeakList(TargetId,1);
    TargetIntensity = PeakList(TargetId,2);
    TargetIntensity = TargetIntensity/max(TargetIntensity);
    ResultId = TargetIntensity > IntensityThreshold;
    IonList = TargetMass(ResultId);
    for i = 1:length(NegList)
        IonList(abs(IonList-NegList(i))<Tolerance) = [];
    end
end
