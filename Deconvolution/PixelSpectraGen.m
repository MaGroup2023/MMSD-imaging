

function [PixelSpectra,TIC] = PixelSpectraGen(PixelData,IonList,dm)
Ni = length(IonList); %Number of ions
PixelSpectra = zeros(1,Ni);
DTScanNumber = length(PixelData);
TIC = 0;
for DTIndex = 1:DTScanNumber
    CurrentSpectra = PixelData{DTIndex};
    PointNumber = size(CurrentSpectra,1);
    TIC = TIC + sum(CurrentSpectra(:,2));
    for MassIndex = 1:PointNumber
        CurrentMass = CurrentSpectra(MassIndex,1);
        for TargetMassIndex = 1:Ni %For each target ion
            CurrentTargetMass = IonList(TargetMassIndex);
            if CurrentMass > CurrentTargetMass-dm
                if CurrentMass <CurrentTargetMass +dm
                    PixelSpectra(TargetMassIndex) = PixelSpectra(TargetMassIndex)+CurrentSpectra(MassIndex,2);
                    break;
                end
            else
                break;
            end
        end
    end
end