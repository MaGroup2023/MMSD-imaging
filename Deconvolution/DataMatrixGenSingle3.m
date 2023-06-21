
%For a single pixel, only parent or fragment

function IonMatrix =  DataMatrixGenSingle3(RawData,DTScanNumber,IonList,dm)

%Initialize parent matrix and fragment matrix
Ni = length(IonList); %Number of target ions
IonMatrix = zeros(DTScanNumber,Ni);
for DTIndex = 1:DTScanNumber
    CurrentMassSpectra = RawData{DTIndex};
    PointNumber = size(CurrentMassSpectra,1); %Number of data points in this scan
    pl = 1;
    ph = PointNumber;
    if ph==0
        continue;
    end
    for i = 1:Ni % For each target mass
        TargetMass = IonList(i);
        pl_new = SearchFirstBinary(CurrentMassSpectra(pl:ph,1),TargetMass-dm)+pl-1;
        ph_new = SearchLastBinary(CurrentMassSpectra(pl:ph,1),TargetMass+dm)+pl-1;
        if (pl_new>=pl)&&(pl_new<=ph_new)&&(ph_new<=ph) % Search find
            IonMatrix(DTIndex,i) = sum(CurrentMassSpectra(pl_new:ph_new,2));
            pl = pl_new;
        end
    end
end
end

% Search first element e in dp such that e>=k
function p = SearchFirstBinary(dp,k)
    l = 1;
    h = length(dp);
    while l<h
        mid = floor((h+l)/2);
        if dp(mid)>=k
            h = mid;
        else
            l = mid+1;
        end
    end
    if dp(l)>=k
        p=l;
    else
        p=0;
    end
end

% Search last element e in dp such that e<=k
function p = SearchLastBinary(dp,k)
    l = 1;
    h = length(dp);
    while l<h
        mid = ceil((h+l)/2);
        if dp(mid)>k
            h = mid-1;
        else
            l = mid;
        end
    end
    if dp(l)<=k
        p=l;
    else
        p=0;
    end
end




