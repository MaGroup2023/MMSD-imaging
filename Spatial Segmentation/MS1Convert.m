% Convert a MS1 imaging file (output *.txt file from HDI software) to a 3D
% array
% mz : m/z of extract ions
% ionMaps(r,c,n) : the n th ion intensity at pixel location (r,c)
% TIC : TIC intensity of each pixel

function [mz,ionMapsLE,ionMapsHE,TICLE,TICHE] = MS1Convert(filename,R,C1,C2)
    raw = readtable(filename);
    raw = table2array(raw);
    mz = raw(3,4:end-2);
    N1 = ceil(R/2)*C1+floor(R/2)*C2;
    N2 = ceil(R/2)*C2+floor(R/2)*C1;
    N = R*(C1+C2);
    raw = raw(4:3+N,4:end-2);
    ionMapsLE = zeros(size(raw));
    TICLE = ones(size(raw,1),1);
    ionMapsLE(1:2:end,:) = raw(1:N1,:);
    ionMapsLE(2:2:end,:) = ionMapsLE(1:2:end,:);%
    TICLE(1:2:end,:) = sum(raw(1:N1,:),2);
    TICLE(2:2:end,:) = TICLE(1:2:end,:);%
    ionMapsHE = zeros(size(raw));    
    TICHE = ones(size(raw,1),1);
    ionMapsHE(2:2:end,:) = raw(N1+1:N,:);
    ionMapsHE(1:2:end,:) = ionMapsHE(2:2:end,:);%
    TICHE(2:2:end,:) = sum(raw(N1+1:N,:),2);
    TICHE(1:2:end,:) = TICHE(2:2:end,:);%
    ionMapsLE = reshape(ionMapsLE,C1+C2,R,[]);
    ionMapsLE = permute(ionMapsLE,[2 1 3]);
    ionMapsHE = reshape(ionMapsHE,C1+C2,R,[]);
    ionMapsHE = permute(ionMapsHE,[2 1 3]);
    TICLE = reshape(TICLE,C1+C2,R,[]);
    TICLE = permute(TICLE,[2,1,3]);
    TICHE = reshape(TICHE,C1+C2,R,[]);
    TICHE = permute(TICHE,[2,1,3]);
end