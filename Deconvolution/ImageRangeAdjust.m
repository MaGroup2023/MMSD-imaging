function Iout = ImageRangeAdjust(Iin,options)
arguments
    Iin double
    options.LowerLimit (1,1) double {mustBeNonnegative}=5
    options.UpperLimit (1,1) double {mustBePositive}=95
end
L = prctile(Iin(:),options.LowerLimit);
H = prctile(Iin(:),options.UpperLimit);
Iout = Iin;
Iout(Iin<L)=L;
Iout(Iin>H)=H;

end