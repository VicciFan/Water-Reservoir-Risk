function Phourly=downscaling(Pdaily)

%downscale daily precipitation to hourly precipitation assuming that hourly
%rainfall are exponentially distributed.

%INPUT
%"Pdaily": precipitation at daily timestep

%OUTPUT
%"Phourly": precipitation at horly time step 

%if "Pdaily" is in [mm/day], "Phourly" is in [mm/h] 

Phourly=zeros(length(Pdaily)*24,1);

for i=1:length(Pdaily)
    distr=-log(rand(24,1));
    Phourly((i-1)*24+1:i*24)=Pdaily(i)*distr./sum(distr);
end