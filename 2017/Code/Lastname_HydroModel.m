function [..] = Lastname_HydroModel(day_month,P,n,z,dt,N_years,doTest,...)

%preallocation of variables
Q = zeros(size(P));    %Discharge
R = zeros(size(P));    %Run off
I = zeros(size(P));    %Infiltration
L = zeros(size(P));    %Leaching
ET = zeros(size(P));   %Evapotranspiration
qsup = zeros(size(P)); %Superficial specific discharge
qsub = zeros(size(P)); %Sub-superficial specific discharge
s = zeros(size(P));    %Soil moisture
Vsup = zeros(size(P)); %Superficial storage
Vsub = zeros(size(P)); %Sub-superficial storage
Qsup = zeros(size(P)); %Superficial discharge
Qsub = zeros(size(P)); %Sub-superficial discharge


s(1)=0.5;  %initial conditions

t=0;
for y=1:N_years                     %for loop on the number of years
    for m=1:12                      %for loop on the number of month
        for ts=1:day_month(m)*24;   %for loop on the number of hours of month m              
            t=t+1;                  %counter of hours
             
            % ...to be completed
            % Integrate the system with the euler explicit method:
            % first compute all the water fluxes at time t, then update the
            % water storages. 
            
            
        end
    end
end

%compute total discharge Q    
if doTest     % if doTest=1, check mass balance 

    P_tot=sum(P)*dt;
    R_tot=sum(R)*dt;
    L_tot=sum(L)*dt;
    ET_tot=sum(ET)*dt;

    %"testS" balance for the root zone (input/output). testS close to unity: necessary (but not sufficient) 
    %condition for the implementation to be correct 
    testS=P_tot/(ET_tot+R_tot+L_tot+n*z*(s(end)-s(1)))

    %"testQ" balance for the whole system (input/output). testQ close to unity: necessary (but not sufficient) 
    %condition for the implementation to be correct 
    testQ=sum(P-ET)/(sum(qsup+qsub)+n*z*(s(end)-s(1))+Vsup(end)+Vsub(end))
end
end

