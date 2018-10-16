clear all               %clear all the variables allocated in the workspace
close all               %close all the open figures

edit=1;                 %if "edit"=0 no figures
                        %if "edit"=1 show figures

%load input files
load P.txt              %precipitation data

%*************************************************************************
% PARAMETERS
%*************************************************************************

N_years_gen=100;           %number of years of generated rainfall
N_years=length(P)/24/365;  %number of years of data

%*************************************************************************
%                     RAINFALL GENERATION
%*************************************************************************

%upscaling from hourly rainfall  to daily rainfall 
%the resulting vector "dailyP" should have N_years*365 elements with values in [mm/day].

dailyP=sum(reshape(P,24,365*N_years));        %[mm/day]
dailyP=dailyP(:);                             %transform into column vector



%************************ Montly Statistic ********************************

%"P_day_year": matrix with precipitation values organized in days (rows) and years(columns)
P_day_year=reshape(dailyP,365,N_years);

day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; %"day_month": number of days for each month
month_end=cumsum(day_month);                     %"month_end": last day of each month
month_start=month_end-day_month+1;               %"month_start": first day of each month

monthly_mean_P=zeros(12,1);   %monthly mean precipitation (counting also non rainy days)
monthly_std_P=zeros(12,1);    %monthly precipitation standar deviation (counting also non rainy days)
monthly_lambda=zeros(12,1);   %monthly events arrival rate
monthly_alpha=zeros(12,1);    %monthly mean precipitation (counting only the rainy days)

for m=1:12  %for loop on the month
    P_series=P_day_year(month_start(m):month_end(m),:);
    P_series=P_series(:);                                  %transform into column vector
    monthly_mean_P(m)=mean(P_series);                  
    monthly_std_P(m)=std(P_series);
    monthly_lambda(m)=sum(P_series>0)/length(P_series);   %number of rainy days/total number of days, the probability of that month
    monthly_alpha(m)=mean(P_series(P_series>0));          %all in unit of [mm/day]
end


%***************************** Plot ***************************************
if edit
    figure(1)
    subplot(2,2,1)
    h=plot(1:12,monthly_mean_P,'--ro','markerfacecolor','red');
    set(h,'displayname','data')
    ylabel('mean precipitation [mm]','fontsize',14)
    xlabel('month','fontsize',14)
    hold on
    box off
    subplot(2,2,2)
    h=plot(1:12,monthly_std_P,'--ro','markerfacecolor','red');
    set(h,'displayname','data')
    ylabel('precipitation std [mm]','fontsize',14)
    xlabel('month','fontsize',14)
    hold on
    box off
    subplot(2,2,3)
    h=plot(1:12,monthly_lambda,'--ro','markerfacecolor','red');
    set(h,'displayname','data')
    ylabel('events rate [1/d]','fontsize',14)
    xlabel('month','fontsize',14)
    hold on
    box off
    subplot(2,2,4)
    h=plot(1:12,monthly_alpha,'--ro','markerfacecolor','red');
    set(h,'displayname','data')
    ylabel('mean precipitation (rainy days) [mm]','fontsize',14)
    xlabel('month','fontsize',14)
    hold on
    box off
end

%************************ Rainfall Generation ********************************

%"Pgen_day_year": matrix with generated precipitation values organized in days (rows) and years (columns)
dailyPgen_day_year=zeros(365,N_years_gen);

for y=1:N_years_gen  %for loop on the years
    d=0;  %counter for the day of the year, inizialized every years
    for m=1:12    %for loop on the months
        for i=1:day_month(m)   %for loop on the days of the month
            d=d+1;             %update days counter
            
           %generate the value of rainfall for this day: 
            u=rand;
            if  u<=monthly_lambda(m)       % decide if there is or there isn't a rainfall event during this day    
                dailyPgen_day_year(d,y)=-log(1-rand)*monthly_alpha(m);  %if an evant occurs, assign a precipitation [mm/day] extracted (with the inverse
            %transformation method) from an exponential distribution with mean monthly_alpha(m)
            else
                dailyPgen_day_year(d,y)=0;
            end
           
        end
    end
end

%compute statistics of generated rainfall
monthly_mean_P_gen=zeros(12,1);
monthly_std_P_gen=zeros(12,1);
monthly_lambda_gen=zeros(12,1);
monthly_alpha_gen=zeros(12,1);

for m=1:12
    P_series=dailyPgen_day_year(month_start(m):month_end(m),:);
    P_series=P_series(:);
    monthly_mean_P_gen(m)=mean(P_series);
    monthly_std_P_gen(m)=std(P_series);
    monthly_lambda_gen(m)=sum(P_series>0)/length(P_series);
    monthly_alpha_gen(m)=mean(P_series(P_series>0));
end

%*************************** plot *****************************************
if edit
    figure(1)
    subplot(2,2,1)
    h=plot(1:12,monthly_mean_P_gen,'--bo');
    set(h,'displayname','generated')
    hold off
    box off
    legend('show')
    subplot(2,2,2)
    h=plot(1:12,monthly_std_P_gen,'--bo');
    set(h,'displayname','generated');
    hold off
    box off
    legend('show')
    subplot(2,2,3)
    h=plot(1:12,monthly_lambda_gen,'--bo');
    set(h,'displayname','generated')
    hold off
    box off
    legend('show')
    subplot(2,2,4)
    h=plot(1:12,monthly_alpha_gen,'--bo');
    set(h,'displayname','generated')
    hold off
    box off
    legend('show')
end

%*************************************************************************
%                          HYDROLOGIC MODEL
%*************************************************************************

%downscaling of daily rainfall to hourly rainfall 
dailyPgen=dailyPgen_day_year(:);  %transform in a column vector,[mm/day]
Pgen=downscaling(dailyPgen);      %generated rainfall in [mm/h]

%run the calibrated hydrological model with "Pgen" as forcing rainfall
load temperature.txt        %mean montly temperature
load kc.txt                 %mean monthly crop coefficient

Ksat = 1.45*10^(-6);        %saturated hydraulic conductivity
                            %TO BE CALIBRATED, range 10^(-7)--10^(-5)[m/s];
Ksat = Ksat*1000*3600;            

c=7;           
sw=0.2;        %wilting point []
s1=0.45;       %soil moisture above which plants transpire at kc*ETO []
n=0.3;         %porosity []
z=500;         %root zone depth [mm]
BaseFlow=4;    %base flow [m^3/s]
Ksup=1/15;     %superficial discharge coefficient [h^-1]
Ksub=1/36;     %sub-superficial discharge coefficient [h^-1]
Area=3000;     %area of the basin [km^2]
dt=1;          %timestep of integration [h]
lat=46;        %latitude of the basin [degrees]


%*************************************************************************
%               COMPUTE POTENTIAL EVAPOTRANSPIRATION
%*************************************************************************
N_years=length(Pgen)/24/365;  %number of years with generated pricipitation

day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; %"day_month": number of days for each month
month_end=cumsum(day_month);                     %"month_end": last day of each month
month_start=month_end-day_month+1;               %"month_start": first day of each month


D=1:365;  %day of the years
delta=0.409*sin(2*pi*D/365-1.39);
omega_s=acos(-tan(lat*pi/180)*tan(delta));
N_D=24*omega_s/pi;           %number of daylight hours of day D£¬ 365 values

N_m=zeros(12,1);             %mean monthly number of daylight hours
for m=1:12                   %for loop on the month
    N_m(m)=mean(N_D(month_start(m):month_end(m)));
end
 
mean_monthly_T=temperature(:,1);    
I=sum((mean_monthly_T/5).^1.514);             %heat index
a=6.75e-7*I^3-7.71e-5*I^2+1.79e-2*I+0.49;

ET0=16*(N_m/12).*(10*mean_monthly_T/I).^a;    %potential evapotranspiration in mm for month, 12 values
mean_monthly_ET0=ET0./day_month'/24;          %potential evapotranspiration in mm/hour, 12 values

figure
plot(1:12,mean_monthly_ET0,'--ro','markerfacecolor','red')
xlabel('month','FONTSIZE',14)                                  
ylabel('Potential Evapotranspiration [mm/day]','FONTSIZE',14)   
box off



%*************************************************************************
%                          HYDROLOGIC MODEL
%*************************************************************************

%preallocation of variables
R=zeros(size(P));     %Run off
I=zeros(size(P));     %Infiltration
L=zeros(size(P));     %Leaching
ET=zeros(size(P));    %Evapotranspiration
qsup=zeros(size(P));  %Superficial specific discharge
qsub=zeros(size(P));  %Sub-superficial specific discharge
s=zeros(size(P));     %Soil moisture
Vsup=zeros(size(P));  %Superficial storage
Vsub=zeros(size(P));  %Sub-superficial storage


s(1)=0.5;  %initial conditions
t=0;
for y=1:N_years                     %for loop on the number of years
    for m=1:12                      %for loop on the number of month
        for ts=1:day_month(m)*24;   %for loop on the number of hours of month m              
            t=t+1;                  %counter of hours
             
            I(t) = min(Pgen(t),Ksat);
            R(t) = Pgen(t)-I(t);           
            L(t) = Ksat*s(t)^c;
            
            %compute ET(t): the evaporation at t, [mm/hour] 
            if s(t)<=sw
               ET(t)=0;
            else
                if  s(t)<=s1
                    ET(t)=kc(m)*mean_monthly_ET0(m)*((s(t)-sw)/(s1-sw));
                else 
                    ET(t)=kc(m)*mean_monthly_ET0(m);
                end   
            end
            
            %update s(t),Euler explicit method:
            s(t+1)=s(t)+(I(t)-L(t)-ET(t))*dt/(n*z);
             
            qsup(t)=Ksup*Vsup(t);                      %[mm/h]
            Vsup(t+1)=Vsup(t)+(R(t)-qsup(t))*dt;       %[mm]
            qsub(t)=Ksub*Vsub(t);
            Vsub(t+1)=Vsub(t)+(L(t)-qsub(t))*dt;
            Qsup(t)=Area*10^6*(qsup(t)/(3600*1000));   %units:[m^2 * (m/s)]
            Qsub(t)=Area*10^6*(qsub(t)/(3600*1000));
            
              
        end
    end
end


%compute total discharge Q
Q=Qsup+Qsub+BaseFlow;                   %[m^3/s]


time=2000+(0:length(Pgen)-1)/(24*365);  %time points in year
if edit
    figure(2)
    subplot(2,3,1)
    plot(time,Pgen)
    xlabel('year','FONTSIZE',14)
    ylabel('precipitation [mm/h]','FONTSIZE',14)
    box off
    
    subplot(2,3,2)
    plot(time,R)
    xlabel('year','FONTSIZE',14)
    ylabel('runoff [mm/h]','FONTSIZE',14)       
    box off
    
    subplot(2,3,3)
    plot(time,I)
    xlabel('year','FONTSIZE',14)
    ylabel('infiltration [mm/h]','FONTSIZE',14)  
    box off
    
    subplot(2,3,4)
    plot(time,s(1:end-1))
    xlabel('year','FONTSIZE',14)
    ylabel('soil moisture','FONTSIZE',14)       
    box off
    
    subplot(2,3,5)
    plot(time,L)
    xlabel('year','FONTSIZE',14)
    ylabel('Leakage [mm/h]','FONTSIZE',14)       
    box off
    
    subplot(2,3,6)
    plot(time,ET)
    xlabel('year','FONTSIZE',14)
    ylabel('ET [mm/h]','FONTSIZE',14)     
    box off
    
    figure(3)
    subplot(2,1,1)
    plot(time,Pgen)
    xlabel('time [year]','fontsize',14)
    ylabel('Precipitation(t) [mm/h]','fontsize',14)
    set(gca,'XAxisLocation','top','Ydir','reverse','fontsize',12,'tickdir','out');
    subplot(2,1,2)
    plot(time,Q)
    set(gca,'fontsize',12,'tickdir','out');
    box off
    ylabel('Q  [m^{3}/s]','fontsize',14)
    xlabel('time [years]','fontsize',14)
end


%*************************************************************************
%                     RESERVOIR ROUTING
%*************************************************************************

%reservoir parameters
Cq_sluice=0.6;            %discharge coefficient sluice gate [-]
Cq_spill=0.5;             %discharge coefficient spillway [-]
L=100;                    %spillway length [m]
p=18;                     %level of spillway with respect to the empty pool level [m]

%power plant parameters
D=3;                      %pipe diameter [m]
k=0.2*1e-3;               %pipe roughness [m]
eta=0.80;                 %turbines efficiency [-]
Lp=600;                   %pipe length [m]
deltah=50;                %difference of elevation between he empty pool level and the level of the tailrace [m]
minimum_level_for_HU=2;   %minimum level for hydropower production with respect to the empty pool level [m]
Qmax=70;                  %target discharge for the flood protection [m^3/s]

%determination of the minimum flow target(i.e. discharge that is exceeded 95% of the time)
%discharge duration curve 
sort_Q=sort(Q,'descend');
p_exceedence=(1:length(Q))/length(Q);
figure (4)
semilogy(p_exceedence,sort_Q)
title('Discharge Duration Curve','FONTSIZE',14)
xlabel('fraction of the time equalled or exceeded','FONTSIZE',14)                   
ylabel('Discharge  [m^{3}/s]','FONTSIZE',14) 
box off
%determination of the Q347 (i.e. discharge that is exceeded 95% of the time)
Q347=sort_Q(find(p_exceedence>0.95,1,'first'))
hold on
plot([0 0.95],[Q347 Q347],'--r')

Q347=sort_Q(round(length(Q)*0.95))

%************************  Volume Rating Curve *****************************

%area rating curve
load area_rating_curve.txt
level=area_rating_curve(:,1);          %"level": levels within the reservoir [m] above see level
level=level-level(1);                  %"level": levels with respect to the empty pool level (level(1))
lake_area=area_rating_curve(:,2);      %"lake_area": lake area corresponding to the "level" values

clear area_rating_curve                %clear the variable volume_rating_curve

%determination of the Volume rating curve (volume corresponding to the "level" values)
Volume=cumtrapz(level,lake_area); 


%plot area and volume rating curve
figure(1001)
subplot(1,2,1)
plot(level,lake_area,'ro')
xlabel('level  [m]','FONTSIZE',12)        %insert label x axis
ylabel('Area  [m^{2}]','FONTSIZE',12)     %insert label y axis
box off
subplot(1,2,2)
plot(level,Volume,'ro')
xlabel('level  [m]','FONTSIZE',12)
ylabel('Volume  [m^{3}]','FONTSIZE',12)
box off



%************************  Reservoir Routing ******************************

%variables preallocation
V=zeros(length(Q),1);               % volume time series
l=zeros(length(Q),1);               % level time series
sluice_A=zeros(length(Q),1);        % sluice gate opening time series
Qout=zeros(length(Q),1);            % outflowing dicharge time series
Q_HU=zeros(length(Q),1);            % discharge for hydropower time series
max_level_HU=zeros(9,1);            % maximum level where power plant is running
Vmax_HU=zeros(9,1);                 % volumes corresponds to maximum levels
Qg=zeros(length(Q),1);              % discharge through the gate
g=9.81;
QT=45;                              % turbine discharge
mean_annual_E=zeros(9,1);           % average energy produced in a year
p_flooding=zeros(9,1);              % probability that annual max discharge exceeds 70mÂ³/s

N_days=length(Q)/24;                % number of days of simulation


%initial condition
V(1)=interp1(level,Volume,p);       %initial volume (level equal to the spillway level)


for x=1:9                           %counter for different max_level_HU
    t=0;                            %counter of the number of hours in the simulation

    max_level_HU(x)=x+9;            %maximum level for hydroelectric use
    %the corresponding maximum volume for hydroelectric use:
    Vmax_HU(x)=interp1(level,Volume,max_level_HU(x));   
  
    for d=1:N_days                  %for loop on the number of days
   
        %decide if the power plant will work during this day
        level_midnight(d)=level_volume(Volume,level,V(t+1));
        if level_midnight(d)>=minimum_level_for_HU            flag=1;
        else
            flag=0;
        end

        
        for hT=1:24  %for loop on the number of hour in a day (once every hour)
            t=t+1;  %update counter
        
            %compute the level 
            l(t)=level_volume(Volume,level,V(t));
                    
            %compute Q_HU
            if hT>=12 && hT<=18 && flag==1
                Q_HU(t)=QT;  %[m^3/s]
            else
                Q_HU(t)=0;        
            end
        
            %compute sluice opening 
            Qcrit(t)=(V(t)+(Q(t)-Q_HU(t)).*3600-Vmax_HU(x));
            Qg(t)=max(Q347,min(Qcrit(t)./3600,Qmax));
            sluice_A(t)=Qg(t)/(Cq_sluice*sqrt(g*2*l(t)));
        
            %compute Qout
            if l(t)<=p
                Qout(t)=Cq_sluice*sluice_A(t).*sqrt(2*g*l(t));
            else
                Qout(t)=(Cq_sluice*sluice_A(t).*sqrt(2*g*l(t)))+(Cq_spill*L*sqrt(2*g*(l(t)-p).^3));
            end
        
            %update V (integration of the storage equation with the Euler scheme)
            V(t+1)=V(t)+((Q(t)-Qout(t)-Q_HU(t)).*3600);
        end
    end
V(end)=[];


%************************average annual energy production**********************

nu=1e-6;                %viscosity of water [L^2/]        
gamma_w=9806;           %specific weight of water[N/m^3]
K_in=0.5;               %assume entrance head losses factor is .5

v=QT/((D^2*pi)/4);      %velocity in pipe
Re=(v*D)/nu;            %Reynolds number

%calculation of friction factor
f=1;                    
for i=1:10
    f(i+1)=(1/(-2*log10(((k/D)/3.7)+(2.51/(Re*sqrt(f(i)))))))^2;
end
f=f(10);

%calculation of head for every hour
hf=(f*Lp*v^2)/(D*2*g);                     % friction head losses      
h0=K_in*(v^2/(2*g));                       % entrance head losses
hT=deltah+l-hf-h0;                                 

%calculation of power and energy
power=eta*gamma_w*Q_HU.*hT;                    % power for every hour, equal to energy produced in that hour [Wh]
mean_annual_E(x)=sum(power)/N_years_gen/1e9;   % average annual energy for given maximum HU level


%********************************flooding probality************************

%flooding probability is the ANNUAL probability where the max discharge exceeds 70m^3/s
Qout_reshape=reshape(Qout,8760,N_years_gen);
Floods=sort(Qout_reshape,'descend');
p_flooding(x)=sum(Floods(1,:)>71)./N_years_gen;   
end

% ***************************** Plot Results *******************************

if edit
    figure(1002)
    subplot(3,1,1)
    plot(time,Q,time,Qout,'-r');
    box off
    legend('Qin','Qout')
    ylabel('discharge [m^{3}/s]','fontsize',12)
    
    subplot(3,1,2)
    plot(time,V)
    box off
    ylabel('volume  [m^{3}/s]','fontsize',12)
    
    subplot(3,1,3)
    plot(time,l,time,p*ones(size(time)),time,minimum_level_for_HU*ones(size(time)))
    box off
    xlabel('t  [years]','fontsize',12)
    ylabel('level  [m]','fontsize',12)
    
    figure(1003)
    plot(max_level_HU,mean_annual_E,'o-r','markerfacecolor','red');
    xlabel('max level for HU [m]','FONTSIZE',12)    
    ylabel('Annual Energy [GWh]','FONTSIZE',12)
    box off
    
    figure(1004)
    plot(max_level_HU,p_flooding,'o-r','markerfacecolor','red');
    xlabel('max level for HU [m]','FONTSIZE',12)    
    ylabel('Flooding probability','FONTSIZE',12)
    box off
end
