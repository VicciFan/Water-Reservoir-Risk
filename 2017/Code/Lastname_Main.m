clear all               %clear all the variables allocated in the workspace
close all               %close all the open figures
clc                     %clear command window

edit=1;                 %if "edit"=0 no figures
                        %if "edit"=1 show figures
                        
%load input files
load P.txt              %precipitation data

%*************************************************************************
% PARAMETERS
%*************************************************************************

N_years_gen=100;           %number of years of generated rainfall
N_years=length(P)/365/24;  %number of years of precipitation data

%*************************************************************************
% RAINFALL GENERATION
%*************************************************************************

%upscaling from hourly rainfall to daily rainfall
dailyP= %...TO BE COMPLETED


%************************ Monthly Statistics ******************************

day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; %"day_month": number of days for each month
month_end=cumsum(day_month);                     %"month_end": last day of each month
month_start=month_end-day_month+1;               %"month_start": first day of each month

% preallocation of variables
monthly_mean_P=zeros(12,1);   %monthly mean precipitation (counting also non rainy days)
monthly_std_P=zeros(12,1);    %monthly precipitation standar deviation (counting also non rainy days)
monthly_lambda=zeros(12,1);   %monthly events arrival rate
monthly_alpha=zeros(12,1);    %monthly mean precipitation (counting only the rainy days)


%... TO BE COMPLETED: evaluate monthly_mean_P, monthly_std_P, monthly_lambda, monthly_alpha
% for the measured time series of precipitation

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
deltat = 1;

%generate the value of rainfall
for y=1:N_years_gen            %for loop on the years
    d=0;                       %counter for the day of the year, inizialized every years
    for m=1:12                 %for loop on the months
        for i=1:day_month(m)   %for loop on the days of the month
            d=d+1;             %update days counter
            
            %here you need to generate the value of rainfall for this
            %day: dailyPgen_day_year(d,y) with model parameters of the
            %corresponding month (monthly_lambda(m) and monthly_alpha(m))
            
            %first you need to decide whether or not there is  a rainfall event during this
            %day: an event occurs with probability monthly_lambda(m)*deltat.
            %deltat is equal to 1 day
            
            %if an event occurs, you need to assign a precipitation [mm/day] extracted (with the inverse
            %transformation method) from an exponential distribution with mean monthly_alpha(m)
        end
    end
end

%compute statistics of generated rainfall
monthly_mean_P_gen=zeros(12,1);
monthly_std_P_gen=zeros(12,1);
monthly_lambda_gen=zeros(12,1);
monthly_alpha_gen=zeros(12,1);

%... TO BE COMPLETED: evaluate monthly_mean_P_gen, monthly_std_P_gen, 
% monthly_lambda_gen, monthly_alpha_gen for the generated time series of precipitation

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
% HYDROLOGIC MODEL
%*************************************************************************

%downscaling of daily rainfall to hourly rainfall 
dailyPgen=dailyPgen_day_year(:);  %transform in a column vector
Pgen=downscaling(dailyPgen);      %generated rainfall in [mm/h]

% run the calibrated hydrological model with "Pgen" as forcing rainfall
%...TO BE COMPLETED

%**************************************************************************
%**************************** END OF PART 1 *******************************
%**************************************************************************


%*************************************************************************
% RESERVOIR ROUTING
%*************************************************************************

%reservoir parameters
Cq_sluice=...            %discharge coefficient sluice gate [-]
Cq_spill=...             %discharge coefficient spillway [-]
L=...                    %spillway length [m]
p=...                     %level of spillway with respect to the empty pool level [m]

%power plant parameters
D=...                       %pipe diameter [m]
k=...                       %pipe roughness [m]
eta=...                     %turbines efficiency [-]
Lp=...                      %pipe length [m]
deltah=...                  %difference of elevation between the empty pool level and the level of the tailrace [m]
minimum_level_for_HU=...    %minimum level for hydropower production with respect to the empty pool level [m]
QT = ...                    % m^3/s turbine discharge
Qlim = ...                  % target discharge for the flood protection [m^3/s]

% determination of the minimum flow target(i.e. discharge that is exceeded 95% of the time)
%... TO BE COMPLETED

% ************************  Volume Rating Curve *****************************

%area rating curve
load area_rating_curve.txt
level=area_rating_curve(:,1);          %"level": levels within the reservoir [m] above see level
level=level-level(1);                  %"level": levels with respect to the empty pool level (level(1))
lake_area=area_rating_curve(:,2);      %"lake_area": lake area corresponding to the "level" values

clear area_rating_curve                %clear the variable volume_rating_curve


%determination of the Volume rating curve (volume corresponding to the "level" values)
Volume= %... TO BE COMPLETED

% ****************************** plot *************************************
if edit
    figure(1001)
    subplot(1,2,1)
    plot(level,lake_area,'ro')
    xlabel('level  [m]','FONTSIZE',14)
    ylabel('Area  [m^{2}]','FONTSIZE',14)
    box off
    subplot(1,2,2)
    plot(level,Volume,'ro')
    xlabel('level  [m]','FONTSIZE',14)
    ylabel('Volume  [m^{3}]','FONTSIZE',14)
    box off   
end

%************************  Reservoir Routing ******************************

%variables preallocation
V=zeros(length(Q_gen),1);           % volume time series
l=zeros(length(Q_gen),1);           % level time series
sluice_A=zeros(length(Q_gen),1);    % sluice gate opening time series
Qout=zeros(length(Q_gen),1);        % outflowing dicharge time series
Q_HU=zeros(length(Q_gen),1);        % discharge for hydropower time series

N_days=length(Q_gen)/24;            % number of days of simulation

%initial condition
V(1)=interp1(level,Volume,p);   %initial volume (level equal to the spillway level)

t=0; %counter of the number of hours in the simulation
for d=1:N_days  %for loop on the number of days
    %decide if the power plant will work during this day
    
    for h=1:24  %for loop on the number of hour in a day
        t=t+1;  %update counter
        
        %compute the level (the function level_volume is faster than the interp1 function)
        %... TO BE COMPLETED
        
        %compute Q_HU
        %... TO BE COMPLETED
        
        %compute sluice opening
        %(for this week, we decide the opening so that the flow through the
        %gate equals the minumum flow target. Next week we will see how to
        %control the gate opening to implement the flood control practice)
        %... TO BE COMPLETED
        
        %compute Qout
        %... TO BE COMPLETED
        
        %update V (integration of the storage equation with the Euler scheme)
        %... TO BE COMPLETED
    end
end
V(end)=[];

% plot for a maximum level of hydroelectric use of 15m
if %max_level_hu is equal to 15 [m]
	if edit
	    figure(1002)
	    subplot(3,1,1)
	    plot(time,Q_gen,time,Qout,'-r');
	    box off
	    legend('Qin','Qout')
	    ylabel('discharge [m^{3}/s]','fontsize',14)
	    
	    subplot(3,1,2)
	    plot(time,V)
	    box off
	    ylabel('volume  [m^{3}/s]','fontsize',14)
	    
	    subplot(3,1,3)
	    plot(time,l,time,p*ones(size(time)),time,minimum_level_for_HU*ones(size(time)))
	    box off
	    xlabel('t  [years]','fontsize',14)
	    ylabel('level  [m]','fontsize',14)
	end
end

% compute mean annual energy production and the probability that the hourly
% discharge Qout is larger than Qmax = 150 m3/s 
% (for numerical reasons, the condition should read: Qout > 151 m3/s)
%... TO BE COMPLETED

% *************************  Plot Results  ********************************
if edit
    figure(1003)
    subplot(1,2,1)
    plot(max_lev_hu,E,'r-o')
    xlabel('Max level for HU [m]','fontsize',14)
    ylabel('Annual energy [GWh]','fontsize',14)
    subplot(1,2,2)
    plot(max_lev_hu,p_flood,'r-o')
    xlabel('Max level for HU [m]','fontsize',14)
    ylabel('Flooding probability [%]','fontsize',14) 
end
