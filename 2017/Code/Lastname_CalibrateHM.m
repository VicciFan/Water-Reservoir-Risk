clear all               %clear all the variables allocated in the workspace
close all               %close all the open figures
clc                     % clear command window
doTest = 1;             % if doTest=1, check on mass balance for the hydrological model is performed;
                        % otherwise, doTest=0 

% load input files (see detailed description of the headers of the txt files)
load P.txt              %precipitation data
load Q_obs.txt          %discharge data
load temperature.txt    %mean montly temperature
load kc.txt             %mean monthly crop coefficient

%*************************************************************************
% PARAMETERS
%*************************************************************************

sw=...          %wilting point []
s1=...          %soil moisture above which plants transpire at kc*ETO []
n=...           %porosity [-]
BaseFlow=...    %base flow [m^3/s]
tsup=...        %superficial residence time [h]
Area=...        %area of the basin [km^2]
lat=...         %latitude of the basin [degrees]

N_years=length(P)/24/365;  %number of years of data

day_month=[31 28 31 30 31 30 31 31 30 31 30 31]; %"day_month": number of days for each month
month_end=cumsum(day_month);                     %"month_end": last day of each month
month_start=month_end-day_month+1;               %"month_start": first day of each month

%plot precipitation and discharge time series
time=2000+(0:length(P)-1)/(24*365);  %time points in years

%*************************************************************************
% COMPUTE POTENTIAL EVAPOTRANSPIRATION
%*************************************************************************

% ...to be completed, use the Thornthwaite method, as in exercise 1


%*************************************************************************
% TEST HYDROLOGICAL MODEL
%*************************************************************************

% give arbitrary values to parameters (Ksat, c, tsub, z) in order to test
% the function Lastname_HydroModel
Ksat=...; % to be completed

Q= Lastname_HydroModel(..);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%

%*************************************************************************
% SIMULATED ANNEALING
%*************************************************************************
N_iter=...; % number of iterations for the simulated annealing algorithm
cooling_rate=...; % decay of probability of jumping

% prealocate vectors
Ksat_MC = zeros(N_iter,1); c_MC = zeros(N_iter,1);
tsub_MC = zeros(N_iter,1); z_MC = zeros(N_iter,1);
NS = zeros(N_iter,1);

% initial values for parameters that require calibration
% Choose meaningful values
Ksat=...; 
c=...; 
tsub=...; 
z=...;
Ksat_MC(1)=Ksat; c_MC(1)=c; 
tsub_MC(1)=tsub;  z_MC(1)=z; 
% standard deviations for parameters that require calibration 
sigmaKsat=...;
sigmac=...;
sigmatsub=...;
sigmaz=...;

Q = ...; %use Lastname_HydroModel.m

NS_fun = @(Q) ...;

NS(1)= NS_fun(Q);
NS_current=NS(1);
% simulated annealing temperature
Temp_SA=@(x) ...;

indAccept=0; OldAccept=0;
tic
for indAll = 1:N_iter % iteration count
    % use truncated normal distribution (function TruncNormRnd) to obtain new parameter set
    % use realistic under and upper limits for the parameters
    Ksat_new = ...;  
    c_new    = ...;
    tsub_new = ...;
    z_new = ...;
    
    % run model and evaluate fit
    Q = ...;
    NS_new = ...;
    % accept or reject new parameter set
    if  % accept if new NS is higher, otherwise accept with probability based on simulated annealing
            indAccept=indAccept+1; % update index counting accepted values in prealocated vector
            
            % complete
            
    end
    
    % print information on MC
    if mod(indAll,1000)==0;
        fprintf('\n');
        fprintf('%0.2f %% COMPLETE\n',indAll/N_iter*100);
        fprintf('Iterations: %d    Chain length: %d   Acceptance rate %0.2f%%\n',indAll,indAccept,indAccept/indAll*100);
        fprintf('Total elapsed time = %.1f s',toc);
        fprintf(' Local accept. rate = %.2f%%',(indAccept-OldAccept)/1000*100);
        fprintf('Last point:\n');
        fprintf('Ksat = %0.2f mm/h   c = %0.2f    tsub = %0.1f h    z = %0.0f m\n',Ksat,c,tsub,z);
        fprintf('NS = %0.4f    Temperature = %0.4f',NS(indAccept),Temp_SA(indAll));     
        OldAccept=indAccept;
        disp(' ')
    end
end
       
% cut end
Ksat_MC(indAccept+1:end)=[]; c_MC(indAccept+1:end)=[]; tsub_MC(indAccept+1:end)=[]; z_MC(indAccept+1:end)=[]; 
NS(indAccept+1:end)=[];    

%*************************************************************************
% PLOT RESULTS
%*************************************************************************
% compare time series of discharge
figure(101);clf
set(gcf,'units','centimeters','position',[0 0 20 11]);
hold on; l1=plot(time,Q); l2=plot(time,Q_obs);
legend([l1 l2],'modelled','observed','location','northeast'); legend boxoff
ylabel('Discharge [m^3/s]')

% plot Markov Chains
figure(102);
set(gcf,'units','centimeters','position',[0 0 20 15]);
subplot(2,3,1); plot(Ksat_MC)
set(gca,'xlim',[1 indAccept])
ylabel('K_{sat} [mm/h]'); 
subplot(2,3,2); plot(c_MC)
set(gca,'xlim',[1 indAccept])
ylabel('c [-]');
subplot(2,3,3); plot(tsub_MC)
set(gca,'xlim',[1 indAccept])
ylabel('t_{sub} [h]');
subplot(2,3,4); plot(z_MC)
set(gca,'xlim',[1 indAccept])
ylabel('z [m]');
subplot(2,3,5); plot(NS)
set(gca,'xlim',[1 indAccept])
ylabel('NS [-]');
subplot(2,3,6); 
% plot values of Temp_SA corresponding to accepted parameter sets
ylabel('NS [-]');

% Plot time series of P, R, I, s, L, ET for the calibrated model
figure(103)
set(gcf,'units','centimeters','position',[0 0 20 15]);
subplot(2,3,1)
plot(time,P)
xlabel('year','FONTSIZE',14)
ylabel('precipitation [mm/h]','FONTSIZE',14)
xlim([2000 2006]);
set(gca,'XTick',2000:3:2006);
box off

subplot(2,3,2)
plot(time,R)
xlabel('year','FONTSIZE',14)
ylabel('runoff [mm/h]','FONTSIZE',14)    
xlim([2000 2006]);   
set(gca,'XTick',2000:3:2006);
box off

subplot(2,3,3)
plot(time,I)
xlabel('year','FONTSIZE',14)
ylabel('infiltration [mm/h]','FONTSIZE',14)  
xlim([2000 2006]);
set(gca,'XTick',2000:3:2006);
box off

subplot(2,3,4)
plot(time,s(1:end-1))
xlabel('year','FONTSIZE',14)
ylabel('soil moisture','FONTSIZE',14)     
xlim([2000 2006]);  
set(gca,'XTick',2000:3:2006);
box off

subplot(2,3,5)
plot(time,L)
xlabel('year','FONTSIZE',14)
ylabel('Leakage [mm/h]','FONTSIZE',14)      
xlim([2000 2006]); 
set(gca,'XTick',2000:3:2006);
box off

subplot(2,3,6)
plot(time,ET)
xlabel('year','FONTSIZE',14)
ylabel('ET [mm/h]','FONTSIZE',14)     
xlim([2000 2006]);
set(gca,'XTick',2000:3:2006);
box off