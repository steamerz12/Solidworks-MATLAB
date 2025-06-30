%% MAE Midsizing Report CODE
clear
clc
%close all
%%%%NAME AND COMMENT EACH OF THESE TO debug
load('DeltaMdivCritical'); %  critical Mach number variations - digitized plots
load('MdivCO.mat'); %  Mach number division interpolation function
load('CritDeltaMdivV3.mat'); %  updated critical Mach number division data
load('tc_curve.mat'); %  thickness-to-chord ratio data for various airfoils
load('tcConv.mat'); %  thickness-to-chord ratios for conventional airfoils
load('tcCrti.mat'); %  thickness-to-chord ratios for supercritical airfoils
load('tcCrit.mat'); %  additional supercritical airfoil data
load('tcCritt.mat'); %  more supercritical airfoil data
load('ClMaxLand.mat'); %  maximum lift coefficient data for landing configurations
load('ClMTakeoff.mat'); %  maximum lift coefficient data for takeoff configurations
load('CLMaxClean.mat'); %  maximum clean configuration lift coefficients
load("ClClean0Sweep.mat"); %  clean configuration lift coefficient for 0-degree sweep wings
load('ClClean35Sweep.mat'); %  clean configuration lift coefficient for 35-degree sweep wings
load('ClCleanSweep15.mat'); %  clean configuration lift coefficient for 15-degree sweep wings
load('FuelFractionJT8D.mat'); %  fuel fraction data for JT8D engine configurations
load('Thrust.mat'); %  engine thrust capabilities data
load('MaxThrustSL.mat'); %  maximum sea level thrust 
load('MaxDryThrustSealevel'); %  maximum dry thrust data at sea level
load('SLTHRUST.mat'); %  sea level thrust 
load('KFactorNonWing.mat'); %  non-wing aerodynamic factors for drag calculation
load('Cfcoef.mat'); %  skin friction coefficient 
load('CRUISE.mat'); %  cruise condition parameters
load("CDPChangeLand.mat"); %  parasitic drag coefficient upon landing
load("CDChange.mat"); %  drag coefficient change data
load('TakeOffCdpdelta.mat'); %  change in drag coefficient during takeoff

%% SPECS/Variation
wingtype=0; % 0 for critical wingtype, 1 for conventional
AdvanceEngines=1; 
AluminummMat=0; 
CompositeMat=1; 

if CompositeMat==1
   WCM=.7; % Weight correction factor for wings 
   FCM=.85; % Fuselage correction factor 
   FECM=.9; % Fixed equipment correction factor 
   NacelePylonCM=.8; % Nacelle and pylon correction 
   c=1.2; % General correction coefficient 
else
   WCM=1;
   FCM=1;
   FECM=1;
   NacelePylonCM=1;
end
if AluminummMat==1
   wingAL=.94; % Correction factor for wing weight aluminumm
   FAL=.954; % Correction factor for fuselage weight 
   c=1.356; % 
else
   wingAL=1;
   FAL=1;
end
if AdvanceEngines==0
   ATSFC=1; % Standard thrust-specific fuel consumption
   ATENGINE=1; 
else
   ATSFC=.9; % Reduced fuel consumption for advanced engines
   ATENGINE=1.1; % Increased performance for advanced engines
end
M=.80; % Cruise Mach number
Range=3500; % Range in nautical miles
TOFL=6900; % Takeoff field length in feet
V=140; % Stall speed in knots
pax=210; % Number of passengers
NA=7; % Number of attendants
NAS=2; % Number of auxiliary staff
Wcargo=8000; % Weight of cargo in pounds
hcl=35000; % Initial cruise altitude in feet
q=1; % Indicates wing-mounted engines
xyz=0.45; % Fraction of fuel used before landing
Z=2; % Number of engines
Cl=.54; % Initial guess for lift coefficient
sweep=35; % Wing sweep angle in degrees
AR=8; % Aspect ratio
TR=.35; % Taper ratio
rho=.0007382; % Air density at cruise height
rhosl=.0023769; % Air density at sea level
sigma=rho/rhosl; % Density ratio SL
sigmaLand=.953; % Density ratio at landing
Talt=394.08; % Temperature at altitude
RQ=1717; % Universal gas constant for air
gam=1.4; % Heat capacity ratio
p=499.34; % Pressure at altitude
pr=p/2116.2; % Pressure ratio
Wfraction=0.32; %Typical Fuel Fraction for a commercial airliner of these specs
if wingtype==1
   jq=4; % Correction factor based on wing type
   ii=7; % Iteration parameter based on wing type
else
   jq=6;
   ii=6;
end
for j=1:4
   sweep=15+5*j;
   sweepcheck(1,j)=sweep;
end

for i=1:ii
   AR=jq+i;
   ARgraph(j,i)=AR;
   adjust=-.018;
   aj=.3;
    if sweep==35
       aj=.8;
    end 
end
while 1

while 1

while 1
    
if wingtype==1   
MdivDelta=MdivCO(Cl)*1.02;
elseif wingtype==0
MdivDelta=CritDeltaMdivV3(Cl);
end
Mdiv=(M+.004)-MdivDelta;
if wingtype==1 %%%%%finding tc loop
 tc=tovercConv(Mdiv, sweep);
    if sweep==20
     tc=(((sweep-20)/(25-20)*(-17.0631*Mdiv^3 + 43.9452*Mdiv^2 - 38.1861*Mdiv + ...
            11.2626-(-14.7915*Mdiv^3 + 37.4607*Mdiv^2 - 32.0814*Mdiv + 9.3546))+(-14.7915*Mdiv^3 +37.4607*Mdiv^2 - 32.0814*Mdiv + 9.3546)))*.93;
    elseif sweep==25
     tc=((sweep-25)/(30-25)*(-16.4307*Mdiv^3 + 44.1728*Mdiv^2 - 40.0666*Mdiv +12.3264- ...
            (-18.2575*Mdiv^3 + 47.0214*Mdiv^2 - 40.8591*Mdiv + 12.0510))+(-18.2575*Mdiv^3 + 47.0214*Mdiv^2 - 40.8591*Mdiv + 12.0510))*.9;
    end 

% Define sweep range adjustments
SWEEP_ADJUSTMENT_FACTOR = 0.95;
sweep_ranges = [20, 25; 25, 30; 30, 35]; % Define ranges in a matrix
correction_factors = [-12.8960, 32.6167, -27.9454, 8.1774; -15.4403, 39.6188, -34.3375, 10.1208];

% Calculate thickness to chord ratio based on sweep range
if wingtype == 0
    for idx = 1:size(sweep_ranges, 1)
        if sweep >= sweep_ranges(idx, 1) && sweep < sweep_ranges(idx, 2)
            range_factor = (sweep_ranges(idx, 2) - sweep) / (sweep_ranges(idx, 2) - sweep_ranges(idx, 1));
            tc = sum(correction_factors(idx, :) .* [Mdiv^3, Mdiv^2, Mdiv, 1]) * (1 - range_factor) * SWEEP_ADJUSTMENT_FACTOR;
            break;
        end
    end
end
weightLoadng=(cos(sweep*pi/180))^2*tc^2*AR; % Calc wing loading
CLmaxT=ClMTakeoff(weightLoadng); % Max takeoff lift coef.
CLmaxL=ClMLand(weightLoadng); % Max landing lift coef.
WSL=(V/1.3)^2*sigmaLand*CLmaxL/296; % Stall wing loading
Vcruise=M*sqrt(RQ*Talt*gam)*.5924; % Compute cruise speed
RangeAO=Range+200+.75*Vcruise; % Adjust nominal range

WSto=WSL/(1-xyz*Wfration); % Takeoff wing loading
WIC=.965*WSto; % Weight index coeff.
Clfinal=WIC/(1481*pr*M^2); % Final lift coef.

 if abs(Clfinal-Cl)<.001
     break
 else
     if Cl>Clfinal
         Cl=Cl-.0001;
     else
         Cl=Cl+.0001;
     end
 end
end
if Z==2
   TOFLVAR=0.0285*TOFL - 10.893;
elseif Z==3
   TOFLVAR=0.0316*TOFL - 8.8235;
elseif Z==4
   TOFLVAR=0.0327*TOFL - 0.8809;
end
WeightoverT70V=TOFLVAR/WSto*(sigmaLand*CLmaxT);
Vlo=1.2*(296*WSto/sigmaLand/CLmaxT)^.5;
Mlo=Vlo/661/sqrt(sigmaLand);
Mlo7Vl=.7*Mlo;
y=[45500 39120 34820 317501];
xq=[0 .15 .30 .45]; %AR correction
TMlo7Vl=interp1(xq,y,Mlo7Vl);
TSLST=interp1(xq,y,0);
WoverTy=WeightoverT70V*TMlo7Vl/TSLST-aj;
%% WEIGHT CALCULATIONS
% Define load factor for design safety and structural integrity considerations
n=1.5*2.5; % load factor

% Determine weight correction factors based on engine configuration
if q==1
    Kw=1.0; % Weight factor for wing-mounted engines
    Kts=.17; % Tail size factor for wing-mounted engines
end

%  fuselage weight based on passenger numbers and attendants
lfus=(3.76*pax/NA+33.2); % Length of fuselage adjusted by passenger count 
dfus=(1.75*NA+1.58*NAS+1); % Fuselage diameter adjustment based on crew and service staff

% thickness chord ratio
tca=tc+.03; % Adjusted thickness-to-chord ratio for airfoil

% Wing weight calculation
Ww=.00945*(AR)^.8*(1+TR)^.25/(tca^.4*cosd(sweep)*WSto^.695)*Kw*n^.5;

% Fuselage weight calculation incorporating fuselage length and diameter
Wfus=.6727*11.5*lfus^.6*dfus^.72*n^.3;

% Weight of landing gear assumed constant
WLG=.040; % Constant landing gear weight factor

% Nacelle and pylon weight calculated based on total weight over thrust ratio
Wnp=.0555/WoverTy*NacelePylonCM; % Weight of nacelles and pylons adjusted by thrust factor and material

% Tail weight determined by tail size factor and wing weight
Wtail=(Kts)*Ww; 
if roff==1
    Wtail=(Kts+.08/3)*Ww; % Additional tail weight adjustment for profile 1
end

% Weight of propulsion for using advanced specific fuel consumption factor
Wpp=1/(3.58*WoverTy)*ATSFC; 

% Fuel weight calculation w/ fuel fraction
Wfuel=1.0275*Wfration; 

% Payload weight sum of passenger and cargo weights
Wpayload=215*pax+Wcargo; 

% Fixed equipment weight considering passenger capacity, number of engines, and equipment factor
Wfe=(132*pax+300*Z*ATENGINE+260*2+170*ceil(pax/50))*FECM;
% Total weight equation to be solved using iterative root finding method
WEIGHTequation = @(wto) (Ww+Wtail)*WCM*wingAL * wto.^1.195 + (Wfus)*FCM*FAL * wto.^0.235+(WLG+Wnp+Wpp+Wfuel+.035*FECM-1)* wto + Wpayload+Wfe;
% Initial guess for the total weight optimization
initialGuess = 350000; 
% Solve the total weight equation using fzero to find optimal weight
wo = fzero(WEIGHTequation, initialGuess);
% wing area based on optimized weight
St=wo/WSto;
Bw=(AR*St)^.5;
mac=St/Bw;
% Thrust requried
Thrust=wo/WoverTy;
% Calculate thrust per engine
ThrustE=Thrust/Z;

%% DRAG
RnL=(.5*994.85)/.00034884;
maco=mac*1.3;
%%Wings CFI== 0.0826*(Rn/L*L)^-0.196
Mo=.5;
Zm=(2-Mo^2)*cosd(sweep)/sqrt(1-Mo^2*(cosd(sweep)^2));
Kwing=1+Zm*tc+100*tc^4;
CfiW= 0.0826*(RnL*mac)^-0.196;
SwetWing=2*(St-dfus*maco)*1.02; %%fix correction factor to 
Fwing=Kwing*CfiW*SwetWing;
%%Fuseluge
Swetfus=.9*pi*lfus*dfus;
cfifus=0.0826*(RnL*lfus)^-0.196;
FFfus=lfus/dfus;
Kfus=KFactorNonWing(FFfus);
Ffus=cfifus*Swetfus*Kfus;
%%Tail
Ftail=.38*Fwing;
%Nacelles
SwetNAC=2.1*(ThrustE)^.5*Z;
cfinac= 0.0826*(RnL*mac)^-0.196;
Knac=1.25;
Fnac=Knac*cfinac*SwetNAC;
%Pylons
Fpylon=.20*Fnac;
Ftotalsum=(Fwing+Ffus+Ftail+Fnac+Fpylon)*1.06;
Cdop=Ftotalsum/St;
eez=1/(1.035+.38*Cdop*pi*AR);
%% CLIMB
% air density ratio for climb phase
sigma=.53317;

% average weight 
Wcl=(1+.965)/2*wo;

% Compute the climb speed (VcCl) using the drag (Ftotal) and Oswald efficiency number (e)
VcCl=1.3*12.9/(Ftotalsum*eez)^.25*(Wcl/(sigma*Bw))^0.5;

% Mach number during climb
Mcll=VcCl/614.3464;  % Converts climb speed to Mach number assuming standard sea level speed of sound= knots

% Required thrust for climb (Trcl) based onspeed, weight, and drag 
Trcl=sigma*Ftotalsum*VcCl^2/296+94.1*(Wcl/Bw)^2/(sigma*eez*VcCl^2);
TaG=(-9993.8*Mcll^3 + 23895*Mcll^2 - 25308*Mcll + 27030 -16043*Mcll^4 + 35397*Mcll^3 - 24488*Mcll^2 + 3521.4*Mcll + 15536)/2; F% Combined equation for available thrust

% Calculate specific fuel consumption (SFC) for the climb phase, averaged for different Mach numbers
SFC=(0.3664*Mcll+0.344+0.4238*Mcll+0.3235)/2; % Average SFC acrossthe bored

% Calculate available thrust (Ta) adjusted by actual thrust and available gradient
Ta=Thrust/45500*TaG;  % Normalize thrust factor


RC=101*(Z*Ta-Trcl)*VcCl/Wcl;  % RATE OF CLIMB Factor 101 converts to appropriate units, Z is number of engines

%time required to reach cruise altitude (hcl)
TimeCL=hcl/RC;  % hcl is cruise altitude in feet

% Convert climb time to range covered during climb
RangeCL=VcCl*TimeCL/60;  % Convert time to hours and multiply by climb speed

% fuel used during climb
WfCL=Z*Ta*SFC*TimeCL/60;  % Calculate total fuel consumption based on time, thrust, and SFC

%% Total Range
wo=wo-WfCL;
w1=(1-Wfration)*wo;
CLavg=(wo+w1)/(2*St)/(1481*.2360*M^2);
Cdi=CLavg^2/(pi*AR*eez);
CDR=Cdi+Cdop+.0010;
LoD=CLavg/CDR;
TrR=(wo+w1)/2/LoD;
TrJ9=TrR*45500/ThrustE/Z;
SFCR=CatCruiseH(TrJ9)*ATSFC;
Rcruise=Vcruise/SFCR*LoD*log(wo/w1);
R=RangeCL+Rcruise;
if abs(R-RangeAO)>30 %%ADJUSTMENTS TO MAKE GRAPHS WORK< 
   if R>RangeAO
       adjust=adjust-.0005;
   elseif RangeAO>R
       adjust=adjust+.0005;
   end
else
end
end
% Thrust Check
CLIC=wo/St/1481/.2360/M^2;
CdiC=CLIC^2/pi/8/eez;
CDC=CdiC+Cdop+.0010;
LoDC=CLIC/CDC;
TreqC=Wcl/LoDC/Z;
TreqCJ9=TreqC*45500/ThrustE;
% 1 Climb Gradient
CL1=CLmaxT/(1.2)^2;
CLtoOCLmaxT1=1/1.2^2;
deltaCDo1=TakeOffCdpdelta(CLtoOCLmaxT1);
CD1=2*Cdop+deltaCDo1+CL1^2/(pi*AR*eez);
LoD1=CL1/CD1;
Treq1=wo/LoD1;
Ta1eng=ThrustE/45500*ThrustSL(Mlo);
Grad1=(((Z-1)*Ta1eng-Treq1)/wo)*100;
a=0;
%Second Gradient
CD2=CD1-Cdop;
LoD2=CL1/CD2;
Treq2=wo/LoD2;
Grad2=((Z-1)*Ta1eng-Treq2)/wo*100;
b=2.4;
%3rd Gradient
% if sweep>=0 &&sweep<15
%     CLClean=ClClean0Sweep(tc)*((15-sweep)/(15-0))+ClCleanSweep15(tc)*(1-(15-sweep)/(15-0));
% elseif sweep>=15 &&sweep<35
%     CLClean=ClCleanSweep15(tc)*((35-sweep)/(35-15))+ClClean35Sweep(tc)*(1-(35-sweep)/(35-15));
% elseif sweep==35
% CLClean=ClClean35Sweep(tc);
% end
CLClean=CLMaxClean(tc,sweep);
V3=1.2*(296*WSto/.925/CLClean)^.5;
M3=V3/659;
CL3=CLClean/1.2^2;
CD3=Cdop+CL3^2/pi/AR/eez;
LoD3=CL3/CD3;
Treq3=wo/LoD3;
Taeng3=ThrustE/45500*MaxThrustSL(M3);
Grad3=((Z-1)*Taeng3-Treq3)/wo*100;
% Approach Gradient
CL4=CLmaxT/1.3^2;
ClAp_Clmax=1/1.3^2;
deltaCDp4=TakeOffCdpdelta(ClAp_Clmax);
CD4=deltaCDp4+Cdop+CL4^2/pi/AR/eez;
LoD4=CL4/CD4;
Treq4=Wlanding4/LoD4;
V4=(296*WSL/(.953*CL4))^.5;
M4=V4*1.6878/sqrt(gam*1717*543.67);
Ta4=ThrustE/45500*MaxThrustSL(M4);
Grad4=((Z-1)*Ta4-Treq4)/Wlanding4*100;
d=2.1;
%Landing Grad
CL5=CLmaxL/1.3^2;
CL_CLmax5=1/1.3^2;
deltaCDp5=CDPChangeLandingN(CL_CLmax5);
CD5=deltaCDp5+2*Cdop+CL5^2/(pi*AR*eez);
LoD5=CL5/CD5;
Treq5=Wlanding4/LoD5;
V5=140;
M5=V5*1.6878/sqrt(gam*1717*543.67);
Ta5=ThrustE/45500*ThrustSL(M5);
Grad5=((Z)*Ta5-Treq5)/Wlanding4*100;
eez=3.2;
c=1.2;
%% Direct Operating COSTS
   D=Range*1.15;
   Tgm=.25;
   Tcl=.18;
   Td=0;
   Tam=.10;
   Tcr=(D*1.02+20-RangeCL*1.15)/Vcruise/1.15;
   TB=(Tgm+Tcl+Td+Tcr+Tam);
   VB=D/TB;
   Fcl=WfCL;
   Fcr_Fam=TrR*SFCR*(Tcr+Tam);
   FB=Fcl+Fcr_Fam;
%Flight Crew
   P=Wpayload/2000;
   dollarbhour=17.849*(Vcruise*1.15*wo/10^5)^.3+40.83;
   Ctma=dollarbhour/(VB*P);
%Fuel and Oil
   CF=.28*(1/6.4);
   Cott=2.15;
   Ctmb=(1.02*FB*CF+Z*Cott*TB*.135)/(D*P);
 %Hull Insurance
   Wa=wo*(1-Wfration)-Wpayload-Wpp*wo;
   Ca=2.4E6+87.5*Wa;
   Ce=590000+16*ThrustE;
   CT=Ca+Z*Ce;
   IRa=.01;
   U=630+4000/(1+1/(TB+.5));
   Ctmc=IRa*CT/(U*VB*P);
  %Maintence
   Kfha=4.9169*log10(Wa/10^3)-6.425;
   Kfca=.21257*(log10(Wa/10^3))^3.7375;
   Tf=TB-Tgm;
   Rl=8.60;
   Ctmd=(Kfha*Tf+Kfca)*Rl/(VB*TB*P);
  %Airframe Mat
   Cfha=1.5994*Ca/10^6+3.4263;
   Cfca=1.9229*Ca/10^6+2.2504;
   Ctme=(Cfha*Tf+Cfca)/(VB*TB*P);
  %Labor
   Kfhe=Z*(ThrustE/10^3)/(.82715*(ThrustE/10^3)+13.639);
   Kfce=.20*Z;
   Ctmf=(Kfhe*Tf+Kfce)*Rl/(VB*TB*P);
  %Engine mats
   Cfhe=(28.2352*Ce/10^6-6.5716)*Z;
   Cfce=(3.6698*Ce/10^6+1.3685)*Z;
   Ctmg=(Cfhe*Tf+Cfce)/(VB*TB*P)*ATENGINE;
  % Maintence sum
   CtmT=(Ctmg+Ctmf+Ctme+Ctmd)*2;
   %%Depreciation
   CtmD=(CT+.06*(CT-Z*Ce)+.3*Z*Ce)/(14*U*VB*P);
   %%%Total DOC
   CtmTot=(CtmT+Ctma+Ctmb+Ctmc+CtmD)*P/pax;
   DOC1(j,i)=CtmTot;
   MaxW(j,i)=wo;
if wingtype==1
   if sweep==25 && AR==8
       finalWeight=wo;
       LCFinal=LoD;
       FinalVelocity=Vcruise;
       Cruise_SFC=SFCR;
       Final_FuelFraction=Wfration;
       Final_PayloadWeight=61750;
       OEW=wo-Wpayload;
   end
elseif wingtype==0
       if sweep==25 && AR==9
       finalWeight=wo;
       LCFinal=LoD;
       FinalVelocity=Vcruise;
       Cruise_SFC=SFCR;
       Final_FuelFraction=Wfration;
       Final_PayloadWeight=61750;
       OEW=wo-Wpayload;
       end
end
end
end
%% Graphs
% Constants for display
SWEEP_INCREMENT = 5;
BASE_SWEEP = 15;
NUMBER_OF_SWEEPS = 4;

% Plotting Direct Operating Costs
figure;
hold on;
for sweepIndex = 1:NUMBER_OF_SWEEPS
    currentSweep = BASE_SWEEP + SWEEP_INCREMENT * sweepIndex;
    sweepLabel = ['Sweep Angle = ', num2str(currentSweep), ' degrees'];
    plot(ARgraph(sweepIndex, :), DOC1(sweepIndex, :), 'DisplayName', sweepLabel);
end
xlabel('Aspect Ratio (AR)');
ylabel('Direct Operating Cost ($/ton mile)');
airplaneType = 'Conventional Airplane Design';
if wingtype == 0 && AdvanceEngines == 1 && CompositeMat == 1
    airplaneType = 'Advanced Technology Airplane design';
end
title(airplaneType);
legend show;
hold off;

% Plotting Take-Off Weight
figure;
hold on;
for sweepIndex = 1:NUMBER_OF_SWEEPS
    currentSweep = BASE_SWEEP + SWEEP_INCREMENT * sweepIndex;
    sweepLabel = ['Sweep Angle = ', num2str(currentSweep), ' degrees'];
    plot(ARgraph(sweepIndex, :), MaxW(sweepIndex, :), 'DisplayName', sweepLabel);
end
xlabel('Aspect Ratio (AR)');
ylabel('Take-Off Weight (lbs)');
title(airplaneType);
legend show;
hold off;

% Displaying configuration and results
configurationInfo = ['Using correction factor x=', num2str(xyz), ...
                     '. Take Off Weight= ', num2str(wo), ' lbs.'];
disp(configurationInfo);

engineType = 'Regular Engines Used';
if AdvanceEngines == 1
    engineType = 'Advanced Engines Used';
end
disp(engineType);
