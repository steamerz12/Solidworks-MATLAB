 % clc;clear all;
load('AR_CP_E.mat')  %Aspect Ratio vs E accounting for CP   (inputs: AR,CDP) outputs E
load('K_LD.mat')     %L/D vs K           (input:Fineness ratio LD) output K
load('RNvsCf.mat')   %Reynold # vs Skin Frict. coef. Cf      (input true Reynold number) gives K 
% Z=((2-0.5^2)*cos(sweepbangle))/(sqrt(1-0.5^2 *cos(sweepbangle)^2));                 
% K=(1+Z*tc+100*tc^4);
% use this formula ^^^ to get K when only t/c and sweep back angle is given

V0_values = [];  % Initialize an array to store V0 values
LDratio_values = [];  % Initialize an array to store L/D ratio values
Di_values=[] 
Dp_values=[]
Dtotal=[]
%% wing
%define unique constants given by 
hold on
for V0=250:15:800
rho=0.0008754;
MU=3.025*10^-7;
b=93;
Cr=17;
taper=0.19;
sweepbangle=24;
tc=0.108;

q0=0.5*(rho)*V0^2;  %Dynamic 
RE_l=rho*V0/MU ;  %reynold per unit length

Ct=taper*Cr;
mac=(2/3)*(Cr+Ct-Cr*Ct/(Cr+Ct));
RE=RE_l*mac;
Cf=RNvsCf(RE)*1.0798 *10^-3 ; %Correction factor for Reynolds number (very slighly off still)
 Z=((2-0.5^2)*cos(sweepbangle))/(sqrt(1-0.5^2 *cos(sweepbangle)^2));
 K=(1+Z*tc+100*tc^4);

 Swet=1000*(0.83)*2*1.02;   %recacluate this once solidworks model complete  Sexp*2*1.03
 F_w= K*Cf*Swet; 
CDP_w=F_w/1000;
%% Fuselage
Lf=105;
Df=12;

Swet=0.8*pi*Lf*Df;  

LD=Lf/Df;
kf=K_LD(LD);
RE_f=RE_l*Lf;
Cf_f=RNvsCf(RE_f)*1.0598 *10^-3 ;

f_f=kf*Cf_f*Swet ;
CDP_f=f_f/1000;

%% Horizonal tail 
tc=0.09;
sweepbangle=31.6;
tr=0.35;
Cr=11.1;
Z=((2-0.5^2)*cos(sweepbangle))/(sqrt(1-0.5^2 *cos(sweepbangle)^2));
K=(1+Z*tc+100*tc^4);
Swet=261*2*1.02;
mac=(2/3)*Cr*(1+tr-tr/(1+tr));
RE_Ht=RE_l*mac;

Cf_Ht=RNvsCf(RE_Ht)*1.0598 *10^-3 ;
f_Ht=K*Cf_Ht*Swet;
CDP_Ht=f_Ht/1000;

%% Vertical Tail

Swet=161*2*1.02;
tc=0.09;
sweepbangle=43.5;
tr=0.8;
mac=(2/3)*15.5*(1+0.8-0.8/(1+0.8));
RE_Vt=RE_l*mac;
Cf_Vt=RNvsCf(RE_Vt)*1.0598 *10^-3 ;
Z=((2-0.5^2)*cos(sweepbangle))/(sqrt(1-0.5^2 *cos(sweepbangle)^2));
K=(1+Z*tc+100*tc^4);
f_Vt=K*Cf_Vt*Swet;
CDP_Vt=f_Vt/1000;

%% Pylon
Swet=117;
tc=0.06;
sweepbangle=0;
Z=((2-0.5^2)*cos(sweepbangle))/(sqrt(1-0.5^2 *cos(sweepbangle)^2));
K=(1+Z*tc+100*tc^4);
tr=1;
mac=16.2; %Taper ratio of 1 means mac does not change
RE_py=RE_l*mac;

Cf_Py=RNvsCf(RE_py)*1.0598 *10^-3 ;
f_py=K*Cf_Py*Swet;
CDP_py=f_py/1000;

%% Nacelles

Swet=455;
LD=5;
l=168;

kf=K_LD(LD);
RE_N=RE_l*l;

Cf_N=RNvsCf(RE_N)*1.462 *10^-3 ;%correction factor changes when reynold number reaches past RE>10^8  (plot digitizer not 100% accurate)

f_N=kf*Cf_N*Swet;
CDP_N=f_N/1000;

%% Flap Hinges
fh=0.15;
CDP_fh=0.15/1000;

Ftotal=(f_N+f_Vt+f_Ht+f_py+f_f+F_w+fh)*1.125;%multiply by 1.125 to get answer similar to HW 5 

% CDP=Ftotal/1000;    %METHOD 1
CDP=CDP_N+CDP_py+CDP_Vt+CDP_Ht+CDP_f+CDP_w+CDP_fh ;  %Method 2 (slight different return for e)

AR=b^2/1000;

e=AR_CP_E(AR,CDP);


Cl=98000/(q0*1000);
Cdi=Cl^2 /(pi *AR*e);

Di=q0*Cdi*1000;
Dp=q0*CDP*1000;
D=Di+Dp;
CD=Dp+Di;
LDratio=98000/D
V0_values = [V0_values, V0];  %stores calulated values in each iteration of loop 
    LDratio_values = [LDratio_values, LDratio];
    Di_values=[Di_values,Di];
    Dp_values=[Dp_values,Dp];
    Dtotal=[Dtotal,D];
end
plot(V0_values,LDratio_values)
title('Flight speed vs Lift to Drag ratio')
xlabel('Airspeed ft/s')
ylabel('Lift/Drag')
hold off
figure
plot(V0_values,Dtotal)
hold on
plot(V0_values,Dp_values)
hold on
plot(V0_values,Di_values)
hold off
legend('Total Drag','Profile Drag','Induced Darg')
xlabel('Airspeed ft/s')
ylabel('D (lb)')