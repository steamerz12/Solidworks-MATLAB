%% HW 5 p2
Pa=101300 ;       % ambient kpa  (Mach 0)
Ta=288.2 ; 
P1= 18750 ;      % Pressure at altitude 1 (Mach 0.85) 
T1= 216.75 ; 
P2= 7170 ;       %Mach 2
T2=T1;
P3= 2097 ;
T3=T2 ; 

%Adiabetic Efficiencies and Gamma 
Nd=0.97 ;
Gd=1.4 ;

Nc=0.85;
Gc=1.37;

Nb=1;
Gb=1.35;
Nt=0.9;
Gt=1.33;
Nn=0.98;
Gn=1.36; 

Nf=0.85;
Gf=Gd;
Nfn=0.97;
B=5;

%% M=0 
M0=0;
U= sqrt(1.4*287*Ta) * M0;
figure
% Compressor inlet 
T02=Ta*(1+(1.4-1)/2*M0^2) ;
P02=Pa*(1+Nd*(T02/Ta-1))^(Gd/(Gd-1));


Pcr= linspace(2,100); % Compressor pressure ratio (X-axis value) 
P03=P02.*Pcr ;
T03=T02.*(1+(1/Nc).*(Pcr.^((Gc-1)/Gc) -1) );
T04=1500 ;% Tmax according to textbook
P04=P03;
% Turbo Fan START
Q=45000000; % KJ/KG
Cp1=1107; %KJ/KG 
f=(T04./T03-1)./(Q./(Cp1.*T03)-T04./T03);

Prf=1.5;
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1)).*(287.*T08).*(1-(Pa./P08).^((Gf-1)/Gf)) ).^(1/2);

T05=T04-(T03-T02)-B*(T08-T02);
P05=P04.*(1-(1/Nt).*(1-T05./T04)).^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1)).*(287.*T06).*(1-(Pa./P06).^((Gn-1)/Gn)) ).^(1/2);

ST=(1+f).*Ue + B*Uef-(1+B)*U;
STreal=(ST(3:42));
TSFC= 1000.*f./ST;
TSFC=TSFC(3:42);
PCR=(Pcr(3:42));
semilogx(PCR,STreal/1000,'b')
xlim([2,100])
ylim([0,2])
yyaxis right
semilogx(PCR,TSFC,'b')
ylim([0 0.04])
hold on

T04=1600 ;% Tmax according to textbook
P04=P03;
% Turbo Fan START
Q=45000000; % KJ/KG
Cp1=1107; %KJ/KG 
f=(T04./T03-1)./(Q./(Cp1.*T03)-T04./T03);

Prf=1.5;
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1)).*(287.*T08).*(1-(Pa./P08).^((Gf-1)/Gf)) ).^(1/2);

T05=T04-(T03-T02)-B*(T08-T02);
P05=P04.*(1-(1/Nt).*(1-T05./T04)).^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1)).*(287.*T06).*(1-(Pa./P06).^((Gn-1)/Gn)) ).^(1/2);

ST=(1+f).*Ue + B*Uef-(1+B)*U;
STreal=(ST(3:60));
TSFC= 1000.*f./ST;
TSFC=TSFC(3:60);
PCR=(Pcr(3:60));
yyaxis left
semilogx(PCR,STreal/1000.,'r-');
yyaxis right
semilogx(PCR,TSFC,'r-')

T04=1700 ;% Tmax according to textbook
P04=P03;
Q=45000000; % KJ/KG
Cp1=1107; %KJ/KG 
f=(T04./T03-1)./(Q./(Cp1.*T03)-T04./T03);

Prf=1.5;
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1)).*(287.*T08).*(1-(Pa./P08).^((Gf-1)/Gf)) ).^(1/2);

T05=T04-(T03-T02)-B*(T08-T02);
P05=P04.*(1-(1/Nt).*(1-T05./T04)).^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1)).*(287.*T06).*(1-(Pa./P06).^((Gn-1)/Gn)) ).^(1/2);

ST=(1+f).*Ue + B*Uef-(1+B)*U;
STreal=(ST(2:70));
TSFC= 1000.*f./ST;
TSFC=TSFC(2:70);
PCR=(Pcr(2:70));
yyaxis left
semilogx(PCR,STreal/1000,'g-');
xlabel('Compressor Pressure Ratio')
ylabel('ST KNs/Kg')
yyaxis right
semilogx(PCR,TSFC,'g-')
xlim([2 100])
ylabel('Kg/Kns')
title('Turbofan Static thrust and Fuel Consumption')
hold off

%% M= 0.85

M0=0.85;
U= sqrt(1.4*287*T1) * M0;
figure
% Compressor inlet 
T02=T1*(1+(1.4-1)/2*M0^2) ;
P02=P1*(1+Nd*(T02/T1-1))^(Gd/(Gd-1));


Pcr= linspace(2,100); % Compressor pressure ratio (X-axis value) 
P03=P02.*Pcr ;
T03=T02.*(1+(1/Nc).*(Pcr.^((Gc-1)/Gc) -1) );
T04=1500 ;% Tmax according to textbook
P04=P03;
% Turbo Fan START
Q=45000000; % KJ/KG
Cp1=1107; %KJ/KG 
f=(T04./T03-1)./(Q./(Cp1.*T03)-T04./T03);

Prf=1.5;
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1)).*(287.*T08).*(1-(P1./P08).^((Gf-1)/Gf)) ).^(1/2);

T05=T04-(T03-T02)-B*(T08-T02);
P05=P04.*(1-(1/Nt).*(1-T05./T04)).^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1)).*(287.*T06).*(1-(P1./P06).^((Gn-1)/Gn)) ).^(1/2);

ST=(1+f).*Ue + B*Uef-(1+B)*U;
STreal=(ST(2:90));
TSFC= 1000.*f./ST;
TSFC=TSFC(2:90);
PCR=(Pcr(2:90));
semilogx(PCR,STreal/1000,'b')
xlim([2,100])
ylim([0,2])
yyaxis right
semilogx(PCR,TSFC,'b')
ylim([0 0.04])
hold on

T04=1600 ;% Tmax according to textbook
P04=P03;
% Turbo Fan START
Q=45000000; % KJ/KG
Cp1=1107; %KJ/KG 
f=(T04./T03-1)./(Q./(Cp1.*T03)-T04./T03);

Prf=1.5;
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1)).*(287.*T08).*(1-(P1./P08).^((Gf-1)/Gf)) ).^(1/2);

T05=T04-(T03-T02)-B*(T08-T02);
P05=P04.*(1-(1/Nt).*(1-T05./T04)).^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1)).*(287.*T06).*(1-(P1./P06).^((Gn-1)/Gn)) ).^(1/2);

ST=(1+f).*Ue + B*Uef-(1+B)*U;
STreal=(ST(2:100));
TSFC= 1000.*f./ST;
TSFC=TSFC(2:100);
PCR=(Pcr(2:100));
yyaxis left
semilogx(PCR,STreal/1000.,'r-');
ylabel("KNs/Kg")
yyaxis right
semilogx(PCR,TSFC,'r-')
xlabel('Compressor Pressure Ratio')
ylabel('Kg/KNs')
T04=1700 ;% Tmax according to textbook
P04=P03;
Q=45000000; % KJ/KG
Cp1=1107; %KJ/KG 
f=(T04./T03-1)./(Q./(Cp1.*T03)-T04./T03);

Prf=1.5;
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1)).*(287.*T08).*(1-(P1./P08).^((Gf-1)/Gf)) ).^(1/2);

T05=T04-(T03-T02)-B*(T08-T02);
P05=P04.*(1-(1/Nt).*(1-T05./T04)).^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1)).*(287.*T06).*(1-(P1./P06).^((Gn-1)/Gn)) ).^(1/2);

ST=(1+f).*Ue + B*Uef-(1+B)*U;
STreal=(ST(2:100));
TSFC= 1000.*f./ST;
TSFC=TSFC(2:100);
PCR=(Pcr(2:100));
yyaxis left
semilogx(PCR,STreal/1000,'g-');
ylabel('ST KNs/Kg')
yyaxis right
semilogx(PCR,TSFC,'g-')
ylim([0 0.06])
xlim([2 100])
xlabel('Compressor pressure ratio')
ylabel('Kg/KNs')
legend('1500', '1600','1700')
title('Turbofan cruise thrust and Fuel Consumption')
hold off

figure
Nth=((1+f).*(Ue.^2/2)-0.5*U^2)./(f.*Q);
semilogx(Pcr,Nth)
Np=2*U./(Ue+U)+0.15;
Nth=((1+f).*(Ue.^2).*0.5-U^2/2)./(f.*Q)+0.1;
No=Np.*Nth;
semilogx(Pcr,Np)
hold on

semilogx(Pcr,Nth)
hold on
semilogx(Pcr,No)
xlim([2 100])
ylim([0 1])
legend('N propulsion','N thermal','N overall')
title('Turbofan Efficiency vs Compressor p ratio')
xlabel('Compressor Pressure Ratio')