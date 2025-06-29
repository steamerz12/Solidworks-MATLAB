clear
clc

%% Takeoff M=0 
%Define efficiencies
Nd=0.95;
Gd=1.4;
Nc=0.9;%compressor
Gc=1.37;
Nf=0.92; %Fan efficiency
Gf=1.4;
Nb=0.97;
Pb=0.95; %burner pressure recovery factor
Gb=1.35;
Nt=0.92;%turbine
Gt=1.33;
Nn=0.98;%core exit nozzle
Gn=1.36;
Nfn=0.99; %Fan nozzle efficiency
Gfn=1.4;
Q=45000000;
Pa=101300;
%ambient P and T
P1=7170;
Ta=288.2;
T1=216;
M0=0;
M1=1.7;
Cp1=287*Gb/(Gb-1);
syms M0 U T02 T03 T05 T1 P1 Pcr T04  Prf B ST TSFC Ue Uef f
U= sqrt(1.4*287*Ta) * M0 ;
T02=T1*(1+(1.4-1)/2*M0^2) ;
P02=P1*(1+Nd*(T02/T1-1))^(Gd/(Gd-1));
P03=P02*Pcr ;
T03=T02*(1+(1/Nc)*(Pcr^((Gc-1)/Gc) -1) );
P04=P03*0.95;
f=(T04/T03-1)/(Q/(Cp1*T03)-T04/T03);
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1))*(287*T08)*(1-(P1/P08)^((Gf-1)/Gf)) )^(1/2);
T05=T04-(T03-T02)-B*(T08-T02);
P05=P04*(1-(1/Nt)*(1-T05/T04))^(Gt/(Gt-1));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1))*(287*T06)*(1-(P1/P06)^((Gn-1)/Gn)) )^(1/2);
ST=(1+f)*Ue + B*Uef-(1+B)*U;
vpa(ST,5);
TSFC= 1000*f/ST;
Np=2*U/(Ue+U);
Nth=((1+f)*(Ue^2)*0.5-U^2/2)/(f*Q);
No=Np*Nth;

for i=1:9

    for j=1:9 

        T04=1600;
        Pcr=20;
        Ta=217;

        B=0.5+0.2*i ;
        Prf=1+0.05*j ;

        STmatrix(j,i)=subs(ST,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
        TSFCm(j,i)=subs(TSFC,{'B','M0','Prf','U','Pcr','T04','T1','P1'},{B M1 Prf U Pcr T04 Ta P1});


    end
end
figure
  Prfz=[1.05:0.05:1.45]; 
     Bz=[0.7:.2:2.3];
 [x y]=meshgrid(Prfz,Bz);
 [mm nn]=contour(x,y,STmatrix/1000);
nn.ShowText='on';
nn.LineColor='r';
title('ST')
xlabel('Fan pressure ratio')
ylabel('Bypass Ratio')
figure
[mm nn]=contour(x,y,TSFCm);
nn.ShowText='on';
nn.LineColor='b';
title('TSFC')
xlabel('Fan pressure ratio')
ylabel('Bypass Ratio')
% B=1.9
% Fan pressure ratio = 1.25

for i=1:13
    for j=1:14

B=1.9;
Prf=1.25;
Ta=217;
T04=1475+25*i;
Pcr=12+2*j;
Nthmatrix(j,i)=subs(Nth,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
Npmatrix(j,i)=subs(Np,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
NovM(j,i)=subs(No,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});

    end
end
figure
T04z=[1500:25:1800];
Pcrz=[14:2:40];
[x y]=meshgrid(T04z,Pcrz);
[mm nn]=contour(x,y,Nthmatrix);
nn.ShowText='on';
nn.LineColor='r';
title('Nth')
xlabel('Inlet Stagnation Temp')
ylabel('Compressor Pressure Ratio')
figure
[mm nn]=contour(x,y,Npmatrix); 


nn.ShowText='on';
nn.LineColor='r';
title('Np')
xlabel('Inlet Stagnation Temp')
ylabel('Compressor Pressure Ratio')
hold off

figure 
[mm nn]=contour(x,y,NovM);
nn.ShowText='on';
nn.LineColor='r' ;
title('N overall')
xlabel('Inlet Stagnation Temp')
ylabel('Compressor Pressure Ratio')
disp(figure)

%% Carpet Plot
for i=1:13
    for j=1:13

B=1.9;
Prf=1.25;
Ta=217;
T04=1450+50*i;
Pcr=12+2*j;


        STmatrix(i,j)=subs(ST,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
        TSFCm(j,i)=subs(TSFC,{'B','M0','Prf','U','Pcr','T04','T1','P1'},{B M1 Prf U Pcr T04 Ta P1});


    end

end
vpa(STmatrix,5);
for j=1:13; 
    Str=STmatrix(j,:);
    Tsfcr=TSFCm(j,:);
    plot(Str/1000,Tsfcr);
    hold on
    Str=STmatrix(:,j);
    Tsfcr=TSFCm(:,j);
    plot(Str/1000,Tsfcr);

    
end

%min ST = 0.44
