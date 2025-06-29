clear
clc

%% Takeoff M=0 
%Define efficiencies
Nd=0.95;
Gd=1.4;
Npc=0.9;%compressor
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


syms M1 U T02 T03 T05 T1 P1 Pcr T04  Prf B ST TSFC Ue Uef A

U= sqrt(1.4*287*Ta) * M1;
T02=T1*(1+(1.4-1)/2*M1^2) ;
P02=P1*(1+Nd*(T02/T1-1))^(Gd/(Gd-1));
P03=P02*Pcr ;

mdot=231.8*2*(P02/P1)/((T02/Ta)^(0.5));
mdota=mdot/(1+B);

Nc=((Pcr)^((Gc-1)/Gc)-1)/((Pcr)^((Gc-1)/(Gc*Npc))-1);

T03=T02*(1+(1/Nc)*(Pcr^((Gc-1)/Gc) -1) );
P04=P03*0.95;
f=(T04/T03-1)/(Q/(Cp1*T03)-T04/T03);
P08=P02*Prf;
T08=T02*(1+1/Nf*(Prf^(0.4/1.4)-1));
Uef=( 2*Nf*(Gf/(Gf-1))*(287*T08)*(1-(P1/P08)^((Gf-1)/Gf)) )^(1/2);
T05=T04-(T03-T02)-B*(T08-T02);



P05=P04*(T05/T04)^(Gt/(Nt*(Gt-1)));
Thrustb=mdota*(((1+f)*Ue)+(B*Uef)-((1+B)*U));
Thrust=Thrustb/(1.04+(0.01*B^1.2));
T06=T05;
P06=P05;
Ue=( 2*Nn*(Gn/(Gn-1))*(287*T06)*(1-(P1/P06)^((Gn-1)/Gn)) )^(1/2);
% ST=(1+f)*Ue + B*Uef-(1+B)*U;
ST=Thrust/mdota;
TSFC= 1000*f/ST;
Np=2*U/(Ue+U);
Nth=((1+f)*(Ue^2)*0.5-U^2/2)/(f*Q);
No=Np*Nth;

for i=1:9

    for j=1:9 

        T04=1600;
        Pcr=20;
        Ta=217;

M1=1.7;
A=2;
        B=0.5+0.2*i ;
        Prf=1+0.05*j ;

        STmatrix(j,i)=subs(ST,{'B','M1','Prf','Ue','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
        TSFCm(j,i)=subs(TSFC,{'B','M1','Prf','Ue','Pcr','T04','T1','P1'},{B M1 Prf U Pcr T04 Ta P1});


    end
end
vpa(STmatrix,5)
  Prfz=[1.05:0.05:1.45]; 
     Bz=[0.7:.2:2.3];
 [x y]=meshgrid(Prfz,Bz);
 [mm nn]=contour(x,y,STmatrix/1000);
nn.ShowText='on' ;
nn.LineColor='r' ;
title('ST')
figure
[mm nn]=contour(x,y,TSFCm)
nn.ShowText='on'
nn.LineColor='b'
title('TSFC')
% % B=1.9
% % Fan pressure ratio = 1.25
% 
% for i=1:13
%     for j=1:14
% 
% B=1.9;
% Prf=1.25;
% Ta=217;
% T04=1475+25*i;
% Pcr=12+2*j;
% Nthmatrix(j,i)=subs(Nth,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
% Npmatrix(j,i)=subs(Np,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
% NovM(j,i)=subs(No,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
% 
%     end
% end
% figure
% T04z=[1500:25:1800];
% Pcrz=[14:2:40];
% [x y]=meshgrid(T04z,Pcrz);
% [mm nn]=contour(x,y,Nthmatrix);
% nn.ShowText='on'
% nn.LineColor='r'
% title('Nth')
% figure
% [mm nn]=contour(x,y,Npmatrix); 
% 
% 
% nn.ShowText='on'
% nn.LineColor='r'
% title('Np')
% hold off
% 
% figure 
% [mm nn]=contour(x,y,NovM);
% nn.ShowText='on';
% nn.LineColor='r' ;
% title('N overall')
% disp(figure)
% 
% %% Carpet Plot
% for i=1:12
%     for j=1:12
% 
% B=1.9;
% Prf=1.25;
% Ta=217;
% T04=1450+50*i;
% Pcr=12+2*j;
% 
% 
%         STmatrix(i,j)=subs(ST,{'B','M0','Prf','U','Pcr','T04','T1', 'P1'}, {B M1 Prf U Pcr T04 Ta P1});
%         TSFCm(j,i)=subs(TSFC,{'B','M0','Prf','U','Pcr','T04','T1','P1'},{B M1 Prf U Pcr T04 Ta P1});
% 
% 
%     end
% 
% end
% vpa(STmatrix,5);
% for j=1:12; 
%     Str=STmatrix(j,:);
%     Tsfcr=TSFCm(j,:);
%     plot(Str/1000,Tsfcr)
%     hold on
%     Str=STmatrix(:,j);
%     Tsfcr=TSFCm(:,j);
%     plot(Str/1000,Tsfcr)
% 
% end
