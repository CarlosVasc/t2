close all
clear all

pkg load symbolic
pkg load control
pkg load signal


%tirar dados do ficheiro.txt produzido pelo make
format longg
%abertura do ficheiro
fid=fopen("testi.txt","r");
i=1;
while ~feof(fid)
     linha=fgetl(fid);
   if(linha)~=0 
       if i>6 
        [dado,resto]=strtok(linha," ");
        [lixo,valor]=strtok(resto," ");
        d(i)=str2double(valor);
        
       elseif i==6
        [dado,resto]=strtok(linha," R1 = ");
        [lixoaux,valoraux]=strtok(resto," ");  
        [lixo,valor]=strtok(valoraux," ");
        d(i)=str2double(valor);
        
       end
   else 
       d(i)=0; 
   end
   
    i=i+1;
end
fclose(fid);

R1=d(6);
R2=d(7);
R3=d(8);
R4=d(9);
R5=d(10);
R6=d(11);
R7=d(12);
Vs=d(13);
C=d(14);
Kb=d(15);
Kd=d(16);

G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

%----------  1  ----------------
A = [    0   0   1   0   0   0   0   0   0   0   0   0;
        -1   0   0   1   0   -1  0   0   0   0   0   0; 
        Kb   0   0   0   0   0   0   0   0   -1  0   0;
        0   -1   0   0   0   0   0   0   0   0   0   Kd;
        0   -1   0   0   0   1   0   0  -1   0   0   0;
        0    0   0   0   0   0   0 -G6   0   0   0  -1;
        0    0   0   0   0   0   0   0   0   0   1   0;
        0    0  -G1 G1+G2+G3   -G2   -G3   0   0       0   0    0   0;
        0   0   0  -G2-Kb   G2   Kb      0     0       0   0   0   0;
        0   0   0    Kb     0   -Kb-G5   G5    0       0   0   0   0;
        0   0   0    0      0    0       0   G7+G6   -G7   0   0   0; 
        0   0   G1  -G1     0   -G4      0     0       0   0   0   1];
    
B = [Vs;0;0;0;0;0;0;0;0;0;0;0];

NC=linsolve(A,B);


Vb= NC(1);
Vc= NC(2);
V1= NC(3);
V2= NC(4);
V3= NC(5);
V5= NC(6);
V6= NC(7);
V7= NC(8);
V8= NC(9);
Ib= NC(10);
Ic= NC(11);
Id= NC(12);

%Criar ficheiro octave teorico
fid11=fopen("teorico1.tex","w+");
fprintf(fid11,"V1 & %f\\\\ \\hline \n\
V2 & %f\\\\ \\hline \n\
V3 & %f\\\\ \\hline \n\
V5 & %f\\\\ \\hline \n\
V6 & %f\\\\ \\hline \n\
V7 & %f\\\\ \\hline \n\
V8 & %f\\\\ \\hline \n\
Id & %f\\\\ \\hline", V1,V2,V3,V5,V6,V7,V8,Id)
fclose(fid11)



fid12=fopen("ngspice_1.tex","w");
fprintf(fid12,"Vs V1 0 DC %.11f \n\
R1 V2 V1 %.11f \n\
R2 V3 V2 %.11f \n\
R3 V2 V5 %.11f \n\
R4 0 V5 %.11f \n\
R5 V6 V5 %.11f \n\
R6 V9 V7 %.11f \n\
R7 V7 V8 %.11f \n\
VVc 0 V9 0V \n\
HVc V5 V8 VVc %.11f \n\
GIb V6 V3 V2 V5 %.11f \n\
C1 V6 V8 %.11f", Vs, R1, R2, R3, R4, R5, R6, R7, Kd, Kb, C);


fclose(fid12);

%------  2  -------------------

Vx = V6 - V8;

A2 = [-G1 G1+G2+G3 -G2 -G3 0 0 0 0 0 0 0 ;
        0 -G2 G2 0 0 0 0 -1 0 0 0 ; 
        0 0 0 0 0 G7 -G7 0 -1 0 0 ;
        0 0 0 0 1 0 -1 0 0 0 0 ;
        0 -1 0 1 0 0 0 0 0 1 0 ;
        0 0 0 -1 0 0 1 0 0 0 1 ;
        0 0 0 0 0 0 0 1 0 -Kb 0;
        0 0 0 0 0 0 0 0 -Kd 0 1;
        G1+G4 -G1 0 -G4 0 0 0 0 1 0 0;
        1 0 0 0 0 0 0 0 0 0 0 ;
        -G6 0 0 0 0 G6 0 0 1 0 0];
        
B2 = [0;0;0;Vx;0;0;0;0;0;0;0];

NC2=linsolve(A2,B2);

V1_2 = NC2(1);
V2_2 = NC2(2);
V3_2 = NC2(3);
V5_2 = NC2(4);
V6_2 = NC2(5);
V7_2 = NC2(6);
V8_2 = NC2(7);
Ib_2 = NC2(8);
Id_2 = NC2(9);
Vb_2 = NC2(10);
Vd_2 = NC2(11);





Ix = Ib_2 + (V6_2 - V5_2)*G5;
Req = Vx/Ix;
tau = Req * C;


%Criar ficheiro octave teorico
fid21=fopen("teorico2.tex","w+");
fprintf(fid21,"Vx & %f\\\\ \\hline \n\
Ix & %f\\\\ \\hline \n\
Req & %f\\\\ \\hline \n\
Tau & %f\\\\ \\hline ", Vx,Ix,Req,tau)
fclose(fid21)



%ficheiro para ngspice 2
fid22 = fopen("ngspice_2.txt", "w");
fprintf(fid22, "Vs V1 0 DC 0 \n\
R1 V2 V1 %.11f \n\
R2 V3 V2 %.11f \n\
R3 V2 V5 %.11f \n\
R4 0 V5 %.11f \n\
R5 V6 V5 %.11f \n\
R6 V9 V7 %.11f \n\
R7 V7 V8 %.11f \n\
VVd 0 V9 0V \n\
HVd V5 V8 VVd %.11f \n\
GIb V6 V3 V2 V5 %.11f \n\
Vx V6 V8 DC %.11f", R1, R2, R3, R4, R5, R6, R7, Kd, Kb,Vx); 
fclose(fid22);

%------------------------  3  ------------------------


syms t
syms v6(t)
t=0:1e-6:20e-3;
v6 = (V8_2 + Vx) * exp(-(t/tau));

Teoria_3 = figure ();
plot(t*1000, v6);
xlabel ("t[ms]");
ylabel ("V_{6n}(t) [V]");
print (Teoria_3, "Teoria_3_Fig.eps");


fid3 = fopen("ngspice_3.txt", "w");
fprintf(fid3, "Vs V1 0 DC 0\n\ 
R1 V2 V1 %.11e \n\
R2 V3 V2 %.11e \n\
R3 V2 V5 %.11e \n\
R4 0 V5 %.11e \n\
R5 V6 V5 %.11e \n\
R6 V9 V7 %.11e \n\
R7 V7 V8 %.11e \n\
VVd 0 V9 0V \n\
HVd V5 V8 VVd %.11e \n\
GIb V6 V3 V2 V5 %.11e \n\
C1 V6 V8 %.11e ic = %.11e \n\
.ic v(V6) = %.11e v(V8) = 0", R1, R2, R3, R4, R5, R6, R7, Kd, Kb,C,Vx, Vx); 
fclose(fid3);

%----------4--------------






