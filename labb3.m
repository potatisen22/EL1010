%% 3.1
clc;
clear all;
close all;
[J,umax] = lab3robot(960703);
s = tf('s')
kt=38;
Lm=2;
km=0.5;
n=1/20;
Rm=21;
b=1;
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)
lab3robot(G,960703)


K = 5;
Gc = K*G/(1+G*K);
%i = 1;
step(Gc)
S=stepinfo(Gc);
%lab3robot(G,K,960703)
%k = 1.05;
%step(k*Gc)
%S=stepinfo(k*Gc);s
K = 5.12;
a = 0;
l = 1;
while a < 5
    Gc = K*G/(1+G*K);
    S=stepinfo(Gc);
    a = extractfield(S,'Overshoot');
    kvek(l)=K;
    K = K+0.0001;
    l=l+1;
end
K = kvek(l-2);
Gc = K*G/(1+G*K);
S=stepinfo(Gc)
G0 = K * G;
lab3robot(G,K,[],[],[],[],[],[],960703)
Kp = K;

%% 3.3
clc
close all
[J,umax] = lab3robot(960703);
s = tf('s')
%/definierar G
kt=38;
Lm=2;
km=0.5;
n=1/20;
Rm=21;
b=1;
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)
%/
%%Okompenserat
figure()
bode(Kp*G), grid
title('bode med Kp')
figure()
bode(G), grid
title('bode för G')
Gc = G/(1+G)
S11=stepinfo(Gc)
[GGm, GPm, GWbredd, GWcross] = margin(G)
%% kims test
close all
[J,umax] = lab3robot(960703);
s = tf('s')
kt=38;
Lm=2;
km=0.5;
n=1/20;
Rm=21;
b=1;
wcd = 4*0.222;
s = i*wcd;
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)

beta = 0.17;
Td = 1/(wcd*sqrt(beta))
Flead = (Td*s+1)/(beta*Td*s+1);
k = 1/abs(Flead*G)  %% övning 8?
gamma = 0.0382;
Ti = 10/wcd+10;
s = tf('s')
Flead = k*(Td*s+1)/(beta*Td*s+1);
Flag = (Ti*s+1)/(Ti*s+gamma);
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)

Gc = Flead*Flag*G/(1+Flead*Flag*G); 
figure(66)
bode(Flead * Flag * G), grid
 figure(5)
 step(Gc)
S33=stepinfo(Gc)

lab3robot(G,Kp,Flead*Flag,[],[],[],[],[],960703)
lab3robot(G,960703)
%% sensitivity
S1 = 1/(1+(G*Flead*Flag))
S2 = 1/(1+Kp*G)
figure()
bodemag(S1,S2), grid
legend('LeadLag','Kp')
%% 3.5
clc;
close all;
F=Flead*Flag;
T = G*F/(1+G*F); 
dG1 = (s+10)/40;;
dG2 =(s+10)/(4*(s+0.01));
bodemag(1/T,dG1,dG2), grid
legend('1/T','dG1','dG2')
%% 3.6
clc;
close all;
A = [0 n 0; 0 -b/J kt/J; 0 -km/Lm -Rm/Lm];
B = [0; 0; 1/Lm];
C = [1 0 0];
S = [B A*B (A^2)*B];
O = [C;C*A;C*A^2];
lab3robot(G,Kp,Flead*Flag,A,B,C,[],[],960703)
det(S) %IT IS CONTROLLABlE SINCE THIS VALUE IS NOT 0
det(O) % IT IS OBSERVABLE SINCE THIS VALUE IS NOT 0
%% 3.7
clc

close all;
%beräkningar för att ta fram funktioner för l1,l2,l3 av polerna nedan
% syms s l1 l2 l3;
I = eye(3);
 % L = [l1, l2, l3];
 % det(s*I-A+B*L)
%poler ap+-i*bp och cp samt tillverkning av Gc nedan
s = tf('s');
ap = 1; %reella polen vi fuckade upp polerna, så de måste anges med positiv realdel
bp = 0; %im del av ap
cp = 1; % reella poolen
l1 =cp*(ap^2+bp^2)*J*Lm/(kt*n)
l3 = Lm*(2*ap+cp-Rm/Lm-b/J);
l2 = Lm*J*(ap^2+bp^2+2*ap*cp-Rm*b/(J*Lm)-b*l3/(J*Lm)-km*kt/(J*Lm))/kt;
L = [l1, l2, l3];
L0 = l1;
afb = A-B*L;
bfb = B*L0;
D = 0; %inget beroende av U
Gcfeedback = tf(ss(afb,bfb,C,D))  %C*(s*I-A+B*L)^-1 * B;
stepinfo(Gcfeedback)
step(Gcfeedback)
lab3robot(G,Kp,Flead*Flag,A,B,C,L,L0,960703)
u = Gcfeedback/G
close all;
figure()
step(u)
hold on
%Gcl = (Gcfeedback/(1+Gcfeedback))*(1/Gbild)
step(u)
hold on 
plot(refline([0 110]),refline([0 -110]))
%max(step(u))