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
i = 1;
while a < 5
    Gc = K*G/(1+G*K);
    S=stepinfo(Gc);
    a = extractfield(S,'Overshoot');
    kvek(i)=K;
    K = K+0.0001;
    i=i+1;
end
K = kvek(i-2);
Gc = K*G/(1+G*K);
S=stepinfo(Gc)
G0 = K * G;
%% 3.3
clc
close all
clear all
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
figure()
bode(G)
Gc = G/(1+G)
S11=stepinfo(Gc)

[GGm, GPm, GWbredd, GWcross] = margin(G)
beta = 0.8;
Wcd = GWcross*2;
K = 1;
Td = 1/(Wcd*sqrt(beta));
Ti = 10/Wcd;
gamma = 0;
s = i*Wcd;
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)
Flead = K*(Td*s+1)/(beta*Td*s+1);
Flag = (Ti*s+1)/(Ti*s+gamma);
F = Flead*Flag;
K = 1/10^(abs(F*G)/20);
s = tf('s');
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)
Flead = K*(Td*s+1)/(beta*Td*s+1);
Flag = (Ti*s+1)/(Ti*s+gamma);
F = Flead*Flag;
Go = F*G;
Gc = Go/(1+Go);
S33=stepinfo(Gc)
figure()
bode(Go)
[FGm, FPm, FWbredd, FWcross] = margin(Go)
% [gf, fsa, qwe, mkg] = margin(system)
% Gc = feedback(K*Flead*G*Flag,1)
% figure(5)
figure()
step(Gc)

