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
bode(G)
s = tf('s')
[Gm, Pm, Wbredd, Wcross] = margin(G)
beta = 0.7;
Wcd = 0.95;

Td = 1/(Wcd*sqrt(beta));

Flead = K*(Td*s+1)/(beta*Td+1);

%F = Flead*Flag;
Go = Flead*G;
Gc = Go/(1+Go);
S33=stepinfo(Gc)
bode(Go)
[Gm, Pm, Wbredd, Wcross] = margin(Go)
% [gf, fsa, qwe, mkg] = margin(system)
% Gc = feedback(K*Flead*G*Flag,1)
% figure(5)
% step(Gc)
