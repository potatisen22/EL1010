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
Gc = G/(1+G);
i = 1;
<<<<<<< HEAD
K = 1.05;
step(K*Gc)
S=stepinfo(K*Gc);
lab3robot(G,K,960703)
=======
k = 1.05;
step(k*Gc)
S=stepinfo(k*Gc);s


%lab3robot(G,960703)

[Gm, Pm, Wcg, Wcp] = margin(G)


% bode(system)
% [gf, fsa, qwe, mkg] = margin(system)

% Gc = feedback(K*Flead*G*Flag,1)
% figure(5)
% step(Gc)
>>>>>>> 43d8df7293492bd0ca5a34d6ba528de9cd82282a
