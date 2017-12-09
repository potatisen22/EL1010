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
K = 1.05;
step(K*Gc)
S=stepinfo(K*Gc);
lab3robot(G,K,960703)