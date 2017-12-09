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
K = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (K*T*n/s)/(km*K*T+1)
lab3robot(G,960703)
Gc = G/(1+G);
i = 1;
k = 1.05;
step(k*Gc)
S=stepinfo(k*Gc);s