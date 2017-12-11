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
lab3robot(G,K,[],[],[],[],[],[],960703)
%% 3.3
clc
close all
clear all
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
bode(G), grid
Gc = G/(1+G)
S11=stepinfo(Gc)
[GGm, GPm, GWbredd, GWcross] = margin(G)
%% kims test
close all
%skapar G igen.
[J,umax] = lab3robot(960703);
s = tf('s')
kt=38;
Lm=2;
km=0.5;
n=1/20;
Rm=21;
b=1;
wcd = 4*0.0505;

%beräknar Flead och Flag
s = i*wcd;
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)
%k = 1/abs(G)
beta = 0.3;
Td = 1/(wcd*sqrt(beta))
Flead = (Td*s+1)/(beta*Td*s+1);
k = 1/abs(Flead*G)
gamma = 0.6;
Ti = 10/wcd+10;
s = tf('s')
Flead = (Td*s+1)/(beta*Td*s+1);
Flag = (Ti*s+1)/(Ti*s+gamma);
Kg = kt/(s*Lm+Rm)
T=1/(J*s+b)
G = (Kg*T*n/s)/(km*Kg*T+1)

%konstruerar nya Gc och kontrollerar resultat samt plottar bode av Go
Gc= feedback(k*Flead*Flag*G/(1+k *Flead*Flag*G),1)
figure(66)
bode(k * Flead * Flag * G), grid
 figure(5)
 step(Gc)
S33=stepinfo(Gc)
lab3robot(G,k,Flead*Flag,[],[],[],[],[],960703) % i den översta delen har vi ett annat K,
%k är väll en del av lead/lag, vi får ju inte rätt på proportionallitets
%kontrollern här men det får vi i koden ovan.
lab3robot(G,960703)
