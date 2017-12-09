clc;
clear all;
close all;
%EL101 5.13 9.14ab 6.10
%% 5.13 
s=tf('s');
Go = 725/((s+1)*(s+2.5)*(s+25))
figure()
margin(Go)
[Gm Pm Wc Wp]=margin(Go) %5.13 a) Am=10.7607, Pm = 26.6244, Wc=9.4868, Wp=4.9934
Gm_db=20*log10(Gm)
%rest of 5.13
beta = 0.2899;
wc = 5;
taud =1/(wc*sqrt(beta));
s = i*5;
Flead = (taud*s+1)/(beta*taud*s+1);
Go = 725/((s+1)*(s+2.5)*(s+25));
K=1/(abs(Flead)*abs(Go));
s=tf('s');
Flead = K*(taud*s+1)/(beta*taud*s+1);
Go = 725/((s+1)*(s+2.5)*(s+25));
G=Flead*Go;
figure()
margin(G);
[Gmb Pmb Wcb Wpb]=margin(G)
figure()
step(Go)
hold on
step(G)
legend('Go','G')
figure()
bode(1/(1+Go))
figure()
bode(1/(1+G))
t = 0:0.1:1
ref = t
lsim(G,ref,t)

%% uppgift 6.10
clc;
clear all;
close all;
