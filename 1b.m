%% Project of Automatic Control 
%%%%%%%%
%
%   Project:
%       Type 1.b
%   Team members:
%       Alessio Troffei, Vladyslav Tymofieiev,
%       Francesco Scavello
%
%%%%%%%%

% Project specs:
% - zero steady state error with a step reference signal of w(t)=W1(t)
% - Mf>45°
% - S_%<1%
% - T_a5< 0,2[s]
%% 
clear all, clc
%System parameters definition
b = 0.3; a = 0.2; F_v = -9; 
T_a5 = 0.2;
S_p = 0.01; % Sovraelongazione percentuale
w_n_min = 15 * 10^4;
w_n = 1.5 * 10^4;
A_n = 30; % abbattimento
mu_d=3;
W = 1;
%Matrixs
A = [0 , 1;
    a/4*F_v , -b];

B = [ 0 ;
      a ];
  
C = [1 , 0];
D = 0;

%Transfer function definition
s = tf("s");
[Num,Den] = ss2tf(A,B,C,D); %converts a state-space representation of a system into an equivalent transfer function.
display(Num);
display(Den);
figure();
G=tf(Num,Den);
%step(G);
%title("Step Response G");

zpk(G);
bode(G);
title("G");
%%
%Definition of bode frequency range
w_plot_min=10^(-2);%vedere e nel caso modificare
w_plot_max=10^6;%vedere e nel caso modificare

% Constraints variables
%xi=sqrt((log(S_p)^2)/(pi^2+log(S_p)^2)); %= circa 0.83;
xi=0.85;
Mf=xi*100;

w_c_min=300/(Mf*T_a5);
w_c_max=w_n; % [rad/s]
Xi_dB = 20*log10(1/A_n); %abbattimento di 30 volte

%Regolatore statico us/s^k con k>=0
%Rs = 6/s; % k=1; us=1;
Rs = 1/s;

%Ge = G* Rs
Ge = G * Rs;

zpk(Ge); 
bode(Ge);
title("Ge");
%% Plot Ge
figure();
[MagGe,phaseGe,wGe]=bode(Ge,{ w_plot_min, w_plot_max });
patch([w_n_min,w_plot_max,w_plot_max,w_n_min],[Xi_dB,Xi_dB,200,200],'r','FaceAlpha',0.4,'EdgeAlpha',0);
patch([w_c_min,w_plot_min,w_plot_min,w_c_min],[-300,-300,0,0],'g','FaceAlpha',0.4,'EdgeAlpha',0);
legend('disturbo misura', 'tempo di assest.', 'Ge');
hold on,  zoom on;
margin(MagGe, phaseGe, wGe);
patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-280,-280],'r','FaceAlpha',0.35,'EdgeAlpha',0);
%% REGOLATORE DINAMICO 
%% Uno zero
Rzero = 1 + 100*s;

G_zero = Ge * Rzero;

figure();
bode(G_zero);
title("Ge * zero");
%% Una rete anticipatrice 
omega_cStar = w_c_min+50;
M_star = Mf + 4;
phi_star = M_star - 180 -(-180);
tau= (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star));
a_tau = (cosd(phi_star)-(1/M_star))/(omega_cStar * sind(phi_star));
%tau1 = 1;
%tau2 = 100;
%a1 = 10^(-4);
%a2 = 10^(-5);
%Rd = (1+ tau1*s) * (1+tau2*s) / ((1+a2*tau2*s)*(1 + a1*tau1*s));
%%
R_ant = (1+tau*s) / (1+a_tau*s);
Rd = R_ant*Rzero;
L = Ge*Rd*mu_d;
%%
figure();
[MagGe,phaseGe,wGe]=bode(Ge,{ w_plot_min, w_plot_max });
[MagL,phaseL,wL]=bode(L,{ w_plot_min, w_plot_max });

patch([w_n_min,w_plot_max,w_plot_max,w_n_min],[Xi_dB,Xi_dB,200,200],[1 0.5 0],'FaceAlpha',0.4,'EdgeAlpha',0);
patch([w_c_min,w_plot_min,w_plot_min,w_c_min],[-300,-300,0,0],[1 1 0],'FaceAlpha',0.4,'EdgeAlpha',0);
leg1 = legend('1) disturbo misura', '2) tempo di assest.', '3)Ge');
hold on,  zoom on;

margin(MagGe, phaseGe, wGe);
margin(MagL, phaseL, wL);

leg2 = legend('Ge', 'L', 'M_f');

patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-280,-280],'r','FaceAlpha',0.35,'EdgeAlpha',0);
%%
T_simulation=5;
figure();
%step(Ge, 20, "b");
step(Ge);
title("Step response of Ge");
%% 
figure();
%step(L, 20, "b");
step(L);
title("Step response of L");
%%
figure();
%step(L, 20, "b");
step(L/(1+L), T_simulation, "b"); %F(S)

% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_p),W*(1+S_p),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on;

% add Settling time constraint
patch([T_a5,T_simulation,T_simulation,T_a5],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[W*(1+0.01),W*(1+0.01),W+1,W+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

title("Step response F(s)");
%%
if 1
W=-pi/3;
open("Progetto_1b")
% Define the controller coefficients to run the simulink file
R=Rd*Rs*mu_d;
[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};

display(num_r);
display(den_r);
%initial condition
x_e=[pi/3 0] ;
u_e = 3.8971;
x0=[0 0];
ye= pi/3;
end
return;