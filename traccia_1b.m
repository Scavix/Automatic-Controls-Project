%% Progetto di Controlli Automatici
%%%%%%%%
%
%   Progetto: traccia 1.b
%   Membri del gruppo: Alessio Troffei, Vladyslav Tymofieiev,
%   Francesco Scavello
%
%%%%%%%%

% Specifice di progetto
% - zero steady state error with a step reference signal of w(t)=W1(t)
% - Mf>45°
% - S_%<1%
% - T_a5< 0,2[s]

%% 
clear all  clc
% Definizioni variabili del sistema di controllo
b = 0.3; a = 0.2; F_v = -9; % Coefficienti sistema
T_a5 = 0.2; % Tempo di assestamento
S_p = 0.01; % Sovraelongazione percentuale
w_n_min = 15 * 10^4;
w_n = 1.5 * 10^4;
A_n = 30; % Abbattimento
mu_d=3; % Per regolatore dinamico
W=-pi/3;

% Matrici
A = [0 , 1;
    a/4*F_v , -b];

B = [ 0 ;
      a ];
  
C = [1 , 0];
D = 0;

% Definizione funzione di trasferimento
s = tf("s");

% Converte una rappresentazione dello spazio di stato di un sistema in una funzione di trasferimento equivalente
[Num,Den] = ss2tf(A,B,C,D);

% Stampa sulla console
display(Num); 
display(Den);
figure();

% Ottengo la funzione di trasferimento del sistema in forma N(s)/D(s)
G=tf(Num,Den);
%step(G);
%title("Step Response G");

zpk(G);
bode(G);
title("G");
%%
% Definizione frequenza del Bode
w_plot_min=10^(-2); % vedere e nel caso modificare
w_plot_max=10^6; % vedere e nel caso modificare

% Variabili dei vincoli
%xi=sqrt((log(S_p)^2)/(pi^2+log(S_p)^2)); %= circa 0.83;
xi=0.85;
Mf=xi*100;

w_c_min=300/(Mf*T_a5);
w_c_max=w_n; % [rad/s]
Xi_dB = 20*log10(1/A_n); %abbattimento di 30 volte

% Regolatore statico mu_s/s^k
% Uso questo regolatore statico, per andare a soddisfare la specifica
% sull'errore a regime nullo
%Rs = 6/s; % k=1; us=1;
% Metto la mu_s = 1 per potere lasciare libera la mu_d del regolatore dinamico
Rs = 1/s;


%
% !!! Per ora ho 1 solo polo !!!
%


% La Ge è la funzione di trasferimento estesa
% è ottenuta moltiplicando le funzioni di trasferimento del sistema di 
% controllo (G) e del regolatore statico (Rs)
Ge = G * Rs;

zpk(Ge); 
bode(Ge);
title("Ge");

%% Grafico di Ge
figure();
[MagGe,phaseGe,wGe]=bode(Ge,{ w_plot_min, w_plot_max });
patch([w_n_min,w_plot_max,w_plot_max,w_n_min],[Xi_dB,Xi_dB,200,200],'r','FaceAlpha',0.4,'EdgeAlpha',0);
patch([w_c_min,w_plot_min,w_plot_min,w_c_min],[-300,-300,0,0],'g','FaceAlpha',0.4,'EdgeAlpha',0);
legend('disturbo misura', 'tempo di assest.', 'Ge');
hold on,  zoom on;
margin(MagGe, phaseGe, wGe);
patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-280,-280],'r','FaceAlpha',0.35,'EdgeAlpha',0);

%
% Si nota che devo recuperare una fase maggiore di 90*
% devo avere due zeri o due reti anticipatrici!
%

%% REGOLATORE DINAMICO 
% Questa funzione di trasferimento parziale contiene uno zero
Rzero = 1 + 100*s;

%
% !!! Per ora 1 polo e 1 zero !!!
%

% Questa funzione di trasferimento più evoluta, è stata ottenuta
% moltiplicando la funzione di trasferimento estesa (Ge) per quella dello
% zero scritta sopra
G_zero = Ge * Rzero;

figure();
bode(G_zero);
title("Ge * zero");

%% Una rete anticipatrice 
omega_cStar = w_c_min+50; % Indica la pulsazione di attraversamento desiderata
M_star = Mf + 4; % Aggiungo 4 (perché???)
phi_star = M_star - 180 -(-180); % Indica il margine di fase 
tau = (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star)); % Formula di inversione
a_tau = (cosd(phi_star)-(1/M_star))/(omega_cStar * sind(phi_star)); % Formula di inversione
%tau1 = 1;
%tau2 = 100;
%a1 = 10^(-4);
%a2 = 10^(-5);
%Rd = (1+ tau1*s) * (1+tau2*s) / ((1+a2*tau2*s)*(1 + a1*tau1*s));
%%

% Formula rete anticipatrice
R_ant = (1+tau*s) / (1+a_tau*s); 

%
% !!! Infine, ho 2 poli e 2 zeri !!!
% !!!  GRADO RELATIVO RISPETTATO !!!
%

% Regolatore dinamico, ottenuto andando a moltiplicare
% mu_s * Rete anticipatrice * Zero
%
%
% Il regolatore dinamico mi serve per andare a:
% 1) Imporre wc in un certo range (wcMin, wcMax) dettato dalle specifiche
% 2) Garantire che Mf >= MfMin
% 3) Avere una certa attenuazione e pendenza di L(s) a pulsazioni elevate
%    (w >> wc) così da garantire specifiche su n(t) e fisica realizzabilità
Rd = mu_d*R_ant*Rzero;

% Vado a creare la funzione di trasferimento in anello aperto (L),
% moltiplicando la mia Ge per il regolatore dinamico
L = Ge*Rd;
%%

figure();
[MagGe,phaseGe,wGe]=bode(Ge,{ w_plot_min, w_plot_max });
[MagL,phaseL,wL]=bode(L,{ w_plot_min, w_plot_max });

patch([w_n_min,w_plot_max,w_plot_max,w_n_min],[Xi_dB,Xi_dB,200,200],[1 0.5 0],'FaceAlpha',0.4,'EdgeAlpha',0);
patch([w_c_min,w_plot_min,w_plot_min,w_c_min],[-300,-300,0,0],[1 1 0],'FaceAlpha',0.4,'EdgeAlpha',0);
leg1 = legend('1) disturbo misura', '2) tempo di assest.', '3)Ge');
hold on,  zoom on;

margin(MagGe, phaseGe, wGe); % ????
margin(MagL, phaseL, wL); % ????

leg2 = legend('Ge', 'L', 'M_f');

patch([w_c_min,w_c_max,w_c_max,w_c_min],[-180+Mf,-180+Mf,-280,-280],'r','FaceAlpha',0.35,'EdgeAlpha',0);

%% Risposta al gradino di Ge
T_simulation=5; % ????
figure();
%step(Ge, 20, "b");
step(Ge); % Risposta al gradino di Ge
title("Risposta al gradino di Ge");

%% Risposta al gradino di L
figure();
%step(L, 20, "b");
step(L); % Risposta al gradino di L
title("Risposta al gradino di L");

%%
figure();
%step(L, 20, "b");
[y_step,t_step] = step(L/(1+L), T_simulation, "b"); %F(S)
plot(t_step,y_step)

% Aggiunta vincolo sulla sovraelongazione
patch([0,T_a5,T_a5,0],[(1+S_p),(1+S_p),2,2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on;

% Aggiunta del vincolo sul tempo di assestamento
patch([T_a5,T_simulation,T_simulation,T_a5],[(1-0.01),(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5,T_simulation,T_simulation,T_a5],[(1+0.01),(1+0.01),2,2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

title("F(s) Step Response");
%%

% Mettere 0 per non far avviare Simulink
if 1

open("Progetto_1b")

% Definizione dei coefficienti del controllore per eseguire il file Simulink
R=Rd*Rs*mu_d;
[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};

display(num_r);
display(den_r);

% Condizioni iniziali
x_e=[pi/3 0] ;
u_e = 3.8971;
x0=[0 0];
ye= pi/3;
end
return;
%%
