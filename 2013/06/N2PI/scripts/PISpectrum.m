% This program calculates how a PI spectrum should look like on
% basis of the molecular parameters and a laser linewidth
% Changes made by N Seymour-Smith 05/06/13

clear all;

% fundamental constants

c = 3e8;
kB = 1.380e-23;
h = 6.626*10^(-34);
T = 300;

% Ground state parameters

%BX0 = 1.9895;
DX0 = 5.8E-6; %Literature value
%DX0 = 0;
% wXe = 2358.569;
% wXex = 14.324;

% a1Pig state parameter

%Te = 10521;
% Ea0 = 84371.61; %Matthias' value
Ea0 = 42187.33; %From paper (calc.)
%Ea0 = 42187.49; %From Jack measurement
%Ba0 = 1.47997;
%Ba0 = 1.493; %Matthias' value
Ba0 = 1.423; %From paper (calc.)
%BX0 = Ba0+0.563032; %Matthias' value
BX0 = 0.9788; %from paper (calc.)
%BX0 = 1;
Da0 = 5.8E-6;
%Da0 = 0;
% wae = 1530.307;
% waex = 12.045;
% waey = 0.0351;
% alphab = 0.01656;
% betab = 0.00007;

% vibrational level
% NX = 0;
% Na = 10;
% EX0 = wXe*(NX+1/2)- wXex*(NX+1/2)^2;
% EX0 = 0;
% Ea0 = Te + wae*(Na+1/2) - waex*(Na+1/2)^2 + waey*(Na+1/2)^3;
% Ea0 = Te;
% Ba = Ba0-alphab*(Na+1/2)+betab*(Na+1/2)^2;

Ba = Ba0;
% Measurements

M = [236.99019 237.00868 237.022 237.03708 237.038 237.05605 237.07156 237.078 237.08855 237.104 237.11, 237.094, 237.049]';


% state energies

% number of rot states
N = 5;

EX = zeros(1,N);
Ea = zeros(1,N);
P = zeros(1,N);

for i = 1:N
    j = i-1;
    EX(i) = BX0*j*(j+1) - DX0*(i*(i+1))^2;
    Ea(i) = Ba*j*(j+1) + Ea0 - Da0*(j*(j+1))^2 ;
    P(i)= (2*j+1)*exp(-100*EX(i)*c*h/(kB*T));
    if mod(j,2) == 1
        P(i)= P(i)/2;
    end,
end,

% lambdaS = zeros(N-2,1);
% lambdaR = zeros(N-2,1);
% lambdaQ = zeros(N-2,1);
% lambdaP = zeros(N-2,1);
% lambdaO = zeros(N-2,1);
Pop = zeros(N-2,1);

for i = 1:N-3
    lambdaS(i)= 0.01*1e9/(Ea(i+1)-EX(i));
    lambdaR(i)= 0.01*1e9/(Ea(i)-EX(i));
    lambdaQ(i)= 0.01*1e9/(Ea(i)-EX(i+1));
    lambdaP(i)= 0.01*1e9/(Ea(i)-EX(i+2));
    lambdaO(i)= 0.01*1e9/(Ea(i)-EX(i+3));
    Pop(i)=P(i);
%     if i > 1
%         lambdaP(i)= 0.02/(Ea(i-1)-EX(i));
%     end,
%     if i > 2
%         lambdaO(i)= 0.02/(Ea(i-2)-EX(i));
%     end,
end

%Spectrum

width = 0.000722;
points = 10000;
signalS = zeros(1,points);
signalR = zeros(1,points);
signalQ = zeros(1,points);
signalP = zeros(1,points);
signalO = zeros(1,points);
lambda0 = min(lambdaS)-8*width;
dlambda = (max(lambdaO)-lambda0+8*width)/points;

for i=1:points
      lambda(i) = lambda0+dlambda*i;
      for j=1:N-3
        signalS(i)= signalS(i)+Pop(j)*width^2/((lambda(i)-lambdaS(j))^2+width^2/4);
        signalR(i)= signalR(i)+Pop(j)*width^2/((lambda(i)-lambdaR(j))^2+width^2/4);
        signalQ(i)= signalQ(i)+Pop(j)*width^2/((lambda(i)-lambdaQ(j))^2+width^2/4);
        signalP(i)= signalP(i)+Pop(j)*width^2/((lambda(i)-lambdaP(j))^2+width^2/4);
        signalO(i)= signalO(i)+Pop(j)*width^2/((lambda(i)-lambdaO(j))^2+width^2/4);
    end,
end,

pl = max(Pop)/2*ones(length(M),1);
plT = 0.99*ones(length(lambdaS),1);
plR = 0.99*ones(length(lambdaR),1);
plP = 0.99*ones(length(lambdaP),1);

%signalT = signalS + signalR + signalQ + signalP + signalO;

figure(1)
plot(lambda, signalS, lambda, signalR, lambda, signalQ, lambda, signalP, lambda, signalO, M', pl, '*')
xlim([min(M)-0.02 237.15])
legend('S-Branch', 'R-Branch', 'Q-Branch', 'P-Branch', 'O-Branch', 'Measurement');

% figure(2)
% plot(lambdaS, plS,'*', lambdaQ, plQ, 'x', lambdaO, plO, '*', M', pl, 'o')
% ylim([0.7 1.2])
% xlim([min(M)-0.01 237.15])
% legend('Simulation S-Branch', 'Simulation Q-Branch', 'Simulation O-Branch', 'Mearurement');
% 
% %save('PISpectrum','Pop','lambda*');
% 




