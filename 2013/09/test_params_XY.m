
% Parameters for heat engine experiment
deg = pi/180;
omega = [240e3*2*pi; 1E6*2*pi]; % Trap frequency [axial; vertical]
t_osc = 2*pi/omega(1); % Ion oscillation period

% Laser wavelengths /nm
lambda_1 = 397E-9;
lambda_3 = 866E-9;

% Axes and Angles:   
% x = trap axis       (endcap-endcap)
% y = cavity axis     (horizontal and perp. to x)
% z = vertical axis   (up direction)

% B-field /gauss (rough estimate)
Bx=0.1;
By=0.1;
Bz=0.1;

% Laser direction
alpha = 85*deg;   % angle with trap axis in horizontal plane
beta  = 0*deg;    % angle out of horizontal plane

% laser polarization  (x,y,sigma+ or sigma-) for coherent and incoherent pump
pol1 = 'y';
pol3 = 'y';

% Laser powers and spot sizes (SI)
%P_1 = 4.5E-3;
P_1 = 1E-4;
waist1_0 = 50E-6;
%P_1 = 5E-6;
%waist1_0 = 0.5*68E-6;
waist1_x = waist1_0/sin(alpha);
waist1_y = waist1_0/cos(alpha);

%Position of focus of lasers on the trap axis (x)
x_0 = -waist1_x*(0/4);  

%P_3 = 1E-5;
P_3 = 0;
waist3_0 = 145E-6;
%P_3 = 0.0005E-6; 
%waist3_0 = 0.5*68E-6;
waist3_x = waist3_0/sin(alpha);
waist3_y = waist3_0/cos(alpha);

% Laser detunings (units of Gamma)
%Delta1_0 = -1.0;
Delta1 = -0.0;
Delta3 = +10;

% incoherent pumping 
GammaPump1 = 0;
GammaPump3 = 0.53;

% laser linewidths (in units of Gamma)
gamma_1 = 1E6*2*pi/Gamma; %linewidth of cooling laser (ground to p state)
gamma_3 = 1E6*2*pi/Gamma; %linewidth of repumping laser (excited to p state)
gamma_rel = gamma_1 + gamma_3; %combined linewidth (gnd to exc)

% --Derived variables--
Freq = @(lambda) 2*pi*c/lambda;
Intensity = @(power, waist) 2*power/(pi*waist^2);
E_field = @(intensity) sqrt(2*intensity/(c*epsilon_0));
% Rabi frequency in units of Gamma
Rabi_freq = @(Gamma_t, lambda, power, waist) sqrt(3*pi*epsilon_0*hbar*c^3*Gamma_t/(Freq(lambda)^3))*(E_field(Intensity(power, waist))/hbar)/Gamma;

% cooling laser center rabi frequency (units of Gamma)
Omega1_0 = Rabi_freq(Gamma, lambda_1, P_1, waist1_0);
s_1 = 2*Omega1_0^2;

% repumper center rabi frequency (units of Gamma)
Omega3_0 = Rabi_freq(Gamma/12, lambda_3, P_3, waist3_0);
s_3 = 170*2*Omega3_0^2;

