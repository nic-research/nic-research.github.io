% Function for simulating Hanle effect
% N.Seymour-Smith 31/1/13                      

% full Zeeman level structure with arbitrary magnetic field included
% Polarization of incoming beams converted to frame of reference with
% z-axis in direction of B
% Omega1, Omega2, Delta1, GammaPump1 are in units of Gamma
% Omega1 and Omega2 are the Rabi frequencies of the two ions
% alpha, beta are the angles of the laser with the trap axis and out of the horizontal plane
% theta_out is the angle of the laser beam and the direction of observation 
% with the trap axis 

% wfl 05-04-2003   

tic
% constants
hbar = 6.6e-34/(2*pi); % Reduced Planck constant
lambda = 397e-9;    % Laser wavelength
k = 2*pi/lambda;    % Laser wavenumber
m = 40*1.66053886e-27; % Mass of Ca
Gamma = 22*1e6*2*pi;     % Decay rate
omega = 240e3*2*pi; % Trap frequency
t_osc = 2*pi/omega; % Ion oscillation period
k_b = 1.381E-23;    % Boltzmann constant
epsilon_0 = 8.85e-12;

% B-field /gauss
Bx=0;
By=0;
%Bz=1;

% laser direction
deg = pi/180;
alpha = 0*deg;   % angle with trap axis in horizontal plane
beta  = 0*deg;    % angle out of horizontal plane

% laser polarization  (x,y,sigma+ or sigma-) for coherent and incoherent pump
pol1 = 'x';
pol3 = 'y';

% laser linewidths (in units of Gamma)
gamma_1 = 0E6*2*pi/Gamma; %linewidth of cooling laser (ground to p state)
gamma_3 = 0E6*2*pi/Gamma; %linewidth of repumping laser (excited to p state)
gamma_rel = gamma_1 + gamma_3; %combined linewidth (gnd to exc)

% cooling laser power and spot size (SI)
P_1 = 1E-5;
waist_1 = 200E-6;

% cooling laser rabi frequency (units of Gamma)
lambda_1 = 397E-9;
w_1 = 2*pi*c/lambda_1;
I_1 = 2*P_1/(pi*waist_1^2);
Ef_1 = sqrt(2*I_1/(c*epsilon_0));
Omega1 = sqrt(3*pi*epsilon_0*hbar*c^3*Gamma/(w_1^3))*(Ef_1/hbar)/Gamma;

% repumper power and spot size (SI)
P_3 = 1E-5;
waist_3 = 200E-6;

% repumper rabi frequency (units of Gamma)
lambda_3 = 866E-9;
w_3 = 2*pi*c/lambda_3;
I_3 = 2*P_3/(pi*waist_3^2);
Ef_3 = sqrt(2*I_3/(c*epsilon_0));
Omega3 = sqrt(3*pi*epsilon_0*hbar*c^3*Gamma*(1/12)/(w_3^3))*(Ef_3/hbar)/Gamma;

% saturation parameters
% s_1 = 10;
% s_3 = 144;
% Omega1 = sqrt(s_1/2);
% Omega3 = sqrt(s_3/2)*(1/12);

% Omega1 = 1;
% Omega3 = 1;
Delta1 = -0.75;
Delta3 = 0;

% incoherent pumping 
GammaPump1 = 0;
GammaPump3 = 0;

% B-field
dB = 0.1;  % Resolution
Bz = 0.0:dB:1;        
N_data = length(Bz); 
count = zeros(1,N_data); % array for fluorescence on cooling transition

% Timing variables and arrays
Pre_cool = 1;
Expt_sim = 1;
if Expt_sim == 1
    Pre_cool = 1;
    Eff = 0.001;                % Detector efficiency
    R = 25000;                  % probe/cool rate
    t_cycle = 1/R;              % time per cool/probe cycle
    t_probe = 10e-6;            % probe time per cycle
    t_cool = t_cycle - t_probe; % cool time
    t_point = 100e-3;           % time per data point (real-time)
    t_point = Eff*t_probe*t_point/t_cycle; % total simulation time per point
                                        % NB multiplying by efficiency is
                                        % not quite as nasty as it sounds,
                                        % as long as the result is at least
                                        % one real-time probing cycle long.
                                        % Computation time-scales are just
                                        % not practical if I don't do this
else
    t_cool = 100e-6;
    t_probe = 100e-6;
    N_probe = 1;
    t_point = N_probe*t_probe;
    R = 0;
end

if Pre_cool ~= 1
    t_cool = 0;
    Tmp_cool = 0;
end

dt = 1e-9;                  % simulation resolution initialisation
T = N_data*t_point;         % total simulation time
t = zeros(1,1e7);
dt_h = zeros(1,1e7);
T_lim = t_osc/10; %Time after last emission to search for a new...
                  %emission before recalculating posn/vel.
M = floor(T_lim/dt);         %Iteration limit based on above
Th_2 = 1/100;          % threshold probability of 2-photon events

%Measurement parameters
dt_m = 1e-6;    %PMT measurement bins
t_m = 0:dt_m:T; %fluorescence measurement bins
PMT = zeros(1,length(t_m)); %Pre-allocation for fluorescence rate measurements
Fl_c = zeros(1, N_data); %Calculated fluorescence rate
Fl_s = zeros(1, N_data); %Simulated fluorescence rate
Tmp_avg = zeros(1, N_data); % Average temperature during probe period

% Ion position, velocity and count rates
len = length(t);
% pos = zeros(1,len);      %Ion position array
% vel = zeros(1,len);      %Ion velocity array
% Fl_r = zeros(1,len);     %Ion fluorescence array (realtime)
% Fl_s = zeros(1, N_loop); %Ion fluorescence array (spectrum)
                                 
pos = 0;
vel = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Definitions of constants and operators (don't change)     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/home/nic/q_opt_toolbox')

name = 'PulsedPump';

%lambda = 0.397;     % wavelength in �m

% Calcium system:
Gamma1  = 1;
Gamma3  = 0.076;

J1 = 1/2;    %  ground state      S 1/2    1
N1 = 2*J1+1; %  number of zeeman sub-levels   
J2 = 1/2;    %  excited state     P 1/2    2
N2 = 2*J2+1;    
J3 = 3/2;    %  metastable state  D 3/2    3
N3 = 2*J3+1;    

M1 = -J1:J1; %  array of sub-levels
M2 = -J2:J2;
M3 = -J3:J3;

E1 =            (1:N1); %  numbered energy levels, split into arrays for 
E2 = N1       + (1:N2); %  each manifold
E3 = N1+N2    + (1:N3);

Nat = N1+N2+N3; % total energy levels
idat = identity(Nat); % Identity matrix, N = Nat dimensions

% Zeeman splittings for the three levels per Gauss
% function w = zeemanFS(J,S,L,B,Gamma_phys)

w1 = zeemanFS(1/2,1/2,0,1E-4,Gamma); 
w2 = zeemanFS(1/2,1/2,1,1E-4,Gamma); 
w3 = zeemanFS(3/2,1/2,2,1E-4,Gamma); 

% transition ( |1>  |2>)
A1m = sparse(Nat,Nat); % sigma-minus transition operator (not initialised)
A10 = sparse(Nat,Nat); % pi transition operator
A1p = sparse(Nat,Nat); % sigma-plus transition operator
[am,a0,ap] = murelj(J1,J2); % submatrices of Clebsch-Gordan coefficients
                            % for sigma-minus, pi and sigma-plus
A1m(E1,E2)=am; % Full lowering transition operators
A10(E1,E2)=a0;
A1p(E1,E2)=ap;

A1m = qo(A1m); % Convert to quantum optics toolbox `quantum object'
A10 = qo(A10);
A1p = qo(A1p);

% transition ( |3>  |2>)
A3m = sparse(Nat,Nat); 
A30 = sparse(Nat,Nat); 
A3p = sparse(Nat,Nat); 
[am,a0,ap] = murelj(J3,J2);
A3m(E3,E2)=am;
A30(E3,E2)=a0;
A3p(E3,E2)=ap;

A3m = qo(A3m);
A30 = qo(A30);
A3p = qo(A3p);

% single atom projection operators
% basis(n,index) gives the mth ket vector in n-dim space, with m = index
% E1, E2, E3 are vectors, representing the number of the levels in each
% manifold. Therefore "basis(Nat,E1)" returns a vector of the kets
% associated with those numbers/energy levels.

PE1 = basis(Nat,E1);
PE2 = basis(Nat,E2);
PE3 = basis(Nat,E3);

ProjE1=PE1*PE1'; % This gives us a vector of operators, each element
ProjE2=PE2*PE2'; % associated with a different sublevel
ProjE3=PE3*PE3';

% Zeeman splitting (per Gauss):
% "diag([ w1*M1 w2*M2 w3*M3 ])" is a diagonal matrix of the Zeeman
% splittings with respect to their associated energy levels
HB = qo(diag([ w1*M1 w2*M2 w3*M3 ]));

% Convert above to superoperator (for easier computation)
LB  = -1i*(spre(HB)  - spost(HB));

% New assignement of angles   
kx =  cos(alpha) * cos(beta);           % x = trap axis       (endcap-endcap)
ky = -sin(alpha) * cos(beta);           % y = cavity axis     (horizontal and perp. to above)
kz =               sin(beta);           % z = vertical axis   (up direction)

ke = [ kx ; ky ; kz ];
B  = [ Bx ; By ; Bz ];

if norm(B)==0, 
    B = [0.00001 0 0]'; 
end
    
% Calculate polarization in frame of reference with z-axis in B-direction
[A1Bp,A1B0,A1Bm] = RotAtomPol([A1p,A10,A1m],B,ke);
[A3Bp,A3B0,A3Bm] = RotAtomPol([A3p,A30,A3m],B,ke);

% Polarization in frame of laser (can be sigma+, sigma- or linear combination)

if pol1=='y',
    % y-polarization
    H1 = -1i*(A1Bp-A1Bm)/sqrt(2) ;  
    % incoherent pump in B-frame
    Pu1  = -1i * [  A1Bp, 0 * A1B0,  -A1Bm ] / sqrt(2);
elseif pol1=='x',
    % x-polarization
    H1 = (A1Bp+A1Bm)/sqrt(2) ;  
    Pu1  = [  A1Bp, 0 * A1B0,  -A1Bm ] / sqrt(2);
elseif pol1=='sigma+',
    % sigma+ polarization
    H1 = A1Bp;  
    Pu1  = [  A1Bp, 0 * A1B0,  0 * A1Bm ] ;
elseif pol1=='sigma-',
    % sigma- polarization
    H1 = A1Bm;  
    Pu1  = [  0 * A1Bp, 0 * A1B0,  A1Bm ] ;
else
    warndlg([pol1 ' is invalid polarization for laser 1'],'Solver warning'); 
    return
end


%disp('input polarization in atomic frame');
a1p = -H1(1,4)*sqrt(3/2);      %  A1p contribution
a1m = H1(2,3)*sqrt(3/2);      %  A1m contribution
a10 = -H1(1,3)*sqrt(3);       %  A10 contribution
if ~(a10==H1(2,4)*sqrt(3)), disp('error in polarization analysis'); end;       %  A10 contribution

if pol3=='y',
    % y-polarization
    H3 = -1i*(A3Bp-A3Bm)/sqrt(2) ;  
    % incoherent pump in B-frame
    Pu3  = -1i * [  A3Bp, 0 * A3B0,  -A3Bm ] / sqrt(2);
elseif pol3=='x',
    % x-polarization
    H3 = (A3Bp+A3Bm)/sqrt(2) ;  
    Pu3  = [  A3Bp, 0 * A3B0,  A3Bm ] / sqrt(2);
elseif pol3=='sigma+',
    % sigma+ polarization 
    H3 = A3Bp;  
    Pu3  = [   A3Bp, 0 * A3B0,  0 * A3Bm ];
elseif pol3=='sigma-',
    % sigma- polarization
    H3 = A3Bm;  
    Pu3  = [ 0* A3Bp, 0 * A3B0,  A3Bm];
else
    warndlg([pol3 ' is invalid polarization for laser 3'],'Solver warning'); 
    return;
end


H1d = H1';  H3d = H3';

L1  = -1i*(spre(H1)  - spost(H1));
L1d = -1i*(spre(H1d) - spost(H1d));
L3  = -1i*(spre(H3)  - spost(H3));
L3d = -1i*(spre(H3d) - spost(H3d));

% Loss terms (normal spontaneous decay)

C1   = [A1m,A10,A1p];   
C3   = [A3m,A30,A3p]; 

C1dC1 = C1'*C1;
C3dC3 = C3'*C3;

LC1 = spre(C1)*spost(C1')-0.5*spre(C1dC1)-0.5*spost(C1dC1);
LC3 = spre(C3)*spost(C3')-0.5*spre(C3dC3)-0.5*spost(C3dC3);

% incoherent Pump terms (in B-frame)

Pu3Pu3d = Pu3*Pu3';
LPu3 = spre(Pu3')*spost(Pu3)-0.5*spre(Pu3Pu3d)-0.5*spost(Pu3Pu3d);

Pu1Pu1d = Pu1*Pu1';
LPu1 = spre(Pu1')*spost(Pu1)-0.5*spre(Pu1Pu1d)-0.5*spost(Pu1Pu1d);

% Laser decoherence terms c.f. Pritchard thesis p43
% List of dephasing terms for transitions to one of the ground states:
Dgg = [0 0];
Dgp = -gamma_1*[1 1];
Dge = -gamma_rel*[1 1 1 1];

% Dephasing for transitions to one of the p states
Dpg = -gamma_1*[1 1];
Dpp = [0 0];
Dpe = -gamma_3*[1 1 1 1];

% Dephasing for transitions to one of the excited states (d1/2)
Deg = -gamma_rel*[1 1];
Dep = -gamma_3*[1 1];
Dee = [0 0 0 0];

% Array of dephasing terms for transitions to all ground states
DG_ = cat(2, Dgg, Dgp, Dge);
DG = cat(2, DG_, DG_);

% Dephasing terms for transition to all p states
DP_ = cat(2, Dpg, Dpp, Dpe);
DP = cat(2, DP_, DP_);

% Dephasing terms for transitions to all d states
DE_ = cat(2, Deg, Dep, Dee);
DE = cat(2, DE_, DE_, DE_, DE_);

% Specify the dimensions of the pre/post-multiplication operator ({8 8} is
% the dimension of the Hilbert space)
C = {{8 8} {8 8}};

% Create a matrix with the dephasing terms in the format used by the
% toolbox (see the manual), and convert to quantum object of the dimensions
% of a pre/post-multiplication superoperator
LD = qo(diag(cat(2, DG, DP, DE)'), C); 
                                                                    
% initial condition for densiy matrix
rho_atom=sparse(Nat,Nat);
rho_atom(1,1)=0.5; 
rho_atom(2,2)=0.5; %Equally populated S_1 states
rho_atom=qo(rho_atom);
rho0     = rho_atom; 

% % Atomic detuning terms 
% 
% H0 = Delta1 * sum(ProjE1) +  Delta3 * sum(ProjE3);
% L0  = -1i*(spre(H0)  - spost(H0));

% Build function series for Liouville operator
% (Deleted a term with E*L4 because neither were defined)
L =  norm(B) * LB + ...
     Omega1 * L1 + Omega1' * L1d + ...
     Omega3 * L3 + Omega3' * L3d + ...
     Gamma1  * sum(LC1)  + ...
     Gamma3  * sum(LC3)  + ...
     GammaPump3 * sum(LPu3) + ...
     GammaPump1 * sum(LPu1) + ...
     LD;

E_cool = 0;

%%%%% Initial Doppler temperature limit simulation %%%%%
if Pre_cool == 1
    N = 0; %Sum of calculations 
    t(1) = dt;
    j = 1;

    while t(j) < t_cool
        Delta1 = Delta1c+k*vel/Gamma; % spectrum detuning plus Doppler effect
        H0 = Delta1 * sum(ProjE1) + Delta3 * sum(ProjE3); % Atomic detuning terms
        L0 = -1i*(spre(H0) - spost(H0));
        L_s = L + L0;
        
        % Solve the differential equation
        rhos=steady(qo(L_s));
        % Theoretical fluorescence rate
        gamma = Gamma*real(sum(expect(C1dC1, rhos)));
        
        % Timestep for photon checks, based on limiting 2-photon events
        dt = 2*Th_2/gamma;
        M = floor(T_lim/dt); % Iteration limit based on ion oscillation
        % Probability of fluorescence in time dt
        P_kick = gamma*dt;
        
        % Assume a photon will result, reset if not (i == M)
        photon = 1;
        i = 0;
        U = 1;
        while U > P_kick
            if i < M
                U = rand(1);
            else
                U = 0;
                photon = 0;
            end
            i = i+1;
        end
        if photon == 1
            t_step = i*dt;
            angle = pi*rand(1);   %Random angle of spontaneous emission
            v_k = (hbar*k/m)*(cos(angle)-1); %velocity kick from emission (cos(angle)) and absorption (-1)
        else
            t_step = T_lim;
            v_k = 0;
        end
        
        N = N+1;
        E = (m/2)*(vel*vel + omega*omega*pos*pos); %Total energy of ion in harmonic potential
        E_cool = (E_cool*(N-1) + E)/N;               %Average energy for temperature calc.
        if E ~= 0 %The following needs to be inside this condition to avoid a divide by zero on initialisation
            x_a = sqrt(2*E/(m*omega*omega));           %Maximum amplitude of oscillation
            phi = asin(pos/x_a); %phase of oscillation N.B. not a unique solution...
            if vel < 0           %...
                phi = pi - phi;   %...determines unique phase from direction of motion
            end
            pos = x_a*sin(phi+omega*t_step);           %New position determined from current phase and time til next event
            vel = x_a*omega*cos(phi+omega*t_step);     %New velocity determined from phase and time til next event
        end
        vel = vel + v_k;     %Add velocity kick from photon event (if there was one) N.B. This is applied _after_ working out new state [x, p]
        if photon == 1
            t(j+1) = t(j) + t_step; %Time of next event
            j = j+1;                %Move on to the next event
        else
            t(j) = t(j) + t_step;
        end
    end
    
    Tmp_cool = E_cool/k_b; % Temperature after cooling period
end

%%%%%% Experiment simulation %%%%%%
E_probe = 0;
photon = 0;
N = 0; %Sum of calculations 
t(1) = dt;
j = 1;
l_last = 1;
p_last = 0;

w = waitbar(0, 'Simulating spectrum');
while t(j) < T
    
    % Track progress along spectrum
    l = ceil(t(j)*N_data/T); 
    if l ~= l_last
        N = 0; % Clear the averaging counter
    end
    l_last = l;
    
    % Simulate a cooling cycle before each probe period:
    % Check for the start of a new probe cycle:
    p = floor(t(j)/t_probe); 
    if (p ~= p_last) && (Pre_cool == 1)
        
        % Randomly distribute thermal energy over potential and kinetic
        E_pot = rand(1)*E_cool;
        pos = sqrt(2*E_pot/(m*omega*omega));
        
        E_kin = E_cool - E_pot;
        vel = sqrt(2*E_kin/m);
    end
    p_last = p;
    
    %Calculate Liouvillian for frequency including Doppler effect
    Delta1 = f(l)+k*vel/Gamma; % spectrum detuning plus Doppler effect
    H0 = Delta1 * sum(ProjE1) + Delta3 * sum(ProjE3); % Atomic detuning terms
    L0 = -1i*(spre(H0) - spost(H0));
    L_s = L + L0;
    % Solve the differential equation
    rhos=steady(qo(L_s));
    % Theoretical fluorescence rate
    gamma = Gamma*real(sum(expect(C1dC1, rhos)));
    
    % Running average of the calculated fluorescence at this frequency
    N = N+1;
    Fl_c(l) = (Fl_c(l)*(N-1) + gamma*t_point)/N;
    % Sum the photons at this frequency
    if photon == 1
        Fl_s(l) = Fl_s(l) + 1;
    end

    % Timestep for photon checks, based on limiting 2-photon events
    dt = 2*Th_2/gamma;
    M = ceil(T_lim/dt); % Iteration limit based on ion oscillation
    % Probability of fluorescence in time dt
    P_kick = gamma*dt;
    
    % Assume a photon will result, reset if not (i == M)
    photon = 1;
    i = 0;
    U = 1;
    while U > P_kick
        U = rand(1);
        i = i+1;
        if i == M
            U = 0;
            photon = 0;
        end
    end
    
    if photon == 1
        dt_h(j) = dt;
        t_step = i*dt;
        angle = pi*rand(1);   %Random angle of spontaneous emission
        %angle = pi*(P_kick-U)/P_kick;
        v_k = (hbar*k/m)*(cos(angle)-1); %velocity kick from emission (cos(angle)) and absorption (-1)
    else
        t_step = T_lim;
        v_k = 0;
    end
    
    E = (m/2)*(vel*vel + omega*omega*pos*pos); %Total energy of ion in harmonic potential
    E_probe = (E_probe*(N-1) + E)/N;               %Average energy for temperature calc.
    Tmp_probe = E_probe/k_b;
    Tmp_avg(l) = (Tmp_avg(l)*(N-1) + Tmp_probe)/N; 
    if E ~= 0 %The following needs to be inside this condition to avoid a divide by zero on initialisation
        x_a = sqrt(2*E/(m*omega*omega));           %Maximum amplitude of oscillation
        phi = asin(pos/x_a); %phase of oscillation N.B. not a unique solution...
        if vel < 0           %...
            phi = pi - phi;   %...determines unique phase from direction of motion
        end
        pos = x_a*sin(phi+omega*t_step);           %New position determined from current phase and time til next event
        vel = x_a*omega*cos(phi+omega*t_step);     %New velocity determined from phase and time til next event
    end
    vel = vel + v_k;     %Add velocity kick from photon event (if there was one) N.B. This is applied _after_ working out new state [x, p]
    if photon == 1
        t(j+1) = t(j) + t_step; %Time of next event
        j = j+1;                %Move on to the next event
    else
        t(j) = t(j) + t_step;
    end

    %Timing and progress bar stuff:
    sim_time_left = T - t(j);
    precision = 3;

    waitbar(t(j)/T,w,...
        {['Tmp_{cooled} = ',num2str(Tmp_cool*1e6,3),'uK, ', ...
        'Tmp_{probe} = ',num2str(Tmp_avg(l)*1e6,3),'uK'], ...
        ['Experiment-time left to simulate =', ...
        num2str(sim_time_left*1e3,precision),'ms ']})

end
close(w)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t(1:j-1); %truncate time vector at end of simulation
dt_h = dt_h(2:j-1);
dt_min = min(dt_h);
dt_max = max(dt_h);
dt_mean = mean(dt_h);

% j = 1;
% for l = 1:length(t_m)
%     if j < length(t) - 1
%         while t(j) < t_m(l)
%             PMT(l) = PMT(l) + 1;
%             j = j+1;
%         end
%     end
% end
% PMT = PMT(1:l);
% t_m = t_m(1:l);

%Construct parameter string for plots
s1 = ['[B_x, B_y, B_z] = [', num2str(Bx),', ', num2str(By), ...
    ', ', num2str(Bz),'] gauss, '];

s2 = ['pol1 = ', pol1, ', pol3 = ', pol3',', '];

s3 = ['P_1 = ', num2str(P_1*1E6), ' uW, ',...
    '\Omega_1 = ', num2str(Omega1,2), ' \Gamma, ',...
    's_1 = ', num2str(s_1,2), ', '...
    'P_3 = ', num2str(P_3*1E6), ' uW, ',...
    '\Omega_3 = ', num2str(Omega3,2), ' \Gamma, '...
    's_3 = ', num2str(s_3,2),', '];

s4 = ['\Delta_1 = (', num2str(f(1),2),'...', num2str(f(end),2), ...
    ') \Gamma, \Delta_c = ', num2str(Delta1c), '\Gamma, \Delta_3 = ',...
    num2str(Delta3), '\Gamma, \deltaFreq = ', num2str(dB*Gamma*1E-6/(2*pi)), ...
    ' MHz, '];

s5 = ['R_{probe} = ', num2str(R*1e-3), ' kHz, t_{cool} = ', ...
    num2str(t_cool*1e6), ' us, t_{probe} = ', num2str(t_probe*1e6), ' us, ',...
    't_{point} = ', num2str(t_point*1e3), ' ms, T_{sim} = ', ...
    num2str(T*1e3), ' ms, '];

s6 = ['Temp_{cooled} = ', num2str(Tmp_cool*1e6,3), ' uK, ',...
    'Temp_{max} = ', num2str(max(Tmp_avg)*1e6,3), ' uK, ',...
    'Temp_{min} = ', num2str(min(Tmp_avg)*1e6,3), ' uK, ',...
    'Temp_{avg} = ', num2str(mean(Tmp_avg)*1e6,3), ' uK, '];

s7 = ['Timesteps: \deltat_{max} =', num2str(dt_max*1e9,3), ' ns, ',...
    '\deltat_{min} =', num2str(dt_min*1e9,3),' ns, ',...
    '\deltat_{mean} =', num2str(dt_mean*1e9,3),' ns'];

s8 = ['Laser linewidths: \gamma_{1} = ', ...
    num2str(gamma_1*Gamma*1E-6/(2*pi), 2), ' MHz, \gamma_{3} = ', ...
    num2str(gamma_3*Gamma*1E-6/(2*pi), 2), ' MHz, \gamma_{rel} = ', ...
    num2str(gamma_rel*Gamma*1E-6/(2*pi), 2), ' MHz'];

% figure
% plot(t_m, PMT);
% title({'PMT fluorescence', [s1,s2], s3, s4, s5, s6});
figure
plot(f, Fl_s, f, Fl_c, 'r--');
title({'Simulation parameters', [s1,s2], s3, s4, s8, s5, s6,s7,''});
xlabel('\Delta_1 /\Gamma')
ylabel('N_{\gamma}')
% figure
% plot(f, Tmp_avg*1e6);
% title({'Temperature vs Frequency', [s1,s2], s3, s4, s5, s6,s7,''});
% xlabel('\Delta_1 /\Gamma')
% ylabel('T /uK')