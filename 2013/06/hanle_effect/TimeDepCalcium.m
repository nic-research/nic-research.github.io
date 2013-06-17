function [ occupation_E1, occupation_E2, occupation_E3] = ...
    TimeDepCalcium( alphaUV, betaUV, alphaIR, betaIR, Bx, By, Bz, DeltaB, Omega1, Omega3, Delta1, Delta3, ... 
            pol1, pol3, GammaPump1, GammaPump3, omega, tlist)
        
% full Zeeman level structure with arbitrary magnetic field included
% Polarization of incoming beams converted to frame of reference with
% z-axis in direction of B
% Omega1, Omega2, Delta1, GammaPump1 are in units of Gamma
% Omega1 and Omega2 are the Rabi frequencies of the two ions
% taulist is a list of delay times in units of 1/Gamma
% d_ions is the ion distance in �m
% alpha, beta are the angles of the laser with the trap axis and out of the horizontal plane
% theta_out is the angle of the laser beam and the direction of observation 
% with the trap axis 

name = 'PulsedPump';

% Calcium system:
Gamma1  = 1;
Gamma3  = 0.076;

tpulse = 1000;
wpulse = 300;
tmax = tlist(length(tlist));
tr = tlist(2);

J1 = 1/2;    %  ground state      S 1/2    1
N1 = 2*J1+1;    
J2 = 1/2;    %  excited state     P 1/2    2
N2 = 2*J2+1;    
J3 = 3/2;    %  metastable state  D 3/2    3
N3 = 2*J3+1;    

M1 = -J1:J1;
M2 = -J2:J2;
M3 = -J3:J3;

E1 =            (1:N1);
E2 = N1       + (1:N2);
E3 = N1+N2    + (1:N3);

Nat = N1+N2+N3;
idat = identity(Nat);

% Zeeman splittings for the three levels per Gauss

w1 = zeemanFS(1/2,1/2,0,1E-4,140E6); 
w2 = zeemanFS(1/2,1/2,1,1E-4,140E6); 
w3 = zeemanFS(3/2,1/2,2,1E-4,140E6); 

% transition ( |1>  |2>)
A1m = sparse(Nat,Nat); 
A10 = sparse(Nat,Nat); 
A1p = sparse(Nat,Nat); 
[am,a0,ap] = murelj(J1,J2);
A1m(E1,E2)=am;
A10(E1,E2)=a0;
A1p(E1,E2)=ap;

A1m = qo(A1m);
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

PE1 = tensor(basis(Nat,E1));
PE2 = tensor(basis(Nat,E2));
PE3 = tensor(basis(Nat,E3));

ProjE1=PE1*PE1';
ProjE2=PE2*PE2';
ProjE3=PE3*PE3';

% Zeeman splitting (per Gauss)
HB = tensor(qo(diag([ w1*M1 w2*M2 w3*M3 ])));
HBd = HB';
LB  = -1i*(spre(HB)  - spost(HB));
LBd  = -1i*(spre(HBd)  - spost(HBd));

% New assignement of angles   
kUVx =  cos(alphaUV) * cos(betaUV);           % x = trap axis       (shuttle direction)
kUVy = -sin(alphaUV) * cos(betaUV);           % y = cavity axis     (output to input mirror)
kUVz =                 sin(betaUV);           % z = vertical axis   (up direction)
kUVe = [ kUVx ; kUVy ; kUVz ];

kIRx =  cos(alphaIR) * cos(betaIR);           % x = trap axis       (shuttle direction)
kIRy = -sin(alphaIR) * cos(betaIR);           % y = cavity axis     (output to input mirror)
kIRz =                 sin(betaIR);           % z = vertical axis   (up direction)
kIRe = [ kIRx ; kIRy ; kIRz ];

B  = [ Bx ; By ; Bz ];

if norm(B)==0, 
    B = [0 0 0.0001]'; 
end
    
% Calculate polarization in frame of reference with z-axis in B-direction
[A1Bp,A1B0,A1Bm] = RotAtomPol([A1p,A10,A1m],B,kUVe);
[A3Bp,A3B0,A3Bm] = RotAtomPol([A3p,A30,A3m],B,kIRe);

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

% Loss terms

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

% Atomic detuning terms 

H0 = Delta1 * sum(ProjE1) +  Delta3 * sum(ProjE3);
L0  = -1i*(spre(H0)  - spost(H0));

%initial state

rho_atom=sparse(Nat,Nat);
rho_atom(1,1)=0.5;
rho_atom(2,2)=0.5;
rho0=qo(rho_atom);

% Build function series for Liouville operator

L =  L0 + (norm(B) + DeltaB * fn('pulse',1i*omega,0,tr,tmax,tr)) * LB  + (norm(B) + DeltaB * fn('pulse',-1i*omega,0,tr,tmax,tr)) * LBd  + ...
     Omega1 * L1 + Omega1' * L1d + ...
     Omega3 * L3 + Omega3' * L3d + ...
     Gamma1  * sum(LC1)  + ...
     Gamma3  * sum(LC3)  + ...
     GammaPump3 * sum(LPu3) + ...
     GammaPump1 * sum(LPu1);
  
save 'lioville' L
       
           
% find linked subspace
[Ls, rho0s, is] = findsubspace(L,rho0);

% full subspace
%Ls = L;
%rho0s = qO(rho0(:));
%shL = shape(L);
%is=(1:shL(1))';

disp(['full L: ' int2str(shape(L)) ';  reduced Ls: ' int2str(shape(Ls))]);

% Solve the differential equation

ode2file('setup.dat',L,rho0,tlist);      
odesolve('setup.dat',[ name '.dat' ]);      
fid = fopen([ name '.dat'],'rb');
clear rhos;
rhos = qoread(fid, dims(rho0),size(tlist));
fclose(fid);
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for k=1:N1, occupation_E1(k,:)=abs(expect(ProjE1{k},rhos)); end;   
for k=1:N2, occupation_E2(k,:)=abs(expect(ProjE2{k},rhos)); end;   
for k=1:N3, occupation_E3(k,:)=abs(expect(ProjE3{k},rhos)); end; 





  

    
    
 
     
    
    

