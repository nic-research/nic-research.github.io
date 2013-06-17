function w = zeemanFS(J,S,L,B,Gamma_phys)
%
%  w = zeeman(J,S,L,F,I,B) calculates the angular frequency shift 
%  in terms of the linewidth Gamma_phys for unit change of m 
%  for an atom in magnetic field B (in Tesla)
%
%  e.g. Ca+:   Gamma_phys = 2*pi*22.282E6 = 140E6;
%
%  application:
%  Zeeman-shift for a field of 1 Gauss = 10^-4 Tesla
%  fzeeman = zeemanFS(J,S,L,1E-4,2*pi*1E6)  in MHz 
%
%  zeemanFS(1/2,1/2,0,1E-4,2*pi*1E6) = 2.7993
%  zeemanFS(1/2,1/2,1,1E-4,2*pi*1E6) = 0.9331
%  zeemanFS(3/2,1/2,1,1E-4,2*pi*1E6) = 1.8662
%  zeemanFS(3/2,1/2,2,1E-4,2*pi*1E6) = 1.1197
%  zeemanFS(5/2,1/2,2,1E-4,2*pi*1E6) = 1.6796

%  w1 = zeemanFS(1/2,1/2,0,B,140E6); 
%  w2 = zeemanFS(1/2,1/2,1,B,140E6); 
%  w3 = zeemanFS(3/2,1/2,2,B,140E6); 
%  w4 = zeemanFS(3/2,1/2,1,B,140E6); 
%  Hzeeman = tensor(ida,qo(diag([ w1*M1 w2*M2 w3*M3 w4*M4 ]))); 

muB = 9.274078E-24;
hbar = 1.054572669125E-34;
%
E = (1+(J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1)))* muB*B;
w = E/(hbar*Gamma_phys);

