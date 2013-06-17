function [ABp,AB0,ABm] = RotAtomPol(A,B,k)

% [ABp,AB0,ABm] = RotAtomPol([Ap,A0,Am],B,k) calculates the atomic operators for a beam with z-direction 
% given by k from the original operators [Ap,A0,Am] specified for a qunatization axis given by B

% wfl 29-04-2003  

Ap =A{1};
A0 =A{2};
Am =A{3};

Be = B / norm(B);
ke = k / norm(k);

v     = cross(Be,ke);
theta = acos(dot(Be,ke));
if norm(v)>0, v = v/norm(v); end;

%   Rodrigues formula
%   cross([v v v],eye(3)) = matrix that computes the cross product of v with any other vector
%   R  = (cos(theta)*eye(3,3)+(1-cos(theta))*v*v'-sin(theta)*cross([v v v],eye(3)))';
%   R  = eye(3,3) + cross([v v v],eye(3)) * sin(theta) + cross([v v v],eye(3))^2 * (1-cos(theta));
%   R  = expm(cross([v v v],eye(3))*theta); 

C = cross([v v v],eye(3));
R = eye(3,3) + C * sin(theta) + C^2 * (1-cos(theta));
 
%   Tp0m_xyz' [ax; ay; az]  -->  [ap; a0; am]
Tp0m_xyz = [ 1 0 1 ; i 0 -i ; 0 sqrt(2) 0 ]  / sqrt(2) ;
    
%   Transform from coordinate system with z in beam direction (polarization x,y) to 
%   absolute laboratory frame with z up 
%ez = [0; 0; 1];
%vk = cross(ez,ke);
%thetak = acos(dot(ez,ke));
%if norm(vk)>0, vk = vk/norm(vk); end;
%    
%Ck = cross([vk vk vk],eye(3));
%Rk  =  eye(3,3) + Ck * sin(thetak) + Ck^2 * (1-cos(thetak));
%    
%ex = [1; 0; 0];
%vi = cross(ex,ke);
%thetai = acos(dot(ex,ke));
%if norm(vi)>0, vi = vi/norm(vi); end;
%    
%Ci = cross([vi vi vi],eye(3));
%Ri =  eye(3,3) + Ci * sin(thetai) + Ci^2 * (1-cos(thetai));

% Transform to polarization unity vectors 
epol_z = ke;
ez = [0; 0; 1];
epol_y = cross(ez,ke);
if norm(epol_y)>0, epol_y = epol_y/norm(epol_y); end;
epol_x = cross(epol_y,epol_z);
if norm(epol_x)>0, epol_x = epol_x/norm(epol_x); end;
Rpol = [ epol_x epol_y epol_z ];

% Rpol = Rk*Ri;
% M  = Tp0m_xyz'*Ri'*Rk'*R*Rk*Ri*Tp0m_xyz;

M  = Tp0m_xyz'*Rpol'*R*Rpol*Tp0m_xyz;


% Calculate polarization in frame of reference with z-axis in B-direction

ABp = Ap * M(1,1) + A0 * M(1,2) + Am * M(1,3) ;   
AB0 = Ap * M(2,1) + A0 * M(2,2) + Am * M(2,3) ;   
ABm = Ap * M(3,1) + A0 * M(3,2) + Am * M(3,3) ;   