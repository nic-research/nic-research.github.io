function [rhss, ics, is, isd, Cs] = findsubspace(rhs,ic,collapses)

% finds subspace linked to initial state "ic" by function series "rhs"
% and collection of collapse operators "C1", "C2", ...
% Application:
%       [Ls, rhos, is, isd] = findsubspace(L,rho)
% projects "L" and "rho" into non-structured subspace of Hilbertspace
% "rhos" is vector of density-matrix elements in linked subspace
% "Ls" is Liouvillian in linked subspace
%       [Hs, psis, is, isd, Cs] = findsubspace(H,psi,{C1, C2, C3})
% projects "H", "psi" and collapse operators "C1", "C2", ... into non-structured 
% subspace of Hilbertspace
% "psis" is state vector in linked subspace
% "Hs" is Hamiltonian in linked subspace
% "Cs" is cell array of collapse operators corresponding to {C1, C2, ... }
% "is" is index-vector of linked subspace elements
% "isd" is index-vector of diagonal elements in rhos 

%   Revision history:
%   28-Dec-01   W. Lange    
%   19-Jan-02   W. Lange  support for collapse operators

if nargin<3, coll = 0; else coll = 1; end
        
if ~isa(ic,'qo')
   error('ic parameter should be a quantum object');
end
if ~isequal(prod(size(ic)),1)
   error('ic parameter should have one member');
end   
N = prod(shape(ic));
rhs = fseries(rhs);
if ~isequal([N,N],shape(rhs))
   error('rhs matrices are of the wrong shape');
end

% Check that rhs is a function series
if ~isa(ic,'qo')
   error('ic should be a quantum object');
end
if ~isequal(prod(size(ic)),1)
   error('ic should have one member');
end   
N = prod(shape(ic));

rhs = fseries(rhs);
if ~isequal([N,N],shape(rhs))
   error('rhs matrices are of the wrong shape');
end

% Check that collapse matrices are function series
if coll,
    ncollapses = length(collapses)
    ncoll = zeros(1,ncollapses);
    for k = 1:length(collapses)
        collapses{k} = fseries(collapses{k});
        if ~isequal([N,N],shape(collapses{k}))
   	        error('collapse matrices are of the wrong shape');
        end
        ncoll(k) = nterms(collapses{k});
    end
end

% mark nonzero elements of initial state
mask  =  spones(double(ic));
mask  =  mask(:);
nR = nterms(rhs);

% find mask for operators used in rhs function series
for k = 1:nR
   Omask{k} = spones(double(rhs{k}.ampl));
end

% find mask for operators used in collapse function series
if coll,
    n=1;
    for k = 1:length(collapses)
        C = collapses(k);
        for j=1:ncoll(k);
            Cmask{n} = spones(double(C{j}.ampl));
            n=n+1;
        end
    end
    nC = length(Cmask);  
end
    
% iterate to find matrix elements that will become nonzero
if coll,
    while 1,
        mask1 = mask;
        for k = 1:nR
            mask1 = mask1 | spones( Omask{k} * mask );
        end;
        for k = 1:nC;
            mask1 = mask1 | spones( Cmask{k} * mask );
        end
        if all(mask1 == mask), break; end
        mask = mask1;
    end
else
    while 1,
        mask1 = mask;
        for k = 1:nR
            mask1 = mask1 | spones( Omask{k} * mask );
        end;
        if all(mask1 == mask), break; end
        mask = mask1;
    end
end

% indices of nonzero elements
is = find(mask);

% find rhs operator in subspace
ampl = rhs{1}.ampl;
rhss   =  fseries(ampl(is,is),rhs{1}.type,rhs{1}.params);

if nR>1,
 for k = 2:nR
   ampl = rhs{k}.ampl;  
   rhss   = rhss + fseries(ampl(is,is),rhs{k}.type,rhs{k}.params);
 end   
end

% find collapse operators in subspace
if coll,
for k = 1:length(collapses)
    C = collapses(k);
    ampl = C{1}.ampl;
    Cs{k}   =  fseries(ampl(is,is),C{1}.type,C{1}.params);
    if ncoll(k)>1,
        for j=2:ncoll(k);
            ampl = C{j}.ampl;  
            Cs{k}= Cs{k} + fseries(ampl(is,is),C{j}.type,C{j}.params);
        end 
    end
end
end

% find initial condition in subspace
ics=ic(:);
ics=qo(ic(is));

% find index vector of diagonal elements in new density matrix
diagonal=speye(shape(ic));
diagonal=diagonal(:);
isd = find(diagonal(is));

