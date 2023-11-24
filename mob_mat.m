function [i,j,value] = mob_mat(D,nx,dx,dt)
% This function calculates the local mobility matrix containing the
% diffusive-dispersive fluxes
% - M. Muniruzzaman / August 18, 2019

D(D==0)=1e-100; % add a small value to avoid zero
%(we cannot extract the indices and values if all entries are zero!)

% Preparation of sparse matrix entires [nx,nx]:
P1 = -(D(1:end-1)+D(2:end)); 
P  = [-D(1);P1;-D(end)];     % Entires of the main diagonal
P(1) = -2*D(1);

Q  = [D;0];                  % Entires of the first sub-diagonal
R  = [0;D];                  % Entires of first super-diagonal

% Formation of the sparse mobility matrix
Adt = spdiags(P*dt/dx,0,nx,nx) ...     
    + spdiags(Q*dt/dx,-1,nx,nx) ...  
    + spdiags(R*dt/dx,1,nx,nx);   

[i,j,value]=find(Adt);
end