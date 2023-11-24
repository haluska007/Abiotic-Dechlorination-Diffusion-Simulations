function [D_Interf]=calc_Dmat(Dw,nx,dx,ncomp,poros_w)
% This function calculates mass-transfer coefficients at the
% cell-interfaces
% - M. Muniruzzaman / August 18, 2019

% initialization or vectors and matrices
b        = zeros(nx,ncomp);
b_Interf = zeros(nx-1,ncomp);
D_Interf = zeros(nx-1,ncomp);

% set conc= 0 to very small value to avoid NAN in denom
% c(c<=0)=1e-20;

for i=1:ncomp
% coefficients including invidual cell-properties
b(:,i)    = (2/dx).*( poros_w.*Dw(i) ); 
end

% coefficients at the cell interfaces
b_Interf(:,:) = b(1:end-1,:).*b(2:end,:)./(b(1:end-1,:)+b(2:end,:));

% weighted diffusion/dispersion coefficients at cell interfaces 
D_Interf(:,:) = b_Interf;
end