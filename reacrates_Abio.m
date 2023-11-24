function Mreac = reacrates_Abio(c,k_Abio_TCE,k_Abio_C2H2,k_Abio_C2H4,nx,poros_vec,dt)
% This function calculates the abiotic reaction terms for different compounds
% - M. Muniruzzaman / August 18, 2019

TCE  = c(:,1);
C2H2 = c(:,2);
C2H4 = c(:,3);
C2H6 = c(:,4);

r_TCE      = -k_Abio_TCE.*poros_vec*dt;     
    
r_C2H2_off =  k_Abio_TCE.*poros_vec*dt; 
r_C2H2     = -k_Abio_C2H2.*poros_vec*dt;

r_C2H4_off =  k_Abio_C2H2.*poros_vec*dt;
r_C2H4     = -k_Abio_C2H4.*poros_vec*dt;

r_C2H6_off =  k_Abio_C2H4.*poros_vec*dt;
           
zeromat=spalloc(nx,nx,0);    
% Reaction matrix
Mreac = [spdiags(r_TCE(:),0,nx,nx)       , zeromat                        , zeromat                       , zeromat, zeromat;...
         spdiags(r_C2H2_off(:),0,nx,nx)  , spdiags(r_C2H2(:),0,nx,nx)     , zeromat                       , zeromat, zeromat;...
         zeromat                         , spdiags(r_C2H4_off(:),0,nx,nx) , spdiags(r_C2H4(:),0,nx,nx)    , zeromat, zeromat;...
         zeromat                         , zeromat                        , spdiags(r_C2H6_off(:),0,nx,nx), zeromat, zeromat;...
         zeromat                         , zeromat                        , zeromat                       , zeromat, zeromat];
end
     