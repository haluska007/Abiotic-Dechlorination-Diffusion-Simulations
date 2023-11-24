% =========================================================================
% 1-D Diffusive-Reactive Transport of Chlorinated Hydrocarbons
% Processes:
% - 1D pore diffusion
% - Non-linear (Freundlich) sorption
% - Kinetic abiotic reactions

% Muhammad Muniruzzaman
% Geolgical Survey of Finland
% md.muniruzzaman@gtk.fi
% August 20, 2019
% =========================================================================
%%
clear all
% close all
clc
%======================== TRANSPORT COEFFICIENTS ==========================
L     = 1e-2;%0.5;      % Length of the column [m]
r     = 5.00E-2/2;      % Radius of the column [m]
A     = pi*r^2;         % Cross sectional area of the column [m^2]
poros = 0.10;           % Porosity [-]
te    = 86400*365;      % Total simulation time [s]  


% Discretizations
dx   = 1e-4;             % spatial discretization [m]
dt   = 10e3;              % time discretization [s]
x    = [0.5*dx:dx:L];    % column cell dimensions [m]
nx   = length(x);        % number of cells [-]


% Properties of aqeous/solid phases in a grid cell
Vb_c    = pi*r^2*L;      % Bulk volume of the column [m^3]
rho_b   = 2.2;           % Bulk density of the porous matrix [kg/L]
rho_s   = 2.74;          % Solid density of the porous matrix [kg/L]
RV_cell = Vb_c*1E3/nx;   % Representative volume of a cell [L]

% DERIVED COEFFICIENTS
% Vector of RV for all cells in the domain [L]
RV      = RV_cell*ones(nx,1);
RV(1)   = 5;
RV(end) = 5;

% Vector of porosities in all cells in the domain [-]
poros_w     = ones(nx,1)*poros;
poros_w(1)  = 1;
poros_w(end)= 1;

% Volume of pore space in each cell [L]
RV_w    = poros_w.*RV;

ncomp = 5; % Number of components (TCE, C2H2, C2H4, C2H6, Tracer)
name_comp = {'TCE','Acetylene','Ethene','Ethane','Tracer'};

% Sorption parameters --> [TCE, C2H2, C2H4, C2H6, Tracer]
K_Fr     = [3.2  0 0 0 0];        % Freundlich sorption coefficient [L/Kg]
n_Fr     = [0.64 0 0 0 0];        % Exponent in Freundlich isotherm, n_Fr [-]

% Molecular weight [mg/mol]
MW       = [131.4 26.04 28.05 30.07 1]*1e3 ; 

K_Fr_mol = K_Fr.*MW.^(n_Fr-1);  % K_Fr for molar concentrations [mol/L]

K_Fr_mat = ones(nx,1)*K_Fr_mol; % Matrix of K_Fr (nx x ncomp)
n_Fr_mat = ones(nx,1)*n_Fr;     % Matrix of n_Fr (nx x ncomp)

% Fixing the values at the boundaries: Inlet and outlet reservoirs
K_Fr_mat(1,:)   = 0;
K_Fr_mat(end,:) = 0;
    

% Reaction parameters
k_Abio_TCE  = 0.0032/86400;     % Abiotic TCE degradation rate [1/s]
k_Abio_C2H2 = 0.0093/86400;     % Abiotic C2H2 degradation rate [1/s]
k_Abio_C2H4 = k_Abio_C2H2*0.33; % Abiotic C2H4 degradation rate [1/s]

% Vector for Abiotic rate constants [1/s]
k_Abio_TCE_vec  = zeros(nx,1); 
k_Abio_C2H2_vec = zeros(nx,1);
k_Abio_C2H4_vec = zeros(nx,1);

% Fixing the rate cofficients at the reservoir cells
k_Abio_TCE_vec(2:end-1)  = k_Abio_TCE;
k_Abio_C2H2_vec(2:end-1) = k_Abio_C2H2;
k_Abio_C2H4_vec(2:end-1) = k_Abio_C2H4;


% Diffusion coeffient of the species
Daq = [0.90 1.72 1.30 1.50 0.9]*1e-4/86400; % [m^2/s]

% Pore diffusion coefficients [m^2/s]
Dp = Daq*poros^2; % Archie's law (exponent =2)

% Matrix for cell properties (nx x ncomp)
poros_w_mat = poros_w*ones(1,ncomp); % matrix for porosity [-]
RV_w_mat    = RV_w*ones(1,ncomp);    % matrix for porous volume [L]
RV_mat      = RV*ones(1,ncomp);      % matrix for bulk cell volume [L]

%--------------------------------------------------------------------------
% INITIAL CONDITIONS
c = zeros(nx,ncomp); % Aqueous phase concentrations [mol/L] 
s = zeros(nx,ncomp); % Sorbed concentrations [mol/kg_solid]

% BOUNDARY CONDITIONS
cin_TCE = 1100/131.4e3; % solubility limit [mol/L]
c(1,:) = [cin_TCE 0 0 0 cin_TCE]; 

% Initialization of breakthrough curve
BTC   = zeros(0,ncomp);
BTC_J = zeros(0,ncomp);
tvec  = zeros(1,0);

%--------------------------------------------------------------------------
% CALCULATION OF GLOBAL MATRICES
% Matrix for the Mass-transfer coefficients 
[Dmat]=calc_Dmat(Dp,nx,dx,ncomp,RV_w);
% Inlet Boundary condition
r_bc  = zeros(nx,ncomp);

% AUTOMATIC CONSTRUCTION OF GLOBAL MOBILITY MATRIX [nx*ncomp x nx*ncomp]
% Initialization of indices
ii_loc_max=nx*3-2; 
i_glob=zeros(ii_loc_max*ncomp,1);
j_glob=zeros(ii_loc_max*ncomp,1);
value_glob=zeros(ii_loc_max*ncomp,1); 
count=1;   % counter for local mobility matrices
for i=1:ncomp
    for j=1:ncomp
        if i==j
            [i_loc,j_loc,value_loc]=mob_mat(Dmat(:,i),nx,dx,dt);
            i_glob(ii_loc_max*(count-1)+1:ii_loc_max*count)=i_loc+nx*(i-1);
            j_glob(ii_loc_max*(count-1)+1:ii_loc_max*count)=j_loc+nx*(j-1);
            value_glob(ii_loc_max*(count-1)+1:ii_loc_max*count)=value_loc;
            count=count+1;
        end
    end
end
% Global mobility matrix [nx*ncomp x nx*ncomp]
MMOB   = sparse(i_glob,j_glob,value_glob,nx*ncomp,nx*ncomp);

% Global storage matrix for water phase [nx*ncomp x nx*ncomp]
MSTORE = spdiags(RV_w_mat(:),0,nx*ncomp,nx*ncomp);

% Global reaction matrix [nx*ncomp x nx*ncomp]
MREAC  = reacrates_Abio(c,k_Abio_TCE_vec,k_Abio_C2H2_vec,k_Abio_C2H4_vec,nx,RV_w,dt);
%--------------------------------------------------------------------------
% Extraction spatial profiles at different times
t_sp = [1 5 10 20 50 100 300]*86400; % in [s]
t_sp_count=floor(t_sp/dt);


% TIME LOOP
t=0;
time_counter=1;
tcount =1;
tic
% Open figure and delete its content (clf clear current figure)
figure(1);clf
%==========================================================================
% Loop over all timepoints
%==========================================================================
while t<=te

disp(sprintf('Time = %1.2f[d] ->',[t/86400]));

doagain=true;
  while doagain
  doagain=false;

c_adv=c; 
s_adv=s;

% Overall matrix in the Right Hand Side -->> (Mass in water and solid phases in OLD time step)
rhs = c_adv.*RV_w_mat + r_bc + s_adv.*(1-poros_w_mat).*RV_mat*rho_s;
  
% IMPLICIT SOLUTION OF REACTIVE TRANSPORT PROBLEM 
repeat=true;
iter =0;
while repeat
    iter=iter+1; % Iteration counter

% Matrix for ds/sc
deriv_c = c.^(n_Fr_mat-1);
deriv_c(isinf(deriv_c))=0;

% Matrix for sorbed mass in the solid phase (nx x ncomp)
sorp_cap_mat = (1-poros_w_mat).*RV_mat.*rho_s.*K_Fr_mat.*(deriv_c);

% Global storage Matrix for the sorbed phase (nx*ncomp x nx*ncomp)
MSTORE_sorp = spdiags(sorp_cap_mat(:),0,nx*ncomp,nx*ncomp);

% Overall matrix 
M = -MMOB + MSTORE + MSTORE_sorp - MREAC;


% PICCARD LOOP
    cnew=reshape(M\rhs(:),nx,ncomp);
      
    max_error=norm(cnew(:)-c(:),inf);
    if (max_error<1e-10  || iter>100)
        repeat=false;
    end
    c=cnew;
    disp(sprintf('        Iteration=%1.0f, Max. Error=%1.4g',[iter, max_error]));    
end;
% Update time
t=t+dt;
tcount = tcount +1;    
  end

% Update Sorbed concentrations [mol/kg]
s = K_Fr_mat.*(c.^n_Fr);

% Flux at the cell interfaces [mol]
J_mass = reshape(MMOB*c(:),nx,ncomp);     % in [mol]
J      = J_mass*dx./(RV_w_mat*1e-3)/dt;   % in [mol/m2/s]


%-----------Extract concentration matrices for selected times--------------
    if tcount == t_sp_count(1)
        c_1    =c;
        s_1    =s;
    end

    if tcount == t_sp_count(2)
        c_2    =c;
        s_2    =s;
    end

    if tcount == t_sp_count(3)
        c_3    =c;
        s_3    =s;
    end
    
    if tcount == t_sp_count(4)
        c_4    =c;
        s_4    =s;
    end
    
    if tcount == t_sp_count(5)
        c_5    =c;
        s_5    =s;
    end
    
    if tcount == t_sp_count(6)
        c_6    =c;
        s_6    =s;
    end 
    
    if tcount == t_sp_count(7)
        c_7    =c;
        s_7    =s;
    end

%--------------------------------------------------------------------------
% GRAPHICAL OUTPUT
%--------------------------------------------------------------------------
% Transient Plot
figure(1)
for i=1:ncomp
subplot(2,3,i)
plot(x, c(:,i))
xlabel('x [m]')
ylabel('Concentration [mol/L]')
title(name_comp{i})
end
drawnow

    
BTC=[BTC;c(end,:)]; % Update breakthrough curve
BTC_J = [BTC_J;J(end,:)];
tvec=[tvec t];
end

figure(2)
plot(x, c(:,:))
hold on
xlabel('x [m]')
ylabel('Concentration [mol/L]')
legend('TCE','C2H2','C2H4','C2H6','Tracer')

% Breakthrough curves
figure(3)
subplot(2,2,1)
set(gca,'Fontsize',16)
plot(tvec/86400, BTC(:,[1 5]),'Linewidth',2)
xlabel('Time [d]')
ylabel('Concentration [mol/L]')
legend('TCE','Tracer')
title('(a)')
set(gca,'Fontsize',16,'Linewidth',2)

subplot(2,2,2)
plot(tvec/86400, BTC(:,[2 3 4]),'Linewidth',2)
xlabel('Time [d]')
ylabel('Concentration [mol/L]')
legend('Acetylene','Ethene','Ethane')
title('(b)')
set(gca,'Fontsize',16,'Linewidth',2)

subplot(2,2,3)
plot(tvec/86400, BTC_J(:,[1 5]),'Linewidth',2)
xlabel('Time [d]')
ylabel('Flux [mol/m^2/s]')
legend('TCE','Tracer')
title('(c)')
set(gca,'Fontsize',16,'Linewidth',2)

subplot(2,2,4)
plot(tvec/86400, BTC_J(:,[2 3 4]),'Linewidth',2)
xlabel('Time [d]')
ylabel('Flux [mol/m^2/s]')
legend('Acetylene','Ethene','Ethane')
title('(d)')
set(gca,'Fontsize',16,'Linewidth',2)

%%
figure(4)
subplot(2,2,1)
plot(x, c_1(:,1),'b-','Linewidth',2)
hold on
% plot(x, c_2(:,1),'g-','Linewidth',2)
plot(x, c_3(:,1),'r-','Linewidth',2)
% plot(x, c_4(:,1),'c-','Linewidth',2)
plot(x, c_5(:,1),'g-','Linewidth',2)
plot(x, c_6(:,1),'k-','Linewidth',2)
% plot(x, c_7(:,1),'k-','Linewidth',2)
plot(x, c_1(:,5),'b','Linewidth',2,'Linestyle','--')
% plot(x, c_2(:,5),'g','Linewidth',2,'Linestyle','--')
plot(x, c_3(:,5),'r','Linewidth',2,'Linestyle','--')
% plot(x, c_4(:,5),'c','Linewidth',2,'Linestyle','--')
plot(x, c_5(:,5),'g','Linewidth',2,'Linestyle','--')
plot(x, c_6(:,5),'k','Linewidth',2,'Linestyle','--')
% plot(x, c_7(:,5),'k','Linewidth',2,'Linestyle','--')
text(.005,7e-3,{'Solid Lines: TCE';'Dashed Lines: Tracer'},'Fontsize',16)
hold on
xlabel('x [m]')
ylabel('Concentration [mol/L]')
% legend('TCE','Tracer')
title('(a) TCE and Tracer')
set(gca,'Fontsize',16,'Linewidth',2)

subplot(2,2,2)
plot(x, c_1(:,2),'b-','Linewidth',2)
hold on
% plot(x, c_2(:,1),'g-','Linewidth',2)
plot(x, c_3(:,2),'r-','Linewidth',2)
% plot(x, c_4(:,1),'c-','Linewidth',2)
plot(x, c_5(:,2),'g-','Linewidth',2)
plot(x, c_6(:,2),'k-','Linewidth',2)
hold on
xlabel('x [m]')
ylabel('Concentration [mol/L]')
legend('t = 1 [d]','t = 10 [d]','t = 50 [d]','t = 100 [d]')
title('(b) Acetylene')
set(gca,'Fontsize',16,'Linewidth',2)

subplot(2,2,3)
plot(x, c_1(:,3),'b-','Linewidth',2)
hold on
% plot(x, c_2(:,1),'g-','Linewidth',2)
plot(x, c_3(:,3),'r-','Linewidth',2)
% plot(x, c_4(:,1),'c-','Linewidth',2)
plot(x, c_5(:,3),'g-','Linewidth',2)
plot(x, c_6(:,3),'k-','Linewidth',2)
hold on
xlabel('x [m]')
ylabel('Concentration [mol/L]')
% legend('TCE','Tracer')
title('(c) Ethene')
set(gca,'Fontsize',16,'Linewidth',2)


subplot(2,2,4)
plot(x, c_1(:,4),'b-','Linewidth',2)
hold on
% plot(x, c_2(:,1),'g-','Linewidth',2)
plot(x, c_3(:,4),'r-','Linewidth',2)
% plot(x, c_4(:,1),'c-','Linewidth',2)
plot(x, c_5(:,4),'g-','Linewidth',2)
plot(x, c_6(:,4),'k-','Linewidth',2)
hold on
xlabel('x [m]')
ylabel('Concentration [mol/L]')
% legend('TCE','Tracer')
title('(d) Ethane')
set(gca,'Fontsize',16,'Linewidth',2)

