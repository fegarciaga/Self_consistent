%  A script to set the input parameters and run a CPMC calculation
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% �2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% system parameters:
Lx=4; % The number of lattice sites in the x direction
Ly=4; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction

N_sites=Lx*Ly*Lz;

N_up=5; % The number of spin-up electrons
N_dn=5; % The number of spin-down electrons

kx=+0.0; % The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=+0.0; % The y component of the twist angle in TABC
kz=0; % The z component of the twist angle in TABC

U=4; % The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction

tz2=[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2];
tx2=0.0;
ty2=0.0;
N_run=length(tz2);
%% run parameters:
deltau=0.01; % The imaginary time step
N_wlk=1000; % The number of random walkers
N_blksteps=40; % The number of random walk steps in each block
N_eqblk=30; %The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=100; % The number of blocks used in the measurement phase
N_max=5;
N_iter=20;
itv_modsvd=1; % The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=5; % The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=40; % The interval between two adjacent energy measurements
t_bp=81;
t_pop=1;
suffix=datestr(now,'_yymmdd_HHMMSS'); % time stamp for the saved *.mat filename. Can be changed to any desired string 
%% Create array for potential energy
E_ave_fe=zeros(N_run,1);
E_err_fe=zeros(N_run,1); 
N_up_ave=zeros(N_sites,N_run);
N_dn_ave=zeros(N_sites,N_run);
N_updn=zeros(N_run);

for ii=1:N_run
    suffix=strcat('_try',int2str(tz2(ii)));
    [E_ave_fe(ii),E_err_fe(ii),N_up_ave(:,ii), N_dn_ave(:,ii),N_updn(ii), savedFile]=SCCPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,tx2,ty2,tz2(ii),deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em, t_bp, t_pop, N_iter, N_max, suffix);
end

writematrix(N_updn,'N_updn.txt','Delimiter','tab');
%% post-run:
% load saved data into workspace for post-run analysis:
%load (savedFile);
%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"