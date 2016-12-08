close all
clear all

%%  Pulse parameters
n           = 3;                % n ..... sin^2 cycles on ramp
N           = 0;                % N ..... sin^2 cycles on flat top
I0          = 0.8e14;           % I0 .... peak intensity of pulse
wvlm        = 1800e-9;          % wvlm .. center wavelength

%% Calculation parameters
kmax        = 300;                % kmax ..... over how many intensities do you want to average ?
nsample     = 0.9e6;              % nsample ..... how many electrons do you want to simulate for EACH int.?
int_spacing = 300;             % int_spacing ... how many intensities do you want to have in the linspace?

%% Cross section para meters
atom        = 'Xe';
cutoff_winkel = 20;         % cutoff_winkel ... the angle until which you want to set the cross section constant in DEGREE (used to take care of the divergent forward-scattering)
cross_section_fname = 'xe_cs.mat';

%% define the field
IWcm      = I0;
dt        = 0.1;
CEP       = 0;                % changable in the loop
[Env E I A ALPHA BETA v r tgrid E_fh] = fct_get_EnvEIAAlphBetvr_sin2_nN(dt,wvlm,IWcm,CEP,n,N);
ind = find(I>0.5.*max(I));
tp  = tgrid(max(ind))-tgrid(min(ind)); 
tp_fs =round(tp.*24.2e-18/1e-15,1);

%% Return time tables
save_tables = 1;            % save_tables ... 1 switch on to save tables, 0 tables exists.

%% Main CEP function
nr_ceps = 50; 

CEP_vec = linspace(0, 2*pi, nr_ceps);

for j = 1:nr_ceps
    CEP = CEP_vec(j);
    
    table_name = strcat(num2str(wvlm*10^9),'nm_',num2str(atom),'_',num2str(tp_fs),'fs_',...
                    num2str(j),'_',num2str(nr_ceps),'ceps_','Table.mat'); 
    savename2  = strcat(num2str(wvlm*10^9),'nm_',num2str(atom),'_',num2str(tp_fs),'fs_',...
                    num2str(j),'_',num2str(nr_ceps),'ceps_','CEP.mat'); 
           
    fct_Xe1800nm_PM_withdir_savetables( n, N, CEP, kmax, nsample, I0,...
                 int_spacing, savename2, cross_section_fname, cutoff_winkel, save_tables, table_name )
end


