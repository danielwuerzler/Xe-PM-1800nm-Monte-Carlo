close all 
clear all
PlotOpt         = 1;
start_at_0      = 0;
wvlnm           = 1800;
IWcm            = 1e14;
CEP             = 0;
dt              = 0.1;
cutoff          = 0.000001;
z               = 0;
k2prime         = 0;
fwhm_ftl        = 13.1;
n_c             = 6;
  

type = 'gauss';
subplot(2,2,1)
[ tgrid, A, E_real, Env, ftl ] = fct_master_fields_shared( PlotOpt, type, start_at_0,wvlnm, IWcm,...
                                                  CEP, dt,cutoff, z, k2prime, fwhm_ftl, n_c );
tgrid1 = tgrid;
env1 = Env;                                                    
ftl1 = ftl ;                                                        
                                               
type = 'sin2';
subplot(2,2,2)
[ tgrid, A, E_real, Env, pulselength ] = fct_master_fields_shared( PlotOpt, type, start_at_0,wvlnm, IWcm,...
                                                  CEP, dt,cutoff, z, k2prime, fwhm_ftl, n_c );                                                       
tgrid2 = tgrid;
env2 = Env;
ftl2 = ftl ;  

subplot(2,2,3)
plot(tgrid1, env1,'r') 
hold on
xlim([-1000 1000])

plot(tgrid2, env2,'b') 
hold on
xlim([-1000 1000])




