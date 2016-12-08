clear all
close all



folder_name     = 'n_3_N_0_wvl_1800_I0_7.0e+13_dt_0.1_kmax_300_nsample_1e+06';
cd(strcat(folder_name))

for ii = 1:25
    
    file_name = strcat(num2str(ii),'_25ceps_','ATI.mat')
    load(file_name)
    

    semilogy(E_axis_ati, ati_spec_l, E_axis_ati, ati_spec_r); 

    M(ii)=getframe(gcf);
     
    
end
movie2avi(M,'13.1fs 0.7e14 I0 ATI.avi','FPS',2)