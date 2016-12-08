close all 
clear all

wvlm            = 1800e-9; 
I0              = 8e13;   
nr_ceps         = 25;
cutoff_angle    = 5; % opening of the cone, in DEG
save_ati_switch = 0;

pulse_length    = 13.1; 
folder_name     = 'Xe_1800nm_8.0e+13Wcm2_300kmax_25ceps_1e+06nsample_13.1fs_0mm';
a = 0;
b = 0;

for ii = 1:nr_ceps            %%   n=1, start at ii =2, cause cep 1 only has half of the data.

%     file_name     = strcat(folder_name,'\1800nm_Xe_',...
%                         num2str(pulse_length),'fs_',num2str(ii),'_50ceps_CEP.mat');
    
     file_name     = strcat(folder_name,'\',num2str(ii),'_25ceps_CEP.mat');  % new field
    
    load(file_name)
        
    fct_plot_ati_single_CEP(wvlm, I0, cutoff_angle,ii, save_ati_switch, nr_ceps,...
                            folder_name, hist_vyvz_resc_loop_01, hist_vyvz_dir_loop_01,...
                            vpagrid_01, vpegrid_01)                 
end
%% CEP average
for k = 1:nr_ceps
    
    file_name = strcat(folder_name,'\',num2str(k),'_',num2str(nr_ceps),'ceps_','ATI.mat')
    load(file_name)
    
    a = a + ati_spec_l;
    b = b + ati_spec_r;
    display(num2str(k),'done')
end

ati_spec_l      = a;
ati_spec_r      = b; 

if save_ati_switch == 1
    savename_ati_total  = strcat(folder_name,'\total_ceps_','ATI.mat'); 
    save(savename_ati_total, 'E_axis_ati','ati_spec_l','ati_spec_r');    
end

figure;
semilogy(E_axis_ati, a,'r.'); hold on
semilogy(E_axis_ati, ati_spec_r,'b.');
ylim([1e-5 1e2])
grid on

