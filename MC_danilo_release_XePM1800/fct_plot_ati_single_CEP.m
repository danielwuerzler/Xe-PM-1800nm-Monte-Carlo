function    fct_plot_ati_single_CEP(wvlm, I0, cutoff_angle,ii, save_ati_switch,nr_ceps,...
            folder_name, hist_vyvz_resc_loop_01,hist_vyvz_dir_loop_01, vpagrid_01, vpegrid_01 )


singlet_resc =  hist_vyvz_resc_loop_01;
singlet_dir  =  hist_vyvz_dir_loop_01;
singlet = singlet_resc + singlet_dir;

sum(sum(singlet));

vpagrid = vpagrid_01;
vpegrid = vpegrid_01;
    
%%   Right side   % get the emission angle for each PEMD pixel
for i=1:length(vpagrid)
    for j=1:length(vpegrid)
        if i<= length(vpagrid)/2
            ang(i,j) = atan(vpegrid(j)./vpagrid(i));
        else
            ang(i,j) = 0;
        end
    end
end
    
ang(abs(ang)>=deg2rad(cutoff_angle))=0;   % create angle mask

for i=1:length(vpagrid)
    for j=1:length(vpegrid)
        if ang(i,j) ~= 0
            ang(i,j) = 1;
        end
    end
end
    
singlet_r  = ang.*singlet';   % multiply mask
E_axis_ati = 27.211*0.5*vpagrid.^2;
ati_spec_r = sum(singlet_r')./abs((vpagrid));
ati_spec_r = ati_spec_r./max(ati_spec_r);
      
%     subplot(2,2,1)
%     imagesc(vpegrid,vpagrid,log10(singlet_r)); axis square; shading flat; view(2);
%     subplot(2,2,2) 
%     semilogy(E_axis_ati, ati_spec_r,'r');
%     grid on

%%   Left side      % get the emission angle for each PEMD pixel
for i=1:length(vpagrid)
    for j=1:length(vpegrid)
        if i> length(vpagrid)/2
            ang(i,j) = atan(vpegrid(j)./vpagrid(i));
        else
            ang(i,j) = 0;
        end
    end
end

ang( abs(ang)>= deg2rad(cutoff_angle) )=0;      % create angle mask 

for i=1: length(vpagrid)
    for j=1:length(vpegrid)
        if ang(i,j) ~= 0
            ang(i,j) = 1;
        end
    end
end

singlet_l  = ang.*singlet';   % multiply mask
E_axis_ati = 27.211*0.5*vpagrid.^2;
ati_spec_l = sum(singlet_l')./abs((vpagrid));
ati_spec_l = ati_spec_l./max(ati_spec_l);
      
%     subplot(2,2,1)
%     imagesc(vpegrid,vpagrid,log10(singlet_r)); axis square; shading flat; view(2);
%     subplot(2,2,2) 
%     semilogy(E_axis_ati, ati_spec_r,'r');
%     grid on

    
if save_ati_switch == 1
    
    savename_ati  = strcat(folder_name,'\',num2str(ii),'_',num2str(nr_ceps),'ceps_','ATI.mat'); 
    save(savename_ati, 'E_axis_ati','ati_spec_l','ati_spec_r');    

end


end