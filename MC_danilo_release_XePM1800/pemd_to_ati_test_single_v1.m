close all
% clear all

c0       = 299792458./(2.1876912633e6);        % au, speed of light 
wvlm     = 1800e-9;                            % 
wvl      = wvlm/(5.2917720859e-11);                            % [au]
omega    = (2*pi*c0)/wvl;                                      % [au], angular frequency
T        = 2*pi/omega; 
I0       = 1e14;   
Up       = (0.09337.*I0.*wvlm.^2)./27.211;   % au  
    
singlet_resc =  hist_vyvz_resc_loop_01;
singlet_dir =  hist_vyvz_dir_loop_01;
singlet = singlet_resc + singlet_dir;

sum(sum(singlet));
    
vpagrid = vpagrid_01;
vpegrid = vpegrid_01;
%     subplot(1,2,1) % first not interpolates
%     surf(vpegrid,vpagrid,log10(singlet')); axis square; shading flat; view(2);
i1 = length(vpagrid);     
%%   Right side   
% get the emission angle for each PEMD pixel
for i=1:length(vpagrid)
    for j=1:length(vpegrid)
        if i<= length(vpagrid)/2
            ang(i,j) = atan(vpegrid(j)./vpagrid(i));
        else
            ang(i,j) = 0;
        end
    end
end

cutoff_angle = 5; % opening of the cone, in DEG
ang(abs(ang)>=deg2rad(cutoff_angle))=0;       
% create angle mask
for i=1:length(vpagrid)
    for j=1:length(vpegrid)
        if ang(i,j) ~= 0
            ang(i,j) = 1;
        end
    end
end

singlet_r = ang.*singlet';   % multiply mask

subplot(2,2,1)
imagesc(vpegrid,vpagrid,log10(singlet_r)); axis square; shading flat; view(2);

subplot(2,2,2)
tot_spec_r_ana = sum(singlet_r')./abs((vpagrid));
semilogy(27.211*0.5*vpagrid.^2,tot_spec_r_ana./max(tot_spec_r_ana),'r');
grid on

%%   Left side   
% get the emission angle for each PEMD pixel
for i=1:length(vpagrid)
    for j=1:length(vpegrid)
        if i> length(vpagrid)/2
            ang(i,j) = atan(vpegrid(j)./vpagrid(i));
        else
            ang(i,j) = 0;
        end
    end
end

ang( abs(ang)>= deg2rad(cutoff_angle) )=0;       

% create angle mask
for i=1: length(vpagrid)
    for j=1:length(vpegrid)
        if ang(i,j) ~= 0
            ang(i,j) = 1;
        end
    end
end

singlet_l = ang.*singlet';   % multiply mask

subplot(2,2,3)
imagesc(vpegrid,vpagrid,log10(singlet_l)); axis square; shading flat; %view(2);

subplot(2,2,4)
tot_spec_l_ana = sum(singlet_l')./abs((vpagrid));
semilogy(27.211*0.5*vpagrid.^2,tot_spec_l_ana./max(tot_spec_l_ana),'k');
grid on


ind = find(27.211*0.5*vpagrid.^2 > 1);
figure;
semilogy(27.211*0.5*vpagrid(ind).^2,tot_spec_l_ana(ind)./max(tot_spec_l_ana(ind)),'r.'); hold on
semilogy(27.211*0.5*vpagrid(ind).^2,tot_spec_r_ana(ind)./max(tot_spec_r_ana(ind)),'k.'); hold on
%ylim([-3 0])
grid on



 
