%% here plot the log 10 vmi map after adding up all single dstribution electron loops

clear all
close all

%% Single CEP plot

c0       = 299792458./(2.1876912633e6);        % au, speed of light 
wvlm     = 1800e-9;                            % 
wvl      = wvlm/(5.2917720859e-11);            % [au]
omega    = (2*pi*c0)/wvl;                      % [au], angular frequency
T        = 2*pi/omega; 
IWcm     = 1.2e14;   
Up       = (0.09337.*IWcm.*wvlm.^2)./27.211;   % au     
pmax     = sqrt(2*10*Up);

% resol    = 0.01;                               % 0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV

cd(strcat('n_2_N_0_wvl_1800_I0_1.0e+14_dt_0.1_kmax_120_nsample_5e+06'))

for i = 1 : 25
    
    
    file_name = strcat('1800nm_Xe_8.7fs_',num2str(i),'_25ceps_CEP.mat');
    load(file_name)
    
    singlet =  hist_vyvz_resc_loop_01;
    singlet = singlet/max(max(singlet));       % Normalization
    singlet = singlet/sum(sum(singlet)); %normalize

    vpagrid = vpagrid_01;
    vpegrid = vpegrid_01;
    
    
    surf (vpegrid,vpagrid,log10(( singlet'))); shading interp
    view(2)
    xlabel('vy'); ylabel('vz'); title('resc');
    axis equal
    title(strcat('CEP Step ',num2str(i)))
    colorbar
    drawnow    
%  ylim([-0.5 0.5])   
%  xlim([-lim_correction*y_lim,lim_correction*y_lim]);
%             caxis([-5,0])

    set(gca,'fontsize',16)
    xlhand = get(gca,'xlabel');
    set(xlhand,'string','P_x','fontsize',16) ;
    ylhand = get(gca,'ylabel');
    set(ylhand,'string','P_z','fontsize',16) ;
    title(strcat('Rescatteredt Yield,','CEP steps ',num2str(i)))
    set(gcf,'color','w'); % set figure background to white
    axis tight

    M(i)=getframe(gcf)
end

movie2avi(M,'8.7fs .avi','FPS',2)






