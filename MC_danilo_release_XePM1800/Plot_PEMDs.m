%% here plot the log 10 vmi map after adding up all single dstribution electron loops
%%  first sum up all the k loops
anfang =1;
ende = 61; %largest run number
clearvars -except anfang ende
a=0;
b=0;
c=0;
d=0;
for ii=anfang:ende
    cd(strcat('new_s1.5_n7N0wvl1800I6.0e+12_dt0.1kmax25nsample2e+07_run',num2str(ii)));
    load(strcat('sodium.mat'));
    a = a + hist_vyvz_resc_loop_01;
    b = b + hist_vyvz_dir_loop_01;
cd('..')
end
hist_vyvz_resc_loop_01 = a;
hist_vyvz_dir_loop_01 = b;

    save('sodium.mat');
    
    %% single plot
atom = '1800nm_Xe_8.7fs_2_25_ceps_CEP.mat';
% these parameters are just for fun...
c0       = 299792458./(2.1876912633e6);        % au, speed of light 
wvlm     = 1800e-9;                            % 
wvl      = wvlm/(5.2917720859e-11);                            % [au]
omega    = (2*pi*c0)/wvl;                                      % [au], angular frequency
T        = 2*pi/omega; 
IWcm     = 6e12;   
Up       = (0.09337.*IWcm.*wvlm.^2)./27.211;   % au     
pmax = sqrt(2*10*Up)

resol = 0.01; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV


load(atom);
singlet =  hist_vyvz_resc_loop_01;
vpagrid = vpagrid_01;
vpegrid = vpegrid_01;

%     singlet = singlet/sum(sum(singlet)); %normalize
singlet = singlet/max(max(singlet));


surf (vpegrid,vpagrid,log10(( singlet'))); shading interp
view(2)
   xlabel('vy'); ylabel('vz'); title('resc');
    axis equal
    title('Singlet')
    colorbar
    drawnow    
%  ylim([-0.5 0.5])   
%  xlim([-lim_correction*y_lim,lim_correction*y_lim]);
    caxis([-7,0])

set(gca,'fontsize',16)
xlhand = get(gca,'xlabel');
set(xlhand,'string','P_x','fontsize',16) ;
ylhand = get(gca,'ylabel');
set(ylhand,'string','P_z','fontsize',16) ;
title(strcat('Rescatteredt Yield,','resolution = ', num2str(resol)))
set(gcf,'color','w'); % set figure background to white
axis tight

