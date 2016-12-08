close all
clear all

c0       = 299792458./(2.1876912633e6);        % au, speed of light 
wvlm     = 1800e-9;                            % 
wvl      = wvlm/(5.2917720859e-11);                            % [au]
omega    = (2*pi*c0)/wvl;                                      % [au], angular frequency
T        = 2*pi/omega; 
IWcm     = 1e14;   
Up       = (0.09337.*IWcm.*wvlm.^2)./27.211;   % au  
resol = 0.01

for ii = 1 : 2
    
    
    file_name = strcat('C:\Users\Zhang\Desktop\matlab_git\Xe-PM-1800nm-Monte-Carlo\MC_danilo_release_XePM1800\n_1_N_0_wvl_1800_I0_1.0e+14_dt_0.1_kmax_300_nsample_1e+06\1800nm_Xe_4.4fs_',num2str(ii),'_25ceps_CEP.mat');
    load(file_name)
    
    singlet_resc =  hist_vyvz_resc_loop_01;
    singlet_dir =  0;%hist_vyvz_dir_loop_01;
    singlet = singlet_resc + singlet_dir;

    sum(sum(singlet));
    vpagrid = vpagrid_01;
    vpegrid = vpegrid_01;
%     subplot(1,2,1) % first not interpolates
%     surf(vpegrid,vpagrid,log10(singlet')); axis square; shading flat; view(2);
    
    
    
    
    % get the emission angle for each PEMD pixel
    for i=1:length(vpagrid)
        for j=1:length(vpegrid)
    ang(i,j) = atan(vpegrid(j)./vpagrid(i));
        end
    end

    cutoff_angle = 5; % opening of the cone, in DEG
    ang(abs(ang)>=deg2rad(cutoff_angle))=0;       
    % create angle mask
    for i=1:length(vpagrid)
        for j=1:length(vpegrid)
    if ang(i,j) ~=0
        ang(i,j) = 1;
    end
        end
    end

    singlet = ang.*singlet';   % multiply mask

%     imagesc(vpegrid,vpagrid,log10(singlet)); axis square; shading flat; view(2);
    

    nsample = 5e6;

    % make the 2d histogram to a list 
    P =  singlet(:);
    P = P./sum(P);

    ind  = fct_gen_distr(P',1,nsample);

    [Vpegrid,Vpagrid] = meshgrid(vpegrid,vpagrid);
    %hist_v0xv0y = hist3([vxf vyf], {vpegrid' vpegrid'});

    xsample = Vpagrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1))/2;
    ysample = Vpegrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1))/2;


    %% energy of intersting particles

    Esample = 0.5.*(xsample.^2+ysample.^2) * 27.211;

    Eaxis = 0:0.1:200;
    histEaxis = hist(Esample,Eaxis);

    semilogy(Eaxis,histEaxis)
%     set(ylim([1e2 1e5]))
%     set(gca,'ylim',[1e2 1e5]);
%     
    
    M(ii)=getframe(gcf);
    

end
movie2avi(M,'4.4fs_ati.avi','FPS',2)

