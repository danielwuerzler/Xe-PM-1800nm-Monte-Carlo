close all
clear all

atom = 'C:\Users\Zhang\Desktop\matlab_git\Xe-PM-1800nm-Monte-Carlo\MC_danilo_release_XePM1800\n_2_N_0_wvl_1800_I0_1.0e+14_dt_0.1_kmax_120_nsample_5e+06\1800nm_Xe_8.7fs_25_25ceps_CEP.mat';

c0       = 299792458./(2.1876912633e6);        % au, speed of light 
wvlm     = 1800e-9;                            % 
wvl      = wvlm/(5.2917720859e-11);                            % [au]
omega    = (2*pi*c0)/wvl;                                      % [au], angular frequency
T        = 2*pi/omega; 
IWcm     = 1e14;   
Up       = (0.09337.*IWcm.*wvlm.^2)./27.211;   % au     

resol = 0.01; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV


load(atom);
singlet_resc =  hist_vyvz_resc_loop_01;
singlet_dir =  0;%hist_vyvz_dir_loop_01;
singlet = singlet_resc + singlet_dir;

sum(sum(singlet))
vpagrid = vpagrid_01;
vpegrid = vpegrid_01;
subplot(1,2,1) % first not interpolates
surf(vpegrid,vpagrid,log10(singlet')); axis square; shading flat; view(2);


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

imagesc(vpegrid,vpagrid,log10(singlet)); axis square; shading flat; view(2);


nsample = 1e6;

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