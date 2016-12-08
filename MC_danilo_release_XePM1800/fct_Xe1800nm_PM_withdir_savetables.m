function   fct_Xe1800nm_PM_withdir_savetables( n, N, CEP, kmax, nsample, I0, int_spacing, savename2, cross_section_fname,cutoff_winkel, save_tables, table_name )

% this is the master function to use for the MC Xenon phase meter simulations
% at 18000 nm, some laser parameters need to be changed IN the function and
% aren't passed as a variable, e.g. wvlm, w0, ...

% Input variables:

% n ..... sin^2 cycles on ramp
% N ..... sin^2 cycles on flat top
% CEP ... carrier envelope phase of this pulse
% kmax ..... over how many intensities do you want to average ?
% nsample ..... how many electrons do you want to simulate for EACH int.?
% I0 .... peak intensity of pulse
% int_spacing ... how many intensities do you want to have in the linspace?
% savename2 ... name that you want to give to the intensity averaged 2D PEMD for a single CEP. Needs to be a .mat
% cross_section_fname ... name of the .mat file that contains the cross section.  Needs to be a .mat
% cutoff_winkel ... the angle until which you want to set the cross section constant in DEGREE (used to take care of the divergent forward-scattering)
% save_tables ... do you want to save a table for each new cep ?
% table_name ... which name do you want to give to the table ? It is useful to give it recognizable name, as the table gets moved into the single-cep folder and might be used again later.  Needs to be a .mat
%

% Output:
%
% The function has no direct output. Instead it saves the 2D PEMD histogram and its axes in savename2

% by Max Möller and Danilo Zille

    
% pulse characteristics
wvlm     = 1800e-9;                            %
atom     = 'Xe'; % %atom     = 'He1s2s-singlet';    %'He1s2s-triplet'
s        = 3/2;
pulse    = [n N wvlm];
dt       = 0.1;               % temporaresolution

% focal characteristics
w0    = 0.087;            % w0, [mm]
Zr    = 3.3;              % Zr, [mm]

% get hist of intensity distribution in focus
IWcm_grid_Vavg       = linspace(I0/10,I0,int_spacing);
% IWcm_grid_Vavg       = 1e11:2e11:I0;   % for single intensity!!!
NI_Vavg              = length(IWcm_grid_Vavg);
PlotOpt              =0;
[Iaxis_Vavg,hist_Isample_Vavg] = fct_get_IhistVavg_zR(w0,I0,Zr,IWcm_grid_Vavg,PlotOpt);

% ion yield in igrid from intensity distribution
IonYield = zeros(1,length(hist_Isample_Vavg));
for l=1:1:NI_Vavg-1
    [Yield]  = fct_get_IonYield(atom,Iaxis_Vavg(l),0,pulse);        % CEP = 0 ist a good estimate
    IonYield(l) = Yield;  
end
hist_Isample_Vavg = hist_Isample_Vavg./max(hist_Isample_Vavg);
IonYield          = IonYield./max(IonYield);
Y_Vavg_Iyield     = (hist_Isample_Vavg.*IonYield)./sum(hist_Isample_Vavg.*IonYield);

% prep vol averaged ion yield for big monte carlo
Pnorm       = Y_Vavg_Iyield;
intens_ind  = fct_gen_distr(Pnorm,1,1e3);
Isample     = Iaxis_Vavg(intens_ind); 
hist_Isample= hist(Isample,Iaxis_Vavg);
%   plot(Iaxis_Vavg,hist_Isample)
 
% prepare name1 for this run
savestr =  ['n',num2str(n),'_N',num2str(N),'_',atom,'_',num2str(wvlm./1e-9),'_Ip',num2str(I0,'%10.1e\n')];

resol = 0.02; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV

Up       = (0.09337.*I0.*(wvlm.^2))./27.211;   % aumu
% set up the size of the histogram
y_lim = sqrt(2*10*Up);
x_lim = sqrt(2*3.2*Up);
vpagrid_01 = -y_lim-0.1:resol:y_lim+0.1;
vpegrid_01 = -x_lim-0.1:resol:x_lim+0.1;


hist_vyvz_resc_loop_01 = zeros(length(vpegrid_01),length(vpagrid_01));
hist_vyvz_dir_loop_01 = zeros(length(vpegrid_01),length(vpagrid_01));

% savename2    = atom;

for k=1:kmax % loop over intensities
    
    tic
    Pnorm          = Y_Vavg_Iyield;
    intens_ind     = fct_gen_distr(Pnorm,1,1);

    
    if intens_ind < max(intens_ind)
        Isample(k)  = Iaxis_Vavg(intens_ind) + rand(1,1).*( Iaxis_Vavg(100) - Iaxis_Vavg(99) );
    else
        Isample(k)  = Iaxis_Vavg(intens_ind);
    end
    
%     xsample = Vpagrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1))/2;
%     ysample = Vpegrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1))/2;
    
    
    CEPsample      = CEP;
    savename1      = [num2str(k),'_',savestr,'.mat'];
    
    % get the momenta for one set of CEP and intensity
    
    %nsample = 1e5;
        [vxf_all vyf_all vzf_all vxf_dir vyf_dir vzf_dir vxf_resc vyf_resc vzf_resc T] = fct_singleICEP_v0(k, Isample(k),CEPsample,pulse,atom,nsample,savename1,dt,s,...
                                                                                        cross_section_fname, cutoff_winkel, save_tables, table_name);
    % histogram for current i & cep averaged results
    PlotOpt = 0;
    [hist_vxvy_resc hist_vxvz_resc hist_vyvz_resc_01] = fct_illstr_Projections(vpagrid_01,vpegrid_01,vxf_resc,vyf_resc,vzf_resc,PlotOpt);
    [hist_vxvy_dir hist_vxvz_dir hist_vyvz_dir_01] = fct_illstr_Projections(vpagrid_01,vpegrid_01,vxf_dir,vyf_dir,vzf_dir,PlotOpt);

    hist_vyvz_resc_loop_01 = hist_vyvz_resc_loop_01 + hist_vyvz_resc_01;
    hist_vyvz_dir_loop_01 = hist_vyvz_dir_loop_01 + hist_vyvz_dir_01;
    

display(strcat('Volume intensity ',num2str(k),' done and PEMD saved!'));

   % %     savename2    = [num2str(k),'_all',savestr,'.mat'];
       save(savename2, 'vpagrid_01','vpegrid_01','hist_vyvz_resc_loop_01','hist_vyvz_dir_loop_01');
    toc   
end


toc
foldername = ['n_',num2str(n),'_N_',num2str(N),'_wvl_',num2str(wvlm./1e-9),'_I0_',num2str(I0,'%10.1e\n'),...
             '_dt_',num2str(dt),'_kmax_',num2str(kmax),'_nsample_',num2str(nsample,'%10.0e\n')];
         
mkdir(foldername);
movefile(savename2,foldername)  
movefile(table_name,foldername)  


end
