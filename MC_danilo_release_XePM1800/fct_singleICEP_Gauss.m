function [vxf_all vyf_all vzf_all vxf_dir vyf_dir vzf_dir vxf_resc vyf_resc vzf_resc T] = fct_singleICEP_Gauss(int_index, IWcm,CEP,pulse_,atom,nsample,savename,dt,s,...
                                                                                                           cross_section_fname, cutoff_winkel, save_tables, table_name,...
                                                                                                           type, start_at_0, wvlnm, cutoff, z, kdoubleprime, fwhm_nochirp, n_c) 
% int_index is used to create a new table a the beginning of intensity averaging loop

%% parameters for the laser pulse 
c0       = 299792458./(2.1876912633e6);          % au, speed of light 
n        = pulse_(1);                            % propto pulse duration
N        = pulse_(2);
wvlm     = pulse_(3);                            % 
wvl      = wvlm/(5.2917720859e-11);              % [au]
omega    = (2*pi*c0)/wvl;                        % [au], angular frequency
T        = 2.*pi./omega;

%% define the field
PlotOpt   = 0;
[ tgrid, A, E, Env, pulselength ] = fct_master_fields_shared( PlotOpt, type, start_at_0, wvlnm,...
                                         IWcm, CEP, dt, cutoff, z, kdoubleprime, fwhm_nochirp, n_c );
 
    A_interp = [0 A];
    for i = 1:length(A)
        A_interp(i) = (A_interp(i)+A_interp(i+1))/2;
    end
    A_interp = A_interp(1:1:end-1);
    A        = [0 A(1:1:end-1)];

    % integration vector potential from -infinity to t
    ALPHA = cumsum(A_interp(:))*dt;
    ALPHA = [0 ; ALPHA(1:1:end-1)];

    %integration (vector potential)� from -infinity to t
    BETA = cumsum(A_interp.^2)*dt;
    BETA = [0 BETA(1:1:end-1)];

    %trajectorie of electron released at time tr
    v = @(tr_index) A - A(tr_index);
    r = @(tr_index) ALPHA - ALPHA(tr_index)-A(tr_index).*(tgrid'-tgrid(tr_index));
    
%% function handle for bsi
Fc = @(kappa_,Z_,m_) (kappa_.^4)./(8.*(2.*Z_-kappa_.*(m_+1)));

%%  All atoms, pure rate, hydrogen rate, gen atom, wo_BI, wo_me
SaturationSwitch = 1;
PlotOpt          = 0;
tbiSwitch        = 1;
AtomString       = atom;
[Ip kappa Z Cnl l beta alphaN alphaI] = fct_get_Atom(AtomString);
Up               = (0.09337.*IWcm.*(wvlm.^2))./27.211;   % au
gamma            = sqrt(abs(Ip)/(2*Up));

Ip               = Ip.*ones(1,length(E));
Z                = Z.*ones(1,length(E));
Cnl              = Cnl.*ones(1,length(E));
l                = l*ones(1,length(E));
beta             = beta*ones(1,length(E));
m                = 0*ones(1,length(E));
alphaN           = alphaN.*ones(1,length(E));
alphaI           = alphaI.*ones(1,length(E));
%Ip    = Ip-0.5.*(alphaN-alphaI).*E.^2;
Fc_GenA_H1s_wo_tbi_wo_me  = Fc(kappa,Z,m);
[IonAmp_GenA_H1s_wo_tbi_wo_me Yield_GenA_H1s_wo_tbi_wo_me] = fct_TolRate_GenAtom_TBIcor(E,tgrid,Ip,Z,Cnl,l,m,beta,SaturationSwitch,tbiSwitch,PlotOpt);
IonAmp = IonAmp_GenA_H1s_wo_tbi_wo_me;
IonAmp=IonAmp/max(IonAmp);

%% settings for the monte carlo
% nsample = 1e3;

%% generate starting times
PlotOpt = 0;
[texit_ind texit] = fct_gen_texit(nsample,tgrid,IonAmp,PlotOpt);
hist_texit =  hist(texit,tgrid);
%texit_ind = 1:1:(length(tgrid)-10);

%% all trajectories - set tunneling conditions for all trajectories
Eexit_sample  = E(texit_ind);                                    % field at tunneling
Ip_sample     = Ip(texit_ind);                                   % Ip at tunneling
vpexit        = zeros(1,length(texit));                          % perp velocity
[chi_sample]  = fct_get_angexit(vpexit,PlotOpt);                 % emission angle
vpaxit        = zeros(1,length(texit));                          % para velocity
v0x_sample    = vpexit.*cos(chi_sample);
v0y_sample    = vpexit.*sin(chi_sample);
v0z_sample    = vpaxit;
rpaxit        = zeros(1,length(texit));                          % parallel tunnel exit
rpexit        = zeros(1,length(texit));                          % perp tunnel exit
r0x_sample    = rpexit.*cos(chi_sample);
r0y_sample    = rpexit.*sin(chi_sample);
r0z_sample    = rpaxit;

% all trajectories - condense initial conditions in matrix
[texit v0mat r0mat] = fct_get_initMat(texit,v0x_sample,v0y_sample,v0z_sample,r0x_sample,r0y_sample,r0z_sample);

% all trajectories -  illustrate tunneling 
PlotOpt= 0;
% % v0grid = -0.5:0.005:0.5;
% % r0grid = -15:0.005:15;
%fct_illstr_Tunnel(v0grid,r0grid,v0mat,r0mat,PlotOpt)

% dice an actual return according to the scattering prob and so on
PlotOpt = 0;
%s       = 3/2;
sigwp   = @(ts,tr)   (tr-ts).^(-s); %sigwp   = @(ts,tr)   (tr-ts).^(-s);
% sigtg   = @(q,theta) abs(q);
% sigtg   = @(q,theta) abs(q); %volumeelement for filled circle

%% cross-section for noble gas potential?

% now load the numerical cross sections
load(cross_section_fname);
% the first few angles are divergent, so set them all to some constant
% value ... play around and check how this value effects your pictures

inds = find(Winkel<=cutoff_winkel);
for ind=1:length(cs(:,1))   % energy,winkel
   cs(ind,inds) = cs(ind,max(inds)+1); %für alle energies, set value to constant
end
% flip the winkel vec to go from -pi to 0, then add together
winkel_total(1:length(Winkel)) = -flipud(Winkel);
winkel_total(length(Winkel)+1:2*length(Winkel)-1) = Winkel(2:length(Winkel));
% do the same with the cs
for ind=1:length(cs(:,1))    % für alle energien flippe
cs_total(ind,1:length(Winkel)) = fliplr(cs(ind,:));
cs_total(ind,length(Winkel)+1:2*length(Winkel)-1) = cs(ind,2:length(Winkel));
end
% semilogy(winkel_total,cs_total(80,:)) % plot the cs at energy index 80 to check if everything is cool

      sigtg =@(q,theta)  fct_na_cs_fortran( q,theta,energy,cs_total );    % fct handle for numeric cross section

    rc   = 0.5871;   % 31e-12 m
    rcov = rc;
    
%% %% find return list

PlotOpt = 0;
NrMax   = 5; % nr of returns
texit_ind = 1:1:length(tgrid);%2.65e4:1:2.93e4;%1:1:length(tgrid);

if(save_tables == 1)
    if(int_index == 1)
        
        [trind_mat tr_mat tt_mat rescind dirind] = fct_get_return1Traj1D_v2(tgrid,texit_ind,zeros(1,length(tgrid)),zeros(1,length(tgrid)),r,rcov,NrMax,PlotOpt);
        toc
        save(table_name,'trind_mat','tr_mat','tt_mat','rescind','dirind');
    else
        load(table_name);
    end
end

disp('found return list')
        
%[trind_mat tr_mat tt_mat rescind dirind] = fct_get_return1Traj1D_v1(tgrid,texit_ind,zeros(1,length(tgrid)),zeros(1,length(tgrid)),r,v,NrMax,PlotOpt);
%% scattering & nonscattering
Pne_all        = IonAmp./sum(IonAmp);
te_ind_all     = fct_gen_distr(Pne_all,1,nsample);
te_sample_all  = tgrid(te_ind_all);                         % list of all starting times
Ne_all         = hist(te_sample_all,tgrid);

nsample_resc = sum(Ne_all(texit_ind(rescind)));                        % number of resc electrons
nsample_dir  = sum(Ne_all(texit_ind(dirind)));                         % number of direct electrons

Pne_dir    = zeros(1,length(tgrid));
Pne_dir(texit_ind(dirind))    = Pne_all(texit_ind(dirind));
te_ind_dir = fct_gen_distr(Pne_dir./sum(Pne_dir),1,nsample_dir);
te_sample_dir = tgrid(te_ind_dir);        

%% distribute rescattering
PlotOpt =0;
[vx_resc_list vy_resc_list vz_resc_list te_resc_list tr_resc_list tt_resc_list vr_resc_list th_resc_list trind_resc_list] = fct_gen_resc2D_fromList_v1_fortran(winkel_total,A,IonAmp,tgrid,rescind,texit_ind(rescind),trind_mat,nsample_resc,sigwp,sigtg,PlotOpt);

% calculate final velocity of rescattered
PlotOpt  = 0;                    
[vxf_resc vyf_resc vzf_resc] = fct_get_vfinal(A,trind_resc_list',vx_resc_list,vy_resc_list,vz_resc_list,Up,PlotOpt); 

%
texit_resc = te_resc_list;
tr_resc    = tr_resc_list;
angle_resc = th_resc_list;
vr_resc    = vr_resc_list;
vxsc       = vx_resc_list;
vysc       = vy_resc_list;
vzsc       = vz_resc_list;
vxf_resc   = vxf_resc';  
vyf_resc   = vyf_resc';
vzf_resc   = vzf_resc';

% vpegrid = -3:0.025:3;
% vpagrid = -3:0.025:3;
% [hist_vxvy_all hist_vxvz_all hist_vyvz_all]    = fct_illstr_Projections(vpagrid,vpegrid,vxf_resc',vyf_resc',vzf_resc',1);

%% direct trajectories
% here comes an option, either you take all dirind only or you take
% boths... dirind & rescind for this

texit_ind_dir= te_ind_dir';        % index on tgrid of the starting time
texit_dir    = te_sample_dir;      % starting time of direct 

% dice perp velocities
MeanOpt = 1;
PlotOpt = 0;
[vpexit_ind_dir vpexit_dir vpgrid] = fct_get_vpexit(Eexit_sample(texit_ind_dir),Ip_sample(texit_ind_dir),MeanOpt,PlotOpt);
 
% set/ calculate para velocity 
PlotOpt = 0; 
[vpaxit_dir] = fct_get_vpaexit(vpexit_dir,PlotOpt);

% % dice the emission angle 
PlotOpt = 0;
[chi_dir] = fct_get_angexit(vpexit_dir,PlotOpt);

% set the initial velocities
vsx_dir = (vpexit_dir.*cos(chi_dir))';                 
vsy_dir = (vpexit_dir.*sin(chi_dir))';                 
vsz_dir = (vpaxit_dir)';

% calculate the final velocity of direct electrons
PlotOpt = 0;
[vxf_dir vyf_dir vzf_dir] = fct_get_vfinal(A,texit_ind_dir,vsx_dir',vsy_dir',vsz_dir',Up,PlotOpt); 

%% results direct
% % % % texit_dir;  % starting times of direct
% % % % vxf_dir; 
% % % % vyf_dir; 
% % % % vzf_dir;

%% all thogether
vxf_dir = vxf_dir';
vyf_dir = vyf_dir';
vzf_dir = vzf_dir';


vxf_all = [vxf_dir;vxf_resc];
vyf_all = [vyf_dir;vyf_resc];
vzf_all = [vzf_dir;vzf_resc];

%% save everything for potential use later... 
% % % % save(savename,...
% % % %                 'tgrid','E','I','wvlm',...                                                     % the field somehow
% % % %                 'texit_dir','vxf_dir','vyf_dir','vzf_dir',...                           % directs
% % % %                 'texit_resc','tr_resc','angle_resc','vr_resc','vxsc','vysc','vzsc',...  % rescattering
% % % %                 'vxf_resc','vyf_resc','vzf_resc');                                      % final resc
% % % % 

%save(savename,'vr_resc','angle_resc');
% save(savename,'vyf_all','vzf_all');

end