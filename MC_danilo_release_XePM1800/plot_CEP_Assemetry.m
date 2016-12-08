close all 
clear all


wvlm            = 1800e-9; 
I0              = 8e13;   
nr_ceps         = 50;
cutoff_angle    = 5; % opening of the cone, in DEG
save_ati_switch = 0;

pulse_length    = 13.1; 
folder_name     = 'n_3_N_0_wvl_1800_I0_8.0e+13_dt_0.1_kmax_300_nsample_9e+05';
% CEP_ati_spec_l;
% CEP_ati_spec_r;

PlotOpt = 0;
PlotOpt2= 1;

gate1   = (40: 10: 120);
gate2   = (80: 10: 150);
gate3   = (90: 10: 160);
gate4   = (150: 10: 200);
a       =  1;
b       =  0;

% for ii = 1 : length(gate1)
%     for kk = 1 : length(gate2)
%         for jj = 1 : length(gate3)
%             for ll = 1 : length(gate4)
%                 
%                 gate_com(a,:) = [gate1(ii) gate2(jj) gate3(kk) gate4(ll)];
%                 a = a+1;
%             end
%         end
%     end
% end

tic

Ebins   = [115 120 165 170];


% for iii = 1 : a-1
    
%     if gate_com(iii,1) < gate_com(iii,2) && gate_com(iii,3)< gate_com(iii,4)   
    
%     b            = b+1;
%     Ebins(b,:)   = gate_com(b,:);
     
    
    for k = 1:nr_ceps   
    % load file
    file_name = strcat(folder_name,'\',num2str(k),'_',num2str(nr_ceps),'ceps_','ATI.mat');
    load(file_name);
    
    % divide into left and right on positive e axis 
    ind = find(min(E_axis_ati)==E_axis_ati); 
    Eaxis_r    = fliplr(E_axis_ati(1:ind));
    ati_spec_r = fliplr(ati_spec_r(1:ind));  
    Eaxis_l    = E_axis_ati(ind:length(E_axis_ati));
    ati_spec_l = ati_spec_l(ind:length(E_axis_ati));
    
    if PlotOpt==1
       semilogy(Eaxis_l,ati_spec_l,'.'); hold on
       semilogy(Eaxis_r,ati_spec_r,'.'); hold on
       drawnow
    end
    %force them onto the same axis
    Estep      = 0.25;
    Eaxis      = 0:Estep:(max([max(Eaxis_r) max(Eaxis_l)]));
    ati_spec_r = interp1(Eaxis_r,ati_spec_r,Eaxis);
    ati_spec_l = interp1(Eaxis_l,ati_spec_l,Eaxis);
    
    if PlotOpt==1
        semilogy(Eaxis,ati_spec_l,'-'); hold on
        semilogy(Eaxis,ati_spec_r,'-'); hold on
        drawnow
    end

    
    % calculate asymmetry and set left right cep matrix spectra
    CEP_ati_spec_l(k,:)  = ati_spec_l; 
    CEP_ati_spec_r(k,:)  = ati_spec_r;
    CEP_ati_asy_spec(k,:)= (ati_spec_r - ati_spec_l)./(ati_spec_l + ati_spec_r);
    
    % cep axis
    CEPaxis(k) = k * 2 * pi / nr_ceps;
    CEPaxis(k) = fliplr(CEPaxis(k));
    
    % calculate potato
%     ind             = isnan(ati_spec_l);
%     ati_spec_l(ind) = 0;
%     ind             = isnan(ati_spec_r);
%     ati_spec_r(ind) = 0;
    
%     Alow_l   = sum( ati_spec_l( (b,1)): ati_spec_l(Ebins(b,2)) );
%     Alow_r   = sum( ati_spec_r(Ebins(b,1)): ati_spec_r(Ebins(b,2)) );
%     Alow(k)  = (Alow_r-Alow_l)/(Alow_r+Alow_l);
%     
%     Ahigh_l   = sum( ati_spec_l(Ebins(b,3)): ati_spec_l(Ebins(b,4)) );
%     Ahigh_r   = sum( ati_spec_r(Ebins(b,3)): ati_spec_r(Ebins(b,4)) );
%     Ahigh(k)  = (Ahigh_r-Ahigh_l)/(Ahigh_r+Ahigh_l);
    
    
       ind      = find(Eaxis>Ebins(1) & Eaxis<Ebins(2));
      Alow_l   = sum(ati_spec_l(ind));
    Alow_r   = sum(ati_spec_r(ind));
    Alow(k)  = (Alow_r-Alow_l)/(Alow_r+Alow_l);
    
    ind      = find(Eaxis>Ebins(3) & Eaxis<Ebins(4));
    Ahigh_l  = sum(ati_spec_l(ind));
    Ahigh_r  = sum(ati_spec_r(ind));
    Ahigh(k) = (Ahigh_r-Ahigh_l)/(Ahigh_r+Ahigh_l);   
    
%     end


    [theta,rho] = cart2pol(Alow,Ahigh);
    R = mean(rho);
    
%     end
    
end
  
R_max  = max(R);

toc

    subplot(2,3,1)
    pcolor(CEPaxis,Eaxis,log10(CEP_ati_spec_l'))
    colorbar
    shading flat

    title('ATI left')
    subplot(2,3,2)
    pcolor(CEPaxis,Eaxis,log10(CEP_ati_spec_r'))
    colorbar
    shading flat
    title('ATI right')

    subplot(2,3,3)
    pcolor(CEPaxis,Eaxis,CEP_ati_asy_spec')
    colorbar
    shading flat
    title('Asymmetry')

    subplot(2,3,4)
    plot(Alow,Ahigh,'.-')
    title('Potato')
    grid on

    subplot(2,3,5)
    plot(CEPaxis,Alow,'r.-'); hold on
    plot(CEPaxis,Ahigh,'k.-')
    grid on

    subplot(2,3,6)
    plot(theta,rho,'r.'); hold on
    grid on
        
  figure (2)
      pcolor(CEPaxis,Eaxis,CEP_ati_asy_spec')
    colorbar
    shading flat
    title('Asymmetry')
    
