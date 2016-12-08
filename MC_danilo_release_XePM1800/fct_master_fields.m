function [ tgrid, A, E_real, Env ] = fct_master_fields_shared(type, start_at_0, wvlnm,...
                                                            IWcm, CEP, k2prime, fwhm_ftl,...
                                                            dt, cutoff, z, PlotOpt)
%FCT_MASTER_FIELDS by Danilo Zille
% Input:
%
% PlotOpt ... give out plots ? 1==yes
% type....specifiy which field you want. options: "gauss", "sin2"
% start_at_0 ... do you want the field to start at t = 0 ? 1==yes
% wvlnm ... wavelength in nm
% IWcm ... intensity in W/cm^2
% CEP...carrier-envelope-phase of pulse
% dt ... resolution of time grid
% cutoff ... only for gauss pulse. at what point do you want to truncate? cutoff = 0.1 means thats the intensity envelope decayed down to 0.1*max(Envelope). cutoff=0.0000001 seems like a safe choice (check with plot)
% z ... propagation distance in dispersive medium. only for "gauss". unit is mm
% kdoubleprime ... the group velocity dispersion in fs^2/mm
% fwhm_nochirp ... the pulse length defined as FWHM of the intensity envelope (transform limited). in fs. only for gauss
% n ... number of cycles for "sin2" pulse.

c0           = 299792458./(2.1876912633e6);        % au, speed of light
E0           =  sqrt(IWcm/3.509e16); %convert to au
wvlm         = wvlnm*1e-9;                   
wvl          = wvlm/(5.2917720859e-11);         % [au]
omega        = (2*pi*c0)/wvl;                           % [au], angular frequency
k2prime      = k2prime*(41.34)^2;% 44.929  *41^2 ;  %a.u. GVD  au^2/mm
fwhm_ftl     = fwhm_ftl*41.34; % to a.u.

if (strcmp(type,'gauss'))
    %define the vector potential in analytical form. then derive E-field
    % field from Effects of the carrier-envelope phase of chirped laser pulses in the multiphoton ionization regime - E.Cormier
    %note that the field was slightly altered, since in the paper the pulse lenght is measured in the field envelope. but we measure intensity envelope
    xi =  4*log(2)*k2prime*z / fwhm_ftl^2;
    t_end = sqrt(log(cutoff)/(-4*log(2)) * (1+xi^2) )*fwhm_ftl;
               
    if (start_at_0 == 1)
        tgrid  = 0:dt:2*t_end;
        E = @(t) E0 * exp(-1i*(omega*(t-t_end)+CEP) - 2*log(2).*((t-t_end)/fwhm_ftl).^2*(1/(1-1i*xi)));
    else
        tgrid  = -t_end:dt:t_end;
        E = @(t) E0 * exp(-1i*(omega*t+CEP) - 2*log(2).*(t/fwhm_ftl).^2*(1/(1-1i*xi)));
    end
         
    E       = E(tgrid);
    E_real  = real((E));
    A       = -cumsum(E_real)*dt;
    Env     = @(t) sqrt(conj(E).*E);
    Env     = Env(tgrid);
       
    if (PlotOpt==1)
        plot(tgrid,E_real,tgrid,A*omega,tgrid,Env);
        legend('E','A*omega','Envelope');
    end
end
end

 