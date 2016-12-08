close all
% clear all

c0       = 299792458./(2.1876912633e6);        % au, speed of light 
wvlm     = 1800e-9;                            % 
wvl      = wvlm/(5.2917720859e-11);                            % [au]
omega    = (2*pi*c0)/wvl;                                      % [au], angular frequency
T        = 2*pi/omega; 
IWcm     = 1e14;   
Up       = (0.09337.*IWcm.*wvlm.^2)./27.211;   % au  
    
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


% for i=1:length(vpagrid)/2
%     for j=1:length(vpegrid)
%         ang(i,j) = atan(vpegrid(j)./vpagrid(i));
%     end
% end

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

nsample = 1e7;

    % make the 2d histogram to a list 
P_r =  singlet_r(:);
P_r = P_r./sum(P_r);
    
Eaxis = 0:0.5:350;       % to 11 UP
tot_spec_r = 0;
    
for ii = 1: 1
    tic
    ind  = fct_gen_distr(P_r',1,nsample);

    [Vpegrid,Vpagrid] = meshgrid(vpegrid,vpagrid);
    %hist_v0xv0y = hist3([vxf vyf], {vpegrid' vpegrid'});

    xsample = Vpagrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1));
    ysample = Vpegrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1));

    % energy of intersting particles

    Esample = 0.5.*(xsample.^2+ysample.^2) * 27.211;
    
    histEaxis = hist(Esample,Eaxis);
    tot_spec_r = tot_spec_r + histEaxis;
    
   toc
end

%ind       = find(vpagrid>0);
%Eaxis_ana = 27.211*0.5*vpagrid(ind).^2
tot_spec_r_ana = sum(singlet_r')./abs((vpagrid));


subplot(2,2,2)
semilogy(Eaxis,tot_spec_r./max(tot_spec_r)); hold on
semilogy(27.211*0.5*vpagrid.^2,tot_spec_r_ana./max(tot_spec_r_ana),'r.');
%plot(vpagrid,sum(singlet_l'),'.')
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

P_l =  singlet_l(:);
P_l = P_l./sum(P_l);

tot_spec_l = 0;
    
for ii = 1: 1
    tic
    ind  = fct_gen_distr(P_l',1,nsample);

    [Vpegrid,Vpagrid] = meshgrid(vpegrid,vpagrid);
    %hist_v0xv0y = hist3([vxf vyf], {vpegrid' vpegrid'});

    xsample = Vpagrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1));
    ysample = Vpegrid(ind)+rand(1,nsample).*(vpegrid(2)-vpegrid(1));

    % energy of intersting particles

    Esample = 0.5.*(xsample.^2+ysample.^2) * 27.211;
    
    histEaxis = hist(Esample,Eaxis);
    tot_spec_l = tot_spec_l + histEaxis;
    
   toc
end

subplot(2,2,4)
semilogy(Eaxis,tot_spec_l)

%% 
% figure (2)
% PAP = (tot_spec_r-tot_spec_l)/(tot_spec_l+tot_spec_r);
% plot(PAP)
% 


 
