function [trind_mat tr_mat tt_mat rescind dirind] = fct_get_return1Traj1D_v2(tgrid,texit_ind,r0z_sample,v0z_sample,R,rcov,NrMax,PlotOption)

% Nr    - number of return that is interesting
% NrMax - highest number of return 
% tgrid - time axis
% R     - r vector
% number of return to get
% rcov  - atomis radius, i.e. condition for rescattering

% define the output
trind_mat = zeros(length(texit_ind),NrMax);

for k=1:length(texit_ind)
    ts_indk  = texit_ind(k);
    rk       = R(ts_indk)'+v0z_sample(k).*(tgrid-tgrid(ts_indk))+r0z_sample(k);
    tgridk   = ts_indk:1:length(tgrid);
    rkoi     = rk(tgridk);   
    d_oi     = -1.*abs(rkoi);
    if length(d_oi)>3
        [pks_,locs_]  = findpeaks(d_oi,'NPEAKS',NrMax);
    end
    
%     hold off
%     plot(tgridk,d_oi,'r-'); hold on
%     plot(locs_,pks_,'b.'); hold on
%     drawnow 
    if length(locs_)~=0
       for l=1:length(locs_)
        %         if rkoi(locs_(l)-1).*rkoi(locs_(l)+1)>0
        %             locs_(l) = NaN;
        %         end
        
        if (abs(d_oi(locs_(l)))>rcov) && (rkoi(locs_(l)-1).*rkoi(locs_(l)+1)>0)
            locs_(l) = NaN;
        end
        
       end
    end
    trind_mat(k,1:1:length(locs_))  = locs_+ts_indk.*ones(1,length(locs_));
end

%% post processing
ind = find(isnan(trind_mat)==1);
trind_mat(ind) = 0;

%% resc time
tr_mat = zeros(length(texit_ind),NrMax);
tt_mat = zeros(length(texit_ind),NrMax);
for k=1:NrMax
    trind_mat_k  = trind_mat(:,k);                              % index of starting
    indk         = find(trind_mat_k~=0);                        % find non-zeros, i.e. the actual returns
    trind_mat_k_ = trind_mat(indk,k);                           % get return index of the actual returns
    texit_ind_k = texit_ind(indk);                              % get starting index of the ones with return
    tr_mat(indk,k) = tgrid(trind_mat_k_);
    tt_mat(indk,k) = tgrid(trind_mat_k_)-tgrid(texit_ind_k);
end

%% finally get trajectory type
trind_mat_ = trind_mat(:,1);
% identify direct trajectories
ind    = find(trind_mat_~=0);
rescind = ind;

ind    = find(trind_mat_==0);
dirind = ind;


%% plot
cc = lines((NrMax));
if PlotOption==1
    for m=1:NrMax
        figure;

        indm = find(tr_mat(:,m)>0);
        subplot(2,1,1)
        plot(tgrid(texit_ind(indm)),tr_mat(indm,m),'.','Color',cc(m,:).*0.7); hold on
%         legend('1st','2nd','3rd','4th','5th')#
        grid on
        subplot(2,1,2)
        plot(tgrid(texit_ind(indm)),tt_mat(indm,m),'.','Color',cc(m,:).*0.7); hold on
        grid on
    end
end 

end