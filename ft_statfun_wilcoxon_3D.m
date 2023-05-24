function cluster = ft_statfun_wilcoxon_3D(cfg, struct1, struct2)

% cfg.statistic = 'dep'; % use the dependent samples T-statistic as a measure to
% 
% cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
%                                % will be used for thresholding
% 
% cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
%                                % required for a selected sample to be included
%                                % in the clustering algorithm (default=0). 
% cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
% cfg.clustertail = 0;
% cfg.alpha = 0.025;               % alpha level of the permutation test
% cfg.numrandomization = 1000;      % number of draws from the permutation distribution
% 
% cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
%                                  % which other sensors it can form clusters
% cfg.channel       = labels;     % cell-array with selected channel labels
% cfg.label = labels;


Nsubj1 = length(struct1.trial);
Nsubj2 = length(struct2.trial);
Nchan = length(cfg.channel);
Ntime = length(struct1.time{1});
Nfreq = length(struct1.freq{1});
time = struct1.time{1};
freq = struct1.freq{1};

%% Make neighbours matrix
neigh_matrix = cell(1, Nchan);

for i = 1 : Nchan
    ch_i = cfg.label(i);
    ix = find( strcmp(ch_i, {cfg.neighbours.label}) == 1);
    temp_neigh = cfg.neighbours(ix).neighblabel;   
    for j = 1 : length(temp_neigh)
        ix = find( strcmp(cfg.label,temp_neigh(j)) == 1);
        neigh_matrix{i} = [neigh_matrix{i} ix];
    end
end
        
%% Make preliminary mask
mask = zeros(Nchan, Ntime, Nfreq);
mask_p = zeros(Nchan, Ntime, Nfreq);
mask_n = zeros(Nchan, Ntime, Nfreq);

prob = zeros(Nchan, Ntime, Nfreq);
stat = zeros(Nchan, Ntime, Nfreq);

for fs = 1: Nfreq
    for ts = 1 : Ntime
        % test on individual channels
        for ch = 1 : Nchan

            dat1 = zeros(1, Nsubj1);
            dat2 = zeros(1, Nsubj2);

            if strcmp(cfg.statistic, 'dep')
                for n = 1 : Nsubj1
                    dat1(n) = struct1.trial{n}(ch,ts,fs);
                    dat2(n) = struct2.trial{n}(ch,ts,fs);
                end
                [p, ~, stats] = signrank(dat1, dat2);
                prob(ch, ts, fs) = p;
                stat(ch, ts, fs) = stats.zval;
            end

            if strcmp(cfg.statistic, 'indep')
                for n = 1 : Nsubj1
                    dat1(n) = struct1.trial{n}(ch,ts,fs);
                end
                for n = 1 : Nsubj2
                    dat2(n) = struct2.trial{n}(ch,ts,fs);
                end            
                [p, ~, stats] = ranksum(dat1, dat2);
                prob(ch, ts, fs) = p;
                stat(ch, ts, fs) = stats.zval;            
            end

            if p <= cfg.alpha && stat(ch, ts, fs) > 0
                mask_p(ch, ts, fs) = 1;
                mask(ch, ts, fs) = 1;
            end

            if p <= cfg.alpha && stat(ch, ts, fs) < 0
                mask_n(ch, ts, fs) = 1;
                mask(ch, ts, fs) = 1;
            end        
        end

        % check min neighb mask p
        ix_m = find(mask_p(:,ts,fs) == 1); 
        for i = 1 : length(ix_m)
            [~, pos] = intersect(neigh_matrix{ix_m(i)}, ix_m);
            if length(pos) < cfg.minnbchan
                mask_p(ix_m(i),ts,fs) = 0;
            end
        end        
        % check min neighb mask n
        ix_m = find(mask_n(:,ts,fs) == 1); 
        for i = 1 : length(ix_m)
            [~, pos] = intersect(neigh_matrix{ix_m(i)}, ix_m);
            if length(pos) < cfg.minnbchan
                mask_n(ix_m(i),ts,fs) = 0;
            end
        end        
    end
end

%% Define mask based on neighbours channels
mask_n2 = local_clusters(cfg, mask_n, neigh_matrix);
mask_p2 = local_clusters(cfg, mask_p, neigh_matrix);
%% Define mask based on timings
mask_n3 = combine_cluster(mask_n2);
mask_p3 = combine_cluster(mask_p2);
%% Define mask based on freq
mask_n4 = combine_cluster_freq(mask_n3);
mask_p4 = combine_cluster_freq(mask_p3);
%% Compute stats
[negclusters, mask_n4]  = cluster_stats(cfg, struct1, struct2, mask_n4, stat);
[posclusters, mask_p4]  = cluster_stats(cfg, struct1, struct2, mask_p4, stat);
%% Outputs
cluster.cfg = cfg;
cluster.time = time;
cluster.freq = freq;
cluster.label = cfg.label;
cluster.dimord = 'chan_time_freq';
cluster.prob = prob;
cluster.stat = stat;
cluster.mask = mask;
cluster.posclusters = posclusters;
cluster.posclusterslabelmat = mask_p4;
cluster.negclusters = negclusters;
cluster.negclusterslabelmat = mask_n4;
end

%% Identifies channel-cluster on individual timestamps
function cluster = define_cluster(ch_s, neigh_matrix, minch)

cluster = cell(1,1);
n_cl = 1;
unused_ch = ch_s';
while ~isempty(unused_ch)
    % check first unclassified channel
    ch_1 = unused_ch(1);
    % construct cluster
    cluster_temp = [ch_1];
    % update unused channels
    unused_ch = setdiff(unused_ch, cluster_temp);
    % neighbours of the unused channel
    neigh_cluster =  neigh_matrix{ch_1};
    % verify if neighbours are part of the cluster
    new_ch = intersect(unused_ch, neigh_cluster);
    % add the neighbours if match
    cluster_temp = [cluster_temp new_ch];
    % update unused channels
    unused_ch = setdiff(unused_ch, cluster_temp);    
    
    % iterate till no new neigh
    while ~isempty(new_ch)
        for i = 1 : length(new_ch)
            ch_2 = new_ch(i);
            neigh_plus = neigh_matrix{ch_2};
            % update neigh cluster
            neigh_cluster = [neigh_cluster neigh_plus];
        end
        neigh_cluster = unique(neigh_cluster);
        % verify if neighbours are part of the cluster
        new_temp = intersect(unused_ch, neigh_cluster);
        % update cluster
        cluster_temp = [cluster_temp new_temp];
        cluster_temp = unique(cluster_temp);   
        % update unused ch
        unused_ch = setdiff(unused_ch, cluster_temp); 
        new_ch = new_temp;
    end
    
    % update cluster
    cluster_temp = unique(cluster_temp);
    
    % verify n chans in cluster
    if length(cluster_temp) >= minch
        cluster{n_cl} = sort(cluster_temp);
        n_cl = n_cl + 1;
    end

end
end

%% Makes a new mask with cluster identified in single timestamps
function new_mask = local_clusters(cfg, mask, neigh_matrix)
% define timestamps with clusters

[Nchan, Ntime, Nfreq] = size(mask);
new_mask = zeros(Nchan,Ntime,Nfreq);
for k = 1 : Nfreq
    ix_t = find( sum(mask(:,:,k),1) > 0);
    for i = 1 : length(ix_t)
        t = ix_t(i);
        ch_s = find(mask(:,t,k) == 1);
        if length(ch_s) >= cfg.minnbchan
            clear cluster
            cluster = define_cluster(ch_s, neigh_matrix, cfg.minnbchan);
            for j = 1 : length(cluster)
                ch_cl = cluster{j};
                new_mask(ch_cl,t,k) = j;
            end
        else
            new_mask(ch_s,t,k) = 0;
        end
    end
end

end

%% Combine cluster on time
function new_mask = combine_cluster(mask)

[Nchan,Ntime,Nfreq] = size(mask);
old_mask = mask;
new_mask = zeros(Nchan, Ntime,Nfreq);

for k = 1 : Nfreq
    n_clusters = 1;
    old_mask_f = old_mask(:,:,k);
    while sum(sum(old_mask_f)) ~= 0
        for t = 1 : Ntime-1 % iteration on the timestamps of the clusters
            data = old_mask_f(:,t);
            cl1 = min(data(data ~= 0));
            if ~isempty(cl1)
                ch_1 = find(data == cl1);
                data2 = old_mask_f(:,t+1);
                cl2 = min(data2(data2 ~= 0));
                if ~isempty(cl2)
                    for cl2 = 1 : max(data2)% iteration on cluster inside next time stamp
                        ch_2 = find(data2 == cl2);
                        C = intersect(ch_1, ch_2);
                        if C >= 1
                            new_mask(ch_1,t,k) = n_clusters;
                            new_mask(ch_2,t+1,k) = n_clusters; 
                            old_mask_f(ch_2, t+1) = n_clusters;
                        end
                    end
                end

            old_mask_f(ch_1,t) = 0;           
            end      
        end
        t = Ntime;
        old_mask_f((find(old_mask_f(:,t) == n_clusters)),t) = 0;
        n_clusters = n_clusters + 1;
    end
end

end

%% Combine cluster on freq
function new_mask_final = combine_cluster_freq(mask)

[Nchan,Ntime,Nfreq] = size(mask);
old_mask = mask;
new_mask = zeros(Nchan,Ntime, Nfreq);
msk_false = false(Nchan, Ntime);
msk_false2 = false(Nchan, Ntime, Nfreq);

n_clusters = 1;
while sum(sum(sum(old_mask))) ~= 0
    for fs = 1 : Nfreq-1 % iteration on the timestamps of the clusters
        mask_f = old_mask(:,:,fs);
        cl1 = min(min(mask_f(mask_f ~= 0)));
        if ~isempty(cl1)
            data = [];
            for t = 1 : Ntime
               ch_1 = find(mask_f(:,t) == cl1);
               if ~isempty(ch_1)
                   clear new_data
                   new_data = [reshape(ch_1,length(ch_1),1), repmat(t, length(ch_1), 1)];
                   data = [data; new_data];
               end
            end
            mask_f2 = old_mask(:,:,fs+1);
            cl2 = min(min(mask_f2(mask_f2 ~= 0)));
            mn = cl2;
            if ~isempty(cl2)
                for cl2 = mn : max(max(mask_f2))% iteration on cluster inside next time stamp
                    data2 = [];
                    for t = 1 : Ntime
                       ch_2 = find(mask_f2(:,t) == cl2);
                       if ~isempty(ch_2)
                           clear new_data
                           new_data = [reshape(ch_2,length(ch_2),1), repmat(t, length(ch_2), 1)];
                           data2 = [data2; new_data];
                       end
                    end
                    
                    Int = intersect(data, data2, 'rows');
                    C = size(Int);
                    if C(1) >= 1
                        msk=msk_false; msk((old_mask(:,:,fs)) == cl1) = true; 
                        msk2=msk_false2; msk2(:,:,fs) = msk;
                        new_mask(msk2) = n_clusters; 
                        msk=msk_false; msk((old_mask(:,:,fs+1)) == cl2) = true;
                        msk2=msk_false2; msk2(:,:,fs+1) = msk;
                        new_mask(msk2) = n_clusters; 
                        msk=msk_false; msk((old_mask(:,:,fs+1)) == cl2) = true;
                        msk2=msk_false2; msk2(:,:,fs+1) = msk;
                        old_mask(msk2) = n_clusters;
                    end
                    
                end
            else
                msk=msk_false; msk((old_mask(:,:,fs)) == cl1) = true;
                msk2=msk_false2; msk2(:,:,fs) = msk;
                new_mask(msk2) = n_clusters;
                n_clusters = n_clusters + 1;
                
            end
        
        msk=msk_false; msk(old_mask(:,:,fs) == cl1) = true;
        msk2=msk_false2; msk2(:,:,fs) = msk;
        old_mask(msk2) = 0;           
        end      
    end
    fs = Nfreq;
    msk=msk_false; msk((old_mask(:,:,fs)) == n_clusters) = true;
    msk2=msk_false2; msk2(:,:,fs) = msk;
    old_mask(msk2) = 0;
    n_clusters = n_clusters + 1;
    
    %take care of isolated cluster in the last freq stamp
    mask_f = old_mask(:,:,fs);
    cln = min(min(mask_f(mask_f ~= 0)));
    if ~isempty(cln)
        msk=msk_false; msk((old_mask(:,:,fs)) == cln) = true;
        msk2=msk_false2; msk2(:,:,fs) = msk;
        old_mask(msk2) = 0;
        new_mask(msk2) = n_clusters;
        n_clusters = n_clusters + 1;   
    end
end
new_mask_final = new_mask;
n_cluster_final = 0;
n_clusters = max(max(max(new_mask)));
for i = 1 : n_clusters
    K = find(new_mask == i);
    if ~isempty(K)
        n_cluster_final = n_cluster_final +1;
        new_mask_final(new_mask == i) = n_cluster_final;
    end
end



end

%% Compute cluster stats
function [clusters, new_mask]  = cluster_stats(cfg, struct1, struct2, mask, stat)

[Nchan,Ntime,Nfreq] = size(mask);
msk_false = false(Nchan, Ntime);
msk_false2 = false(Nchan, Ntime, Nfreq);

n_cluster = max(max(max(mask)));  

if n_cluster == 0
    clusters = [];
    new_mask = mask;
else

    clusters = struct;

    Nsubj1 = length(struct1.trial);
    Nsubj2 = length(struct2.trial);
    dat1 = zeros(1, Nsubj1);
    dat2 = zeros(1, Nsubj2);

    for n_clu = 1 : n_cluster
        mask_temp = msk_false2; 
        for fs = 1 : Nfreq
            msk=msk_false; msk(mask(:,:,fs) == n_clu) = true;
            mask_temp(:,:,fs) = msk;
        end
        
        % get length cluster
        t_cluster = [];
        for t = 1 : Ntime
            M = reshape(squeeze(mask_temp(:,t,:)), Nfreq*Nchan, 1);
            M2 = find(M == true);
            if length(M2) >= 1
                t_cluster = [t_cluster t];
            end
        end
        t_range = [min(t_cluster) max(t_cluster)];
        
        if (t_range(2)-t_range(1)) >= cfg.mint
        
            Npoint = sum(sum(sum(mask_temp)));
            if Npoint ~= 0
                for n = 1 : Nsubj1
                    dat1(n) = sum(sum(sum(struct1.trial{n}(mask_temp)))) / Npoint;
                end   
                for n = 1 : Nsubj2
                    dat2(n) = sum(sum(sum(struct2.trial{n}(mask_temp)))) / Npoint;
                end
            else
                dat1 = zeros(1, Nsubj1);
                dat2 = zeros(1, Nsubj2);
            end

            [p, ~, stats] = signrank(dat1, dat2);
            p_real = p;
            if isfield(stats, 'zval')
                stat_real = stats.zval;
            else
                stat_real = 0;
            end

            stat_perm =  zeros(1, cfg.numrandomization);
            p_perm =  zeros(1, cfg.numrandomization);
            for n_perm = 1 : cfg.numrandomization
                dat_all = [dat1 dat2];
                dat_all = dat_all(randperm(length(dat_all)));
                dat_p1 = dat_all(1:Nsubj1);
                dat_p2 = dat_all(Nsubj1+1:Nsubj1+Nsubj2);

                if strcmp(cfg.statistic, 'dep')
                    [p, ~, stats] = signrank(dat_p1, dat_p2);
                    p_perm(n_perm) = p;

                    if isfield(stats, 'zval')
                        stat_perm(n_perm) = stats.zval;
                    else
                        stat_perm(n_perm) = 0;
                    end
                end
                if strcmp(cfg.statistic, 'indep')
                    [p, ~, stats] = ranksum(dat_p1, dat_p2);
                    p_perm(n_perm) = p;

                    if isfield(stats, 'zval')
                        stat_perm(n_perm) = stats.zval;
                    else
                        stat_perm(n_perm) = 0;
                    end
                end        

            end
            p_montecarlo = length(find(p_perm < p_real))/ cfg.numrandomization;
    %         p_montecarlo = length(find(abs(stat_perm) > abs(stat_real)))/ cfg.numrandomization;
            clusters.p_mc(n_clu) = p_montecarlo;
%             clusters.stat_mc(n_clu) = stat_real;
            clusters.stat_mc(n_clu) = max(max(max( abs(stat.*mask_temp)  )));
            clusters.p(n_clu) = p_real;
        else
            clusters.p_mc(n_clu) = 1;
            clusters.stat_mc(n_clu) = 0;
            clusters.p(n_clu) = 1;            
        end

    end

    % sort clusters
    if isfield(clusters, 'p_mc')
        [~ , ord] = sort(clusters.p, 'ascend');
        clusters.p_mc = clusters.p_mc(ord);
        clusters.stat_mc = clusters.stat_mc(ord);
        clusters.p = clusters.p(ord);
    end
    new_mask = mask;
    for i = 1 : n_clu
        new_mask(mask == ord(i)) = i;
    end
end

end
