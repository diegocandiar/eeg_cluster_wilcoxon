function cluster = ft_statfun_cluster(cfg, struct1, struct2)
% This function performs cluster-based permutation analysis in
% two-dimensional data. Data must be organized in struct format with
% subfields trial and time, organized in cell arrays per sample.
%
% IMPORTANT NOTE: This algorithm does not make multiple comparison 
% correction based on the number of samples compared, i.e., the 
% results are not affected by the number of channels in your EEG/MEG
% system, nor the latency length in which the comparisons are done.
% E.g., if your data has a true effect in 0-5 latency, the 
% results should not differ if the comparison is done between 
% 0-5 or 0-10 seconds.
%
% The algorith outputs are:
% cluster.prob contains p-values from univariate tests
% cluster.stat contains statistic from univariate tests
% cluster.mask cluster definition matrix
% cluster.posclusters cluster statistical details for postive clusters
% cluster.negclusters cluster statistical details for negative clusters
%
% cfg struct must contain:
% cfg.statistic = 'dep_param' 'dep_non_param' 'indep_param' or
% 'indep_nonparam': % define if dependent or independent test, and if
% parametrical or non-parametrical statistical test
% 
% cfg.clusteralpha: alpha level of the sample-specific test
% statistic that will be used for thresholding
% 
% cfg.minnbchan: minimum number of neighborhood channels that is
% required for a selected sample to be included in the clustering algorithm
%
% cfg.alpha: alpha level of the permutation test
%
% cfg.numrandomization: number of draws from the permutation distribution
% 
% cfg.neighbours: the neighbours definition 
%
% cfg.label: cell-array with selected channel labels 

%%

Nsubj1 = length(struct1.trial);
Nsubj2 = length(struct2.trial);
Nchan = length(cfg.label);
Ntime = length(struct1.time{1});
time = struct1.time{1};

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
mask = zeros(Nchan, Ntime);
mask_p = zeros(Nchan, Ntime);
mask_n = zeros(Nchan, Ntime);

prob = zeros(Nchan, Ntime);
stat = zeros(Nchan, Ntime);

for ts = 1 : Ntime
    % test on individual channels
    for ch = 1 : Nchan

        dat1 = zeros(1, Nsubj1);
        dat2 = zeros(1, Nsubj2);
        
        if strcmp(cfg.statistic, 'dep_param')
            for n = 1 : Nsubj1
                dat1(n) = struct1.trial{n}(ch,ts);
                dat2(n) = struct2.trial{n}(ch,ts);
            end
            [~,p, ~, stats] = ttest(dat1, dat2);
            prob(ch, ts) = p;
            stat(ch, ts) = stats.tstat;
        end
        
        if strcmp(cfg.statistic, 'dep_non_param')
            for n = 1 : Nsubj1
                dat1(n) = struct1.trial{n}(ch,ts);
                dat2(n) = struct2.trial{n}(ch,ts);
            end
            [p, ~, stats] = signrank(dat1, dat2);
            prob(ch, ts) = p;
            stat(ch, ts) = stats.zval;
        end
        if strcmp(cfg.statistic, 'indep_param')
            for n = 1 : Nsubj1
                dat1(n) = struct1.trial{n}(ch,ts);
            end
            for n = 1 : Nsubj2
                dat2(n) = struct2.trial{n}(ch,ts);
            end            
            [~,p, ~, stats] = ttest2(dat1, dat2);
            prob(ch, ts) = p;
            stat(ch, ts) = stats.tstat;            
        end
        if strcmp(cfg.statistic, 'indep_non_param')
            for n = 1 : Nsubj1
                dat1(n) = struct1.trial{n}(ch,ts);
            end
            for n = 1 : Nsubj2
                dat2(n) = struct2.trial{n}(ch,ts);
            end            
            [p, ~, stats] = ranksum(dat1, dat2);
            prob(ch, ts) = p;
            stat(ch, ts) = stats.zval;            
        end        
        if p <= cfg.alpha && stat(ch, ts) > 0
            mask_p(ch, ts) = 1;
            mask(ch, ts) = 1;
        end
        
        if p <= cfg.alpha && stat(ch, ts) < 0
            mask_n(ch, ts) = 1;
            mask(ch, ts) = 1;
        end        
    end
    
    % check min neighb mask p
    ix_m = find(mask_p(:,ts) == 1); 
    for i = 1 : length(ix_m)
        [~, pos] = intersect(neigh_matrix{ix_m(i)}, ix_m);
        if length(pos) < cfg.minnbchan
            mask_p(ix_m(i),ts) = 0;
        end
    end        
    % check min neighb mask n
    ix_m = find(mask_n(:,ts) == 1); 
    for i = 1 : length(ix_m)
        [~, pos] = intersect(neigh_matrix{ix_m(i)}, ix_m);
        if length(pos) < cfg.minnbchan
            mask_n(ix_m(i),ts) = 0;
        end
    end        
end

%% Define mask based on neighbours channels
mask_n2 = local_clusters(cfg, mask_n, neigh_matrix);
mask_p2 = local_clusters(cfg, mask_p, neigh_matrix);
%% Define mask based on timings
mask_n3 = combine_cluster(mask_n2);
mask_p3 = combine_cluster(mask_p2);
%% Compute stats
[negclusters, mask_n3]  = cluster_stats(cfg, struct1, struct2, mask_n3, stat);
[posclusters, mask_p3]  = cluster_stats(cfg, struct1, struct2, mask_p3, stat);
%% Outputs
cluster.cfg = cfg;
cluster.time = time;
cluster.label = cfg.label;
cluster.dimord = 'chan_time';
cluster.prob = prob;
cluster.stat = stat;
cluster.mask = mask;
cluster.posclusters = posclusters;
cluster.posclusterslabelmat = mask_p3;
cluster.negclusters = negclusters;
cluster.negclusterslabelmat = mask_n3;
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
ix_t = find( sum(mask,1) > 0);
[Nchan, Ntime] = size(mask);
new_mask = zeros(Nchan,Ntime);

for i = 1 : length(ix_t)
    t = ix_t(i);
    ch_s = find(mask(:,t) == 1);
    if length(ch_s) >= cfg.minnbchan
        clear cluster
        cluster = define_cluster(ch_s, neigh_matrix, cfg.minnbchan);
        for j = 1 : length(cluster)
            ch_cl = cluster{j};
            new_mask(ch_cl,t) = j;
        end
    else
        new_mask(ch_s,t) = 0;
    end
end

end

%% Combine cluster on time
function new_mask = combine_cluster(mask)

[Nchan, Ntime] = size(mask);
old_mask = mask;
new_mask = zeros(Nchan,Ntime);
n_clusters = 1;

while sum(sum(old_mask)) ~= 0
    for t = 1 : Ntime-1 % iteration on the timestamps of the clusters
        data = old_mask(:,t);
        cl1 = min(data(data ~= 0));
        if ~isempty(cl1)
            ch_1 = find(data == cl1);
            data2 = old_mask(:,t+1);
            cl2 = min(data2(data2 ~= 0));
            if ~isempty(cl2)
                for cl2 = 1 : max(data2)% iteration on cluster inside next time stamp
                    ch_2 = find(data2 == cl2);
                    C = intersect(ch_1, ch_2);
                    if length(C) >= 1
                        new_mask(ch_1,t) = n_clusters;
                        new_mask(ch_2,t+1) = n_clusters; 
                        old_mask(ch_2,t+1) = n_clusters;
                    end
                end
            else
                n_clusters = n_clusters + 1; % COMMENT IF INF LOOP
            end

        old_mask(ch_1,t) = 0;           
        end      
    end
    t = Ntime;
    old_mask((find(old_mask(:,t) == n_clusters)),t) = 0;
    n_clusters = n_clusters + 1;
end


end

%% Compute cluster stats
function [clusters, new_mask2]  = cluster_stats(cfg, struct1, struct2, mask, stat)

[Nchan,Ntime] = size(mask);
tag_clusters = unique(mask(mask ~= 0));

n_cluster = length(tag_clusters);  

new_mask = zeros(Nchan,Ntime);
new_mask2 = zeros(Nchan,Ntime);

if n_cluster == 0
    clusters = [];
    new_mask2 = mask;
else

    clusters = struct;

    Nsubj1 = length(struct1.trial);
    Nsubj2 = length(struct2.trial);
    dat1 = zeros(1, Nsubj1);
    dat2 = zeros(1, Nsubj2);

    new_n_clusters = 0;
    for n_c = 1 : n_cluster
        n_clu = tag_clusters(n_c);
        mask_temp = zeros(Nchan,Ntime);
        mask_temp(mask == n_clu) = 1;
        
        % get length cluster
        t_cluster = [];
        for t = 1 : Ntime
            M = reshape(squeeze(mask_temp(:,t)), Nchan, 1);
            M2 = find(M == true);
            if length(M2) >= 1
                t_cluster = [t_cluster t];
            end
        end
        t_range = [min(t_cluster) max(t_cluster)];
        
        if (t_range(2)-t_range(1)) >= cfg.mint
        
            Npoint = sum(sum(mask_temp));
            for n = 1 : Nsubj1
                dat1(n) = sum(sum(struct1.trial{n}.*(mask_temp))) / Npoint;
            end   
            for n = 1 : Nsubj2
                dat2(n) = sum(sum(struct2.trial{n}.*(mask_temp))) / Npoint;
            end       
            if strcmp(cfg.statistic, 'dep_param')
                [~,p, ~, stats] = ttest(dat1, dat2);
                p_real = p;
                stat_real = stats.tstat;
            end
            if strcmp(cfg.statistic, 'dep_non_param')
                [p, ~, stats] = signrank(dat1, dat2);
                p_real = p;
                stat_real = stats.zval;
            end            
            
            if strcmp(cfg.statistic, 'indep_param')
                [~,p, ~, stats] = ttest2(dat1, dat2);
                p_real = p;
                stat_real = stats.tstat;
            end
            if strcmp(cfg.statistic, 'indep_non_param')
                [p, ~, stats] = ranksum(dat1, dat2);
                p_real = p;
                stat_real = stats.zval;
            end
            
            stat_perm =  zeros(1, cfg.numrandomization);
            p_perm =  zeros(1, cfg.numrandomization);
            for n_perm = 1 : cfg.numrandomization
                dat_all = [dat1 dat2];
                dat_all = dat_all(randperm(length(dat_all)));
                dat_p1 = dat_all(1:Nsubj1);
                dat_p2 = dat_all(Nsubj1+1:Nsubj1+Nsubj2);

                if strcmp(cfg.statistic, 'dep_non_param')
                    [p, ~, stats] = signrank(dat_p1, dat_p2);
                    p_perm(n_perm) = p;
                    stat_perm(n_perm) = stats.zval;
                end
                if strcmp(cfg.statistic, 'indep_non_param')
                    [p, ~, stats] = ranksum(dat_p1, dat_p2);
                    p_perm(n_perm) = p;
                    stat_perm(n_perm) = stats.zval;
                end        
                if strcmp(cfg.statistic, 'indep_param')
                    [~,p, ~, stats] = ttest2(dat_p1, dat_p2);
                    p_perm(n_perm) = p;
                    stat_perm(n_perm) = stats.tstat;
                end
                if strcmp(cfg.statistic, 'dep_param')
                    [~,p, ~, stats] = ttest(dat_p1, dat_p2);
                    p_perm(n_perm) = p;
                    stat_perm(n_perm) = stats.tstat;
                end
            end
            
            new_n_clusters = new_n_clusters + 1; 
            
            p_montecarlo = length(find(p_perm < p_real))/ cfg.numrandomization;
    %         p_montecarlo = length(find(abs(stat_perm) > abs(stat_real)))/ cfg.numrandomization;
            clusters.p_mc(new_n_clusters) = p_montecarlo;
            clusters.p_perm{new_n_clusters} = p_perm;
            clusters.stat_mc(new_n_clusters) = max(max( abs(stat.*mask_temp)  ));
            clusters.p(new_n_clusters) = p_real;
            
            new_mask = new_mask + new_n_clusters.*mask_temp;
%             new_mask(mask_temp) = new_n_clusters ;
            
        else
%             clusters.p_mc(n_clu) = 1;
%             clusters.stat_mc(n_clu) = 0;
%             clusters.p(n_clu) = 1;     
%             new_mask(mask_temp) = 0;
        end
    end

    % sort clusters
    if isfield(clusters, 'p_mc')
        [~ , ord] = sort(clusters.stat_mc, 'descend');
        clusters.p_mc = clusters.p_mc(ord);
        clusters.stat_mc = clusters.stat_mc(ord);
        clusters.p = clusters.p(ord);
        clusters.p_perm = clusters.p_perm(ord);
    end
    
    for i = 1 : new_n_clusters 
        new_mask2(new_mask == ord(i)) = i;
    end
end

end
