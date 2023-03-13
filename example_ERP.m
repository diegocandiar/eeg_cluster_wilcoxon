%%
clc 
clearvars
close all 

fieldtrip_path = '...';
addpath(fieldtrip_path);

load('data_erp_example.mat')
load('gGamma32_neighb.mat')

cfg = [];
cfg.mint = 0.05;
cfg.statistic = 'dep_non_param'; % 'dep_param' or 'dep_non_param'
cfg.alpha = 0.01/1;
cfg.minnbchan = 4; 
cfg.numrandomization = 10000;  
cfg.neighbours    = neighbours;
cfg.label = {'FP1','FP2','AFz','F7','F3','F4','F8','FC5','FC1',...
    'FC2','FC6','T7','C3','Cz','C4','T8','CP5','CP1','CP2',...
    'CP6','P7','P3','Pz','P4','P8','PO7','O1','Oz','O2','PO8','PO9','PO10'};

stat = ft_statfun_cluster(cfg, ERP_noflash, ERP_flash) ; 

if ~isempty(stat.posclusters)    
    pos = stat.posclusterslabelmat == 1; 
else
    pos = false(32,length(stat.time));
end

if ~isempty(stat.negclusters)     
    neg = stat.negclusterslabelmat == 1;
else
    neg = false(32,length(stat.time));
end

%% plot cluster stats
cfg  = [];
cfg.operation = 'x1-x2';
cfg.parameter = 'trial';
BHIdiff = ft_math(cfg,ERP_flash, ERP_noflash);
BHIdiff.label = ERP_flash.label;

cfg = [];
% cfg.parameter = 'trial';
cfg.removemean = 'no';
raweffect = ft_timelockanalysis(cfg, BHIdiff);
trial_tl = ft_timelockanalysis(cfg, ERP_flash);
trial_tl2 = ft_timelockanalysis(cfg, ERP_noflash);
        
%find time range neg
mask1 = stat.negclusterslabelmat;
t_cluster = [];
for t = 1 : length(stat.time)
    M = mask1(:,t,:);
    mx = M(M~= 0);
    mx = min(mx);
    if mx == 1
        t_cluster = [t_cluster t];
    end
end
t_range = [min(t_cluster) max(t_cluster)];
if isempty(t_range)
    t_range = [1 1];
end 

%find time range pos
mask2 = stat.posclusterslabelmat;
t_cluster = [];
for t = 1 : length(stat.time)
    M = mask2(:,t,:);
    mx = M(M~= 0);
    mx = min(mx);
    if mx == 1
        t_cluster = [t_cluster t];
    end
end
t_range = [min(t_cluster) max(t_cluster)];
if isempty(t_range)
    t_range = [1 1];
end


%% VISUALIZE v2
raweffect.med = -stat.stat;
mask = stat.posclusterslabelmat;
mask2 = stat.negclusterslabelmat;
tw = round(512*0.1);
Nt = floor(length(stat.time)/tw);
stat.time = stat.time-stat.time(1);

figure
set(gcf,'units','points','position',[10,10,300,400])

% figure
for ti = 1 : Nt
    t_range = (ti-1)*tw+1 : (ti-1)*tw+tw ;

    cfg.xlim = [stat.time(t_range(1)) stat.time(t_range(end))];   % time interval of the subplot

    %find channels
    ch_cluster = [];
    for ch = 1 : 32
        M = mask(ch,t_range);
        mx = find(M ~= 0);
        if ~isempty(mx)
            ch_cluster = [ch_cluster ch];
        end
    end

    ch_cluster2 = [];
    for ch = 1 : 32
        M = mask2(ch,t_range);
        mx = find(M ~= 0);
        if ~isempty(mx)
            ch_cluster2 = [ch_cluster2 ch];
        end
    end    
    % plot

    cfg.highlightchannel = unique([ch_cluster ch_cluster2]);
    cfg.highlight = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 14;
    cfg.comment = sprintf('[%d %d ms]', floor(1000*cfg.xlim(1)), ceil(1000*cfg.xlim(2)));
    cfg.commentpos = 'title';
    cfg.layout = 'gGamma32.lay';
    cfg.interactive = 'no';
    cfg.interplimits = 'head';
    cfg.style = 'straight';
    cfg.parameter = 'avg';
    cfg.colorbar = 'yes';
    
    subplot(3,2,ti)  
    ft_topoplotTFR(cfg, raweffect);
    colormap([0 0 0.5625;0 0 0.635416686534882;0 0 0.708333313465118;0 0 0.78125;0 0 0.854166686534882;0 0 0.927083313465118;0 0 1;0.0218803435564041 0.0711111128330231 0.978119671344757;0.0437606871128082 0.142222225666046 0.956239342689514;0.0656410306692123 0.213333338499069 0.934358954429626;0.0875213742256165 0.284444451332092 0.912478625774384;0.109401717782021 0.355555564165115 0.890598297119141;0.131282061338425 0.426666676998138 0.868717968463898;0.150427356362343 0.488888889551163 0.849572658538818;0.169572666287422 0.551111102104187 0.830427348613739;0.18871796131134 0.613333344459534 0.81128203868866;0.207863256335258 0.67555558681488 0.792136788368225;0.227008566260338 0.737777769565582 0.772991478443146;0.246153861284256 0.800000011920929 0.753846168518066;0.321323692798615 0.819943010807037 0.778391420841217;0.396493554115295 0.839886069297791 0.802936673164368;0.471663385629654 0.859829068183899 0.827481925487518;0.546833217144012 0.879772126674652 0.852027177810669;0.622003078460693 0.899715125560761 0.87657243013382;0.697172880172729 0.919658184051514 0.90111768245697;0.77234274148941 0.939601182937622 0.925662934780121;0.804865181446075 0.948229610919952 0.936282515525818;0.837387681007385 0.956857979297638 0.946902096271515;0.86991012096405 0.965486407279968 0.957521677017212;0.90243262052536 0.974114775657654 0.968141257762909;0.934955060482025 0.982743203639984 0.978760838508606;0.967477560043335 0.99137157201767 0.989380419254303;1 1 1;1 1 0.95714282989502;1 1 0.914285719394684;1 1 0.871428549289703;1 1 0.828571438789368;1 1 0.785714268684387;1 1 0.742857158184052;1 1 0.699999988079071;1 1 0.583333313465118;1 1 0.466666668653488;1 1 0.349999994039536;1 1 0.233333334326744;1 1 0.116666667163372;1 1 0;1 0.944444417953491 0;1 0.888888895511627 0;1 0.833333373069763 0;1 0.777777791023254 0;1 0.722222208976746 0;1 0.666666686534882 0;1 0.555555582046509 0;1 0.444444447755814 0;1 0.333333343267441 0;1 0.222222223877907 0;1 0.111111111938953 0;1 0 0;0.916666686534882 0 0;0.833333313465118 0 0;0.75 0 0;0.666666686534882 0 0;0.583333313465118 0 0;0.5 0 0]);
    caxis([-6 6]);  
end
