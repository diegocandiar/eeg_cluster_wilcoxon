clc
clearvars
close all

fieldtrip_path = '...';
addpath(fieldtrip_path);
 
load('data_eeg_example.mat')
load('biosemi32_neighb.mat')

%%
cfg = [];
cfg.mint = 5;
cfg.statistic = 'dep_non_param';
cfg.alpha = 0.01/1;
cfg.minnbchan = 3; 
cfg.tail = 0; 
cfg.numrandomization = 10000;  
cfg.neighbours    = neighbours;
cfg.label = {'Fp1';'AF3';'F7';'F3';'FC1';'FC5';'T7';'C3';'CP1';'CP5';'P7';'P3';'Pz';'PO3';'O1';'Oz';'O2';'PO4';'P4';'P8';'CP6';'CP2';'C4';'T8';'FC6';'FC2';'F4';'F8';'AF4';'Fp2';'Fz';'Cz'};

stat = ft_statfun_wilcoxon(cfg, EEG_neutral_avg, EEG_trial) ; 

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

figure 
cfg  = [];
cfg.operation = 'x1-x2';
cfg.parameter = 'trial';
BHIdiff = ft_math(cfg,EEG_trial,EEG_neutral_avg);
BHIdiff.label = EEG_trial.label;

cfg = [];
cfg.removemean = 'no';
raweffect = ft_timelockanalysis2(cfg, EEG_trial);
    
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
time_range1 = stat.time(t_range) -30;  

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
time_range2 = stat.time(t_range) -30;     

%% VISUALIZE v2
raweffect.med = stat.stat;
mask = stat.posclusterslabelmat;
tw = 12;
Nt = floor(length(stat.time)/tw);
stat.time = stat.time-stat.time(1);

figure
set(gcf,'units','points','position',[10,10,600,150])

% figure
k = 1;
for nclu = [3 1 2 4]
    t_cluster = [];
    mask1 = mask==nclu;
    for t = 1 : length(stat.time)
        M = mask1(:,t);
        if length(find(M==1)) >= 1
            t_cluster = [t_cluster t];
        end
    end
    t_range = [min(t_cluster) max(t_cluster)];

    cfg.xlim = [stat.time(t_range(1)) stat.time(t_range(end))];   % time interval of the subplot

    %find channels
    ch_cluster = [];
    for ch = 1 : 32
        M = mask1(ch,t_range(1):t_range(2));
        mx = find(M == 1);
        if ~isempty(mx)
            ch_cluster = [ch_cluster ch];
        end
    end

    % plot

    cfg.highlightchannel = ch_cluster;
    cfg.highlight = 'on';
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 14;
    cfg.comment = sprintf('[%d %d s]', floor(cfg.xlim(1)), ceil(cfg.xlim(2)));
    cfg.commentpos = 'title';
    cfg.layout = 'biosemi32.lay';
    cfg.interactive = 'no';
    cfg.interplimits = 'head';
    cfg.style = 'straight';
    cfg.parameter = 'med';
    cfg.colorbar = 'yes';
    
    subplot(1,4,k)  
    ft_topoplotTFR(cfg, raweffect);
    colormap([0 0 0.5625;0 0 0.635416686534882;0 0 0.708333313465118;0 0 0.78125;0 0 0.854166686534882;0 0 0.927083313465118;0 0 1;0.0218803435564041 0.0711111128330231 0.978119671344757;0.0437606871128082 0.142222225666046 0.956239342689514;0.0656410306692123 0.213333338499069 0.934358954429626;0.0875213742256165 0.284444451332092 0.912478625774384;0.109401717782021 0.355555564165115 0.890598297119141;0.131282061338425 0.426666676998138 0.868717968463898;0.150427356362343 0.488888889551163 0.849572658538818;0.169572666287422 0.551111102104187 0.830427348613739;0.18871796131134 0.613333344459534 0.81128203868866;0.207863256335258 0.67555558681488 0.792136788368225;0.227008566260338 0.737777769565582 0.772991478443146;0.246153861284256 0.800000011920929 0.753846168518066;0.321323692798615 0.819943010807037 0.778391420841217;0.396493554115295 0.839886069297791 0.802936673164368;0.471663385629654 0.859829068183899 0.827481925487518;0.546833217144012 0.879772126674652 0.852027177810669;0.622003078460693 0.899715125560761 0.87657243013382;0.697172880172729 0.919658184051514 0.90111768245697;0.77234274148941 0.939601182937622 0.925662934780121;0.804865181446075 0.948229610919952 0.936282515525818;0.837387681007385 0.956857979297638 0.946902096271515;0.86991012096405 0.965486407279968 0.957521677017212;0.90243262052536 0.974114775657654 0.968141257762909;0.934955060482025 0.982743203639984 0.978760838508606;0.967477560043335 0.99137157201767 0.989380419254303;1 1 1;1 1 0.95714282989502;1 1 0.914285719394684;1 1 0.871428549289703;1 1 0.828571438789368;1 1 0.785714268684387;1 1 0.742857158184052;1 1 0.699999988079071;1 1 0.583333313465118;1 1 0.466666668653488;1 1 0.349999994039536;1 1 0.233333334326744;1 1 0.116666667163372;1 1 0;1 0.944444417953491 0;1 0.888888895511627 0;1 0.833333373069763 0;1 0.777777791023254 0;1 0.722222208976746 0;1 0.666666686534882 0;1 0.555555582046509 0;1 0.444444447755814 0;1 0.333333343267441 0;1 0.222222223877907 0;1 0.111111111938953 0;1 0 0;0.916666686534882 0 0;0.833333313465118 0 0;0.75 0 0;0.666666686534882 0 0;0.583333313465118 0 0;0.5 0 0]);
    caxis([-4 4]);  
    k=k+1;
end


%% 
figure
imagesc(stat.posclusterslabelmat)
colormap([1 1 1;0.968623220920563 0.991324782371521 0.981530964374542;0.937246441841125 0.982649624347687 0.963061928749084;0.905869662761688 0.973974406719208 0.944592833518982;0.874492883682251 0.965299189090729 0.926123797893524;0.843116044998169 0.956624031066895 0.907654762268066;0.811739265918732 0.947948813438416 0.889185726642609;0.780362486839294 0.939273595809937 0.870716691017151;0.748985707759857 0.930598437786102 0.852247595787048;0.71760892868042 0.921923220157623 0.833778560161591;0.686232149600983 0.913248002529144 0.815309524536133;0.654855370521545 0.90457284450531 0.796840488910675;0.623478591442108 0.895897626876831 0.778371453285217;0.592101812362671 0.887222409248352 0.759902358055115;0.560724973678589 0.878547251224518 0.741433322429657;0.529348194599152 0.869872033596039 0.722964286804199;0.497971415519714 0.86119681596756 0.704495251178741;0.466594636440277 0.852521657943726 0.686026155948639;0.43521785736084 0.843846440315247 0.667557120323181;0.403841078281403 0.835171222686768 0.649088084697723;0.372464269399643 0.826496064662933 0.630619049072266;0.341087490320206 0.817820847034454 0.612150013446808;0.309710711240768 0.809145629405975 0.593680918216705;0.278333932161331 0.800470471382141 0.575211882591248;0.246957138180733 0.791795253753662 0.55674284696579;0.238334774971008 0.783381938934326 0.578905701637268;0.229712411761284 0.774968683719635 0.601068556308746;0.221090033650398 0.766555368900299 0.623231410980225;0.212467670440674 0.758142113685608 0.645394265651703;0.203845307230949 0.749728798866272 0.667557120323181;0.195222944021225 0.741315484046936 0.689719974994659;0.186600565910339 0.732902228832245 0.711882829666138;0.177978202700615 0.724488914012909 0.734045684337616;0.169355839490891 0.716075658798218 0.756208539009094;0.160733476281166 0.707662343978882 0.778371453285217;0.152111113071442 0.699249029159546 0.800534307956696;0.143488734960556 0.690835773944855 0.822697162628174;0.134866371750832 0.682422459125519 0.844860017299652;0.126244008541107 0.674009203910828 0.86702287197113;0.117621637880802 0.665595889091492 0.889185726642609;0.108999274671078 0.657182574272156 0.911348581314087;0.100376904010773 0.648769319057465 0.933511435985565;0.0917545408010483 0.640356004238129 0.955674290657043;0.0831321701407433 0.631942749023438 0.977837145328522;0.0745098069310188 0.623529434204102 1;0.0705882385373116 0.614447891712189 0.986377716064453;0.0666666701436043 0.605366349220276 0.972755432128906;0.062745101749897 0.596284866333008 0.959133148193359;0.0588235333561897 0.587203323841095 0.945510864257812;0.0549019612371922 0.578121781349182 0.931888520717621;0.0509803928434849 0.569040238857269 0.918266236782074;0.0470588244497776 0.559958755970001 0.904643952846527;0.0431372560560703 0.550877213478088 0.89102166891098;0.0392156876623631 0.541795670986176 0.877399384975433;0.0352941192686558 0.532714128494263 0.863777101039886;0.0313725508749485 0.523632645606995 0.85015481710434;0.0274509806185961 0.514551103115082 0.836532533168793;0.0235294122248888 0.505469560623169 0.822910249233246;0.0196078438311815 0.496388047933578 0.809287965297699;0.0156862754374743 0.487306505441666 0.795665621757507;0.0117647061124444 0.478224992752075 0.78204333782196;0.00784313771873713 0.469143450260162 0.768421053886414;0.00392156885936856 0.460061937570572 0.754798769950867;0 0.450980395078659 0.74117648601532]);
box off
set(gcf,'units','points','position',[10,10,400,150])
set(gca, 'FontSize', 14)
yticks(1:32)
yticklabels({})
xticks(0:10:80)
xlabel('seconds'); ylabel('channel')

