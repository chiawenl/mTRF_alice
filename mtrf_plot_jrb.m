%% GitHub version: updated 2020/05/31
% predictor: 1: intensity; 7: logFreq; 8: closebracket
% load data

files = dir('models_20200531/R0*.mat');


%% initialize
% time lag 117: -100~800
x_logfreq   = zeros(33,117,61);  %subject x time x channel
x_intensity = zeros(33,117,61);  %subject x time x channel
x_brackets = zeros(33,117,61);  %subject x time x channel

%time lag 385: -1000~2000
x_logfreq   = zeros(33,385,61);  %subject x time x channel
x_intensity = zeros(33,385,61);  %subject x time x channel
x_brackets = zeros(33,385,61);  %subject x time x channel


%% Load models
label = {};
for j = 1:length(files)
    load(['models_20200531/' files(j).name], 'M_ev', 'e');
    e_num         = str2double(e);
    [~, neworder] = sort(e_num); % order by channel name for grandaverage
    x_logfreq(j,:,:)   = M_ev(7,:,neworder);
    x_intensity(j,:,:) = M_ev(1,:,neworder);
    x_brackets(j,:,:)  = M_ev(8,:,neworder);
end

label  = e(neworder);
time   = load(['models_20200531/' files(1).name], 't');

%% Prep layout

cfg = [];
cfg.layout = 'easycapM10-acti61_elec.sfp';
cfg.center = 'yes';
lay = ft_prepare_layout(cfg);

%% DEBUG
% 
% test = struct();
% test.dimord = 'chan_time';
% test.time   = time.t;
% test.label  = label;
% 
% cfg         = [];
% cfg.layout  = lay;
% cfg.comment = 'no';
% cfg.xlim    = [100 200];
% 
% for i = 1:23
%    subplot(6,4,i)
%    test.avg    = squeeze(x_intensity(i,:,:));
%    ft_topoplotER(cfg, test);
%    title(num2str(i));
% end

%% Grand-average
avg_logfreq   = squeeze(mean(x_logfreq,1));
avg_intensity = squeeze(mean(x_intensity, 1));
avg_brackets  = squeeze(mean(x_brackets, 1));

% test baseline-correct

% Let M = (avg, e*time);
% t = time lag
% baseline = M(:, [tmin:tmax]);?
% baseline_average = mean(M, 2);
% baseline_average_to_subtract = repmat(baseline_average, [1, num_time_points]);?
% M_corrected = M - baseline_average_to_subtract

% avg_lf = mean(x_logfreq,1);
% avg_lf_b = repmat(avg_lf,[33,1]);
% M_c_lf = x_logfreq-avg_lf_b;
% avg_logfreq = squeeze(mean(M_c_lf,1));
% 
% avg_in = mean(x_intensity,1);
% avg_in_b = repmat(avg_in,[33,1]);
% M_c_in = x_intensity-avg_in_b;
% avg_intensity = squeeze(mean(M_c_in,1));

%% Prep grandaveraged data for plotting

load('with_added_trialinfo/R0150.mat', 'dat'); % load processed_data and change the structure to have
%topography
dat_mtrf = dat;
dat_mtrf = rmfield(dat_mtrf, 'trial');
dat_mtrf = rmfield(dat_mtrf, 'dof');
dat_mtrf = rmfield(dat_mtrf, 'var');
dat_mtrf = rmfield(dat_mtrf, 'trialinfo');
dat_mtrf.time = time.t;
dat_mtrf.label = label;
dat_mtrf.dimord = 'chan_time';


%% plot for log freq
dat_mtrf.avg = avg_logfreq';

cfg.zlim = [-.15 .15];

cfg = [];
cfg.layout = lay;
cfg.comment = 'no';
cfg.marker = 'no';
%cfg.zlim = [-.05 .05];
cfg.zlim = [-.15 .15];

subplot(2,6,1);
cfg.xlim = [200 300];
ft_topoplotER(cfg, dat_mtrf);
title ('200-300 ms')

subplot(2,6,2)
cfg.xlim = [300 400];
ft_topoplotER(cfg, dat_mtrf);
title ('300-400 ms')

subplot(2,6,3)
cfg.xlim = [400 500];
ft_topoplotER(cfg, dat_mtrf);
title ('400-500 ms')

subplot(2,6,4)
cfg.xlim = [500 600];
ft_topoplotER(cfg, dat_mtrf);
title ('500-600 ms')

subplot(2,6,5)
cfg.xlim = [600 700];
ft_topoplotER(cfg, dat_mtrf);
title ('600-700 ms')

subplot(2,6,6)
cfg.xlim = [700 800];
ft_topoplotER(cfg, dat_mtrf);
title ('700-800 ms')

subplot(2,6,[7, 12]);
plot(dat_mtrf.time, dat_mtrf.avg)
xlim([-1000 2000])
%xlim([-100 800])
title('33 subjects, Log Word Frequency');



%% plot for intensity 
dat_mtrf.avg = avg_intensity';

cfg.zlim = [-.15 .15];


subplot(2,4,1);
%cfg.xlim = [100 120];
cfg.xlim = [100 100];
ft_topoplotER(cfg, dat_mtrf);
title ('100 ms')

subplot(2,4,2);
%cfg.xlim = [100 120];
cfg.xlim = [140 140];
ft_topoplotER(cfg, dat_mtrf);
title ('140 ms')

subplot(2,4,3)
%cfg.xlim = [155 165];
cfg.xlim = [179 179];
ft_topoplotER(cfg, dat_mtrf);
title ('179 ms')

subplot(2,4,4)
cfg.xlim = [218 218];
ft_topoplotER(cfg, dat_mtrf);
title ('218 ms')

subplot(2,4,[5, 8]);
plot(dat_mtrf.time, dat_mtrf.avg);
%xlim([-100 800])
xlim([-1000 2000])
title('33 subjects, Intensity');

%% plot for brackets
dat_mtrf.avg = avg_brackets';

cfg = [];
cfg.layout = 'easycapM10-acti61_elec.sfp';
cfg.comment = 'no';
cfg.marker = 'no';
cfg.zlim = [-.1 .1];

subplot(2,6,1);
cfg.xlim = [200 300];
ft_topoplotER(cfg, dat_mtrf);
title ('200-300 ms')

subplot(2,6,2)
cfg.xlim = [300 400];
ft_topoplotER(cfg, dat_mtrf);
title ('300-400 ms')

subplot(2,6,3)
cfg.xlim = [400 500];
ft_topoplotER(cfg, dat_mtrf);
title ('400-500 ms')

subplot(2,6,4)
cfg.xlim = [500 600];
ft_topoplotER(cfg, dat_mtrf);
title ('500-600 ms')

subplot(2,6,5)
cfg.xlim = [600 700];
ft_topoplotER(cfg, dat_mtrf);
title ('600-700 ms')

subplot(2,6,6)
cfg.xlim = [700 800];
ft_topoplotER(cfg, dat_mtrf);
title ('700-800 ms')

subplot(2,6,[7, 12]);
plot(dat_mtrf.time, dat_mtrf.avg)
%xlim([-100 800])
xlim([-1000 2000])
title('33 subjects, Close Brackets');