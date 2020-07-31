%% GitHub version: updated 2020/07/30
% predictor: 1: intensity; 7: logFreq; 8: closebracket
% load data

files = dir('models/R0*.mat');


%% initialize
% time lag 117: -100~800
% x_logfreq   = zeros(33,117,61);  %subject x time x channel
% x_intensity = zeros(33,117,61);  %subject x time x channel
% x_brackets = zeros(33,117,61);  %subject x time x channel
% x_onset = zeros(33,117,61);

%time lag 385: -1000~2000
x_logfreq   = zeros(33,385,61);  %subject x time x channel
x_intensity = zeros(33,385,61);  %subject x time x channel
x_brackets = zeros(33,385,61);  %subject x time x channel
x_onset = zeros(33,385,61);

%% Load models
label = {};
for j = 1:length(files)
    load(['models/' files(j).name], 'M_ev', 'e');
    e_num         = str2double(e);
    [~, neworder] = sort(e_num); % order by channel name for grandaverage
    x_logfreq(j,:,:)   = M_ev(7,:,neworder);
    x_intensity(j,:,:) = M_ev(1,:,neworder);
    x_brackets(j,:,:)  = M_ev(8,:,neworder);
    x_onset(j,:,:)  = M_ev(9,:,neworder);
    x_f1(j,:,:) = M_ev(3,:,neworder);
    x_pitch(j,:,:) = M_ev(2,:,neworder);
end

label  = e(neworder);
time   = load(['models/' files(1).name], 't');

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
avg_onset = squeeze(mean(x_onset, 1));
avg_f1 = squeeze(mean(x_f1,1));
avg_pitch = squeeze(mean(x_pitch,1));


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

cfg = [];
cfg.layout = lay;
cfg.comment = 'no';
cfg.marker = 'no';
cfg.zlim = [-.02 .02];

subplot(2,8,1);
cfg.xlim = [200 300];
ft_topoplotER(cfg, dat_mtrf);
title ('200-300 ms')

subplot(2,8,2)
cfg.xlim = [300 400];
ft_topoplotER(cfg, dat_mtrf);
title ('300-400 ms')

subplot(2,8,3)
cfg.xlim = [400 500];
ft_topoplotER(cfg, dat_mtrf);
title ('400-500 ms')

subplot(2,8,4)
cfg.xlim = [500 600];
ft_topoplotER(cfg, dat_mtrf);
title ('500-600 ms')

subplot(2,8,5)
cfg.xlim = [600 700];
ft_topoplotER(cfg, dat_mtrf);
title ('600-700 ms')

subplot(2,8,6)
cfg.xlim = [700 800];
ft_topoplotER(cfg, dat_mtrf);
title ('700-800 ms')

subplot(2,8,7)
cfg.xlim = [800 900];
ft_topoplotER(cfg, dat_mtrf);
title ('800-900 ms')

subplot(2,8,8);
cfg.colorbar = 'EastOutside';
ft_topoplotER(cfg, dat_mtrf);

subplot(2,8,[9, 16]);
plot(dat_mtrf.time, dat_mtrf.avg)
xlim([-1000 2000]) % time lag -1000~2000ms
xlabel('Time lag (ms)')
ylabel('a.u.')
%xlim([-100 800]) % time lag -100~800ms
title('33 subjects, Word Frequency');


%% plot for intensity 
dat_mtrf.avg = avg_intensity';

cfg = [];
cfg.layout = lay;
cfg.comment = 'no';
cfg.marker = 'no';
cfg.zlim = [-.04 .04];

subplot(2,5,1);
cfg.xlim = [100 100];
ft_topoplotER(cfg, dat_mtrf);
title ('100 ms')

subplot(2,5,2);
cfg.xlim = [140 140];
ft_topoplotER(cfg, dat_mtrf);
title ('140 ms')

subplot(2,5,3)
cfg.xlim = [179 179];
ft_topoplotER(cfg, dat_mtrf);
title ('179 ms')

subplot(2,5,4)
cfg.xlim = [218 218];
ft_topoplotER(cfg, dat_mtrf);
title ('218 ms')

subplot(2,5,5);
cfg.colorbar = 'EastOutside';
ft_topoplotER(cfg, dat_mtrf);

subplot(2,5,[6, 10]);
plot(dat_mtrf.time, dat_mtrf.avg);
xlim([-1000 2000])
xlabel('Time lag (ms)');
ylabel('a.u.');
title('33 subjects, Intensity');

%% plot for brackets
dat_mtrf.avg = avg_brackets';

cfg = [];
cfg.layout = 'easycapM10-acti61_elec.sfp';
cfg.comment = 'no';
cfg.marker = 'no';
cfg.zlim = [-.03 .02];

subplot(2,8,1);
cfg.xlim = [200 300];
ft_topoplotER(cfg, dat_mtrf);
title ('200-300 ms')

subplot(2,8,2)
cfg.xlim = [300 400];
ft_topoplotER(cfg, dat_mtrf);
title ('300-400 ms')

subplot(2,8,3)
cfg.xlim = [400 500];
ft_topoplotER(cfg, dat_mtrf);
title ('400-500 ms')

subplot(2,8,4)
cfg.xlim = [500 600];
ft_topoplotER(cfg, dat_mtrf);
title ('500-600 ms')

subplot(2,8,5)
cfg.xlim = [600 700];
ft_topoplotER(cfg, dat_mtrf);
title ('600-700 ms')

subplot(2,8,6)
cfg.xlim = [700 800];
ft_topoplotER(cfg, dat_mtrf);
title ('700-800 ms')

subplot(2,8,7)
cfg.xlim = [800 900];
ft_topoplotER(cfg, dat_mtrf);
title ('800-900 ms')

subplot(2,8,8);
cfg.colorbar = 'EastOutside';
ft_topoplotER(cfg, dat_mtrf);

subplot(2,8,[9, 16]);
plot(dat_mtrf.time, dat_mtrf.avg)
%xlim([-100 800])
xlim([-1000 2000])
xlabel('Time lag (ms)');
ylabel('a.u.');
title('33 subjects, Close Brackets');

%% plot for word onset
dat_mtrf.avg = avg_onset';

cfg = [];
cfg.layout = 'easycapM10-acti61_elec.sfp';
cfg.comment = 'no';
cfg.marker = 'no';
cfg.zlim = [-.04 .04];

subplot(2,8,1);
cfg.xlim = [200 300];
ft_topoplotER(cfg, dat_mtrf);
title ('200-300 ms')

subplot(2,8,2)
cfg.xlim = [300 400];
ft_topoplotER(cfg, dat_mtrf);
title ('300-400 ms')

subplot(2,8,3)
cfg.xlim = [400 500];
ft_topoplotER(cfg, dat_mtrf);
title ('400-500 ms')

subplot(2,8,4)
cfg.xlim = [500 600];
ft_topoplotER(cfg, dat_mtrf);
title ('500-600 ms')

subplot(2,8,5)
cfg.xlim = [600 700];
ft_topoplotER(cfg, dat_mtrf);
title ('600-700 ms')

subplot(2,8,6)
cfg.xlim = [700 800];
ft_topoplotER(cfg, dat_mtrf);
title ('700-800 ms')

subplot(2,8,7)
cfg.xlim = [800 900];
ft_topoplotER(cfg, dat_mtrf);
title ('800-900 ms')

subplot(2,8,8);
cfg.colorbar = 'EastOutside';
ft_topoplotER(cfg, dat_mtrf);

subplot(2,8,[9, 16]);
plot(dat_mtrf.time, dat_mtrf.avg)
%xlim([-100 800])
xlim([-1000 2000])
xlabel('Time lag (ms)');
ylabel('a.u.');
title('33 subjects, Word Onset 200ms');
