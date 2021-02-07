function qc = alice_qualitycheck(dataset)
%% data quality check for alice. 
% updated: 2021/02/03
% can run either raw data or matfile
% impedence plot added

%% load data
if nargin < 1
    [file path] = uigetfile('*.eeg'); 
    dataset = [path file];
end

if contains (dataset, '.mat')
    load(dataset,'proc')
end

if ~exist('refs')
    refs = {'25', '29'};
end

if exist('proc')
    if isempty(proc.impedence.imps)
        hasimps = 0;
    else
        hasimps = 1;
    end
end

if ~exist('proc')
    sidx                        = regexp(dataset, 'R[0-9]{4,4}'); % find a 4 digit channel number preceeded by R
    proc.subject                = dataset(sidx(1):sidx(1)+4);
    
%     % for anonymous data
%     sidx                        = regexp(dataset, 'S[0-9]'); 
%     proc.subject                = dataset(sidx(1):sidx(1)+2);
    
    proc.dataset                = dataset;

    % run these checks on the FIRST file if multiple are supplied
    if iscell(proc.dataset)
        checkdata = proc.dataset{1};
    else
        checkdata = proc.dataset;
    end
    try
        [bads imps labels] = get_high_impedence(dataset, 25);
        hasimps = 1;
    catch
        hasimps = 0;
    end
    
    % Track reference information
    proc.implicitref = '29';
    proc.refchannels = refs;
    
    if hasimps == 1;
    % Track impedence information
        proc.impedence.bads = bads;
        proc.impedence.imps = imps;
        proc.impedence.labels = labels;
    end

end


%% load raw & preprocess
channels                                = {'all', '-VEOG', '-AUD', '-Aux5', '-OPTO'};

cfg                                     = [];
cfg.dataset                             = proc.dataset;
cfg.channel                             = channels;
cfg.reref                               = 'yes';
cfg.refchannel                          = proc.refchannels; % linked mastoids
cfg.implicitref                         = proc.implicitref;
cfg.hpfreq                              = 0.1;
cfg.hpfiltord                           = 6;
cfg.hpinstabilityfix                    = 'split';
cfg.hpfilter                            = 'yes';
cfg.dftfilter                           = 'yes';
cfg.dftfreq                             = [60 120 180];

    dat_raw = ft_preprocessing(cfg);
    
%% Define 1 sec trials
cfg_trl                      = [];
cfg_trl.dataset              = proc.dataset;
cfg_trl.trialdef.triallength = 1;
cfg_trl                      = ft_definetrial(cfg_trl);

    dat_all = ft_redefinetrial(cfg_trl, dat_raw);

    %% Variance matrix
varmat = zeros(length(dat_all.label), length(dat_all.trial));
for e = 1:length(dat_all.label)
    for t = 1:length(dat_all.trial)
        varmat(e,t) = var(dat_all.trial{t}(e,:));
    end
end

% handle noisy data effectively
% 1. find variance per chan
% 2. find variance per trial, rmv'd noisy chans
% 3. if either steps 1 or 2 zero out the data, don't do any rejections

var_chans = mean(varmat, 2);
if sum(var_chans < 1000) 
    pick_chans_idx = var_chans < 1000;
    pick_chans = dat_all.label(pick_chans_idx);
    isexclC = 1;
else
    pick_chans_idx = 1:length(dat_all.label);
    pick_chans = dat_all.label;
    isexclC = 0;
end

var_trials = mean(varmat(pick_chans_idx,:), 1);
if sum(var_trials < 1000) ~= 0
    pick_trials = find(var_trials < 1000);
    isexclT = 1;
else
    pick_trials = 1:length(var_trials);
    isexclT = 0;
end

%% Freq Average

cfg = [];
cfg.preproc.reref = 'yes';
cfg.preproc.refchannel = pick_chans;
cfg.method = 'mtmfft';
cfg.tapsmofrq = 1;
cfg.foilim = [1 80];
cfg.channel = pick_chans;
cfg.trials = pick_trials;
frq = ft_freqanalysis(cfg, dat_all);

%% plots
h(1) = subplot(3,2,1); % impeds

if hasimps
    
    if max(proc.impedence.imps) < 35
        ylimits = [0 35];
    elseif max(proc.impedence.imps) > 100
        ylimits = [0 100];
    else
        ylimits = [0 max(proc.impedence.imps)];
    end

    gndidx = find(strcmp(proc.impedence.labels, 'GND'));
    ref1idx = find(strcmp(proc.impedence.labels, 'REF_29'));
    ref2idx = find(strcmp(proc.impedence.labels, '25'));

    bar(proc.impedence.imps, 'k', 'edgecolor', 'w')
    ylim(ylimits);
    xlim([0 length(proc.impedence.imps)+1]);
    hline(25,'r');
    hold on
    bar(gndidx, proc.impedence.imps(gndidx),  'g');
    bar(ref1idx, proc.impedence.imps(ref1idx),  'b');
    bar(ref2idx, proc.impedence.imps(ref2idx),  'b');
    title({'Impedences', '{\color{blue} REF } {\color{green} GND } {\color{red} high! }'});
    ylabel('kOhm'); xlabel('channel');
    for b = 1:length(proc.impedence.imps)
        if proc.impedence.imps(b) > 25
            %bar(b, imps(b),  'r');
            text(b, proc.impedence.imps(b)+1, proc.impedence.labels{b}, 'color', 'r');
        end
    end

    hold off
    
else  % no impedences!
    plot([0 1], [0 1], 'w');
    axis off
    text(0.2, 0.5, 'no Impedences', 'fontsize', 24);
end


% variance matrix
h(2) = subplot(3,2,2); % variance matrix

imagesc(log(varmat), 'Parent',h(2));
title('log(var) by chan, trial'); ylabel('channel'); xlabel('trial');
colorbar

% variance barplots
h(3) = subplot(3,2,3); % variance:chans
bar(h(3), var_chans);
title('var by chan')
hline(600, 'b'); text(length(var_chans)+1, 600, '600', 'color', 'b');
hline(1000, 'r'); text(length(var_chans)+1, 1000, '1000', 'color', 'r');
xlim([0, length(var_chans)+1]);
set(gca,'YScale','log');

for b = 1:length(var_chans)
    if var_chans(b) > 1000
        text(b, var_chans(b)+1, dat_all.label{b}, 'color', 'r');
    end
end


h(4) = subplot(3,2,4); % variance:trials
bar(h(4), var_trials);
if isexclC
    title('var by trial, excl el var>1k');
else
    title('var by trial');
end
hline(600, 'b'); text(length(var_trials)+1, 600, '600', 'color', 'b');
hline(1000, 'r'); text(length(var_trials)+1, 1000, '1000', 'color', 'r');
xlim([0, length(var_trials)+1]);
set(gca,'YScale','log');

for b = 1:length(var_trials)
    if var_trials(b) > 1000
        text(b, var_trials(b)+1, num2str(b), 'color', 'r');
    end
end

% Frequency butterfply
h(5) = subplot(3,2,5); % freq butterfly
plot(frq.freq, frq.powspctrm);
set(gca,'YScale','log');
if isexclC || isexclT
    title('Spectrum, excl el | tr var>1k');
else
    title('Spectrum, no excl');
end
ylim([0, 1.1 * max(frq.powspctrm(:))])

dim = [0.607142857142857,0.191666666666667,0.1875,0.113095238095238];
str=['subject ID = ',proc.subject,newline,'pick chans = ',num2str(length(pick_chans)), newline,'pick trials = ',num2str(length(pick_trials)) ]
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%% save figs
saveas(gcf, ['figs_qc_alice/' proc.subject '_03.png']);
