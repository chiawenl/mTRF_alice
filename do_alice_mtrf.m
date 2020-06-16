function [M,t] = do_alice_mtrf(dataset)

%% function [M t] = do_alice_mtrf(dataset)
%
% See https://doi.org/10.3389/fnhum.2016.00604
% and https://doi.org/10.1080/23273798.2018.1502458
%
% 

%cd /Volumes/ling-cnl-stor/ana/AliceStory2/ana14-mtrf
%addpath('/Users/jobrenn/Documents/matlab/toolbox/mTRF-Toolbox-2020/mtrf')
%load(dataset, 'dat', 'proc');


%% TODO
% - Finish converting predictors to .txt format (in progress with f1.txt)
% - Decide what to do with stim/predictor GAPS: 0? NaN? other interpolation?

%% GitHub versions
% Updated: 2020/06/15

%% Set Parameters

% predictors are in tab-delimited .txt with three columns
% passage: integer for passage (1-12) 
% time:    times (sec) relative to wav onset
% value:   predictor value
%
% just_content: predictors to zero out for non-lexical words


use_predictors = {'Intensity', 'Pitch', 'F1', 'sentence', 'position', 'LexFunc', 'LogFrqHAL', 'CloseBrackets3', 'WordOnset100ms'};

just_content   = {'LogFrqHAL', 'CloseBrackets3'};

%stim_loc     = '../ana14-mtrf/stim/seg';
%predict_loc  = '../ana11-andrea-composition/alice_word_length_brackets_Nikki.xlsx';
%use_lingpred = {'sentence', 'position', 'LexFunc', 'LogFrqHAL', 'mod_count', 'CloseBrackets3'};
%use_audprped = {''};
%new_fs       = 200;  % not used

load(dataset, 'dat', 'proc');

%load('proc/R0182.mat', 'dat','proc'); % good fit
%load('with_added_trialinfo/R0225.mat', 'dat','proc'); % mTRF doesn't fit -> complex numbers

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

%% Define 10 sec trials
cfg_trl                      = [];
cfg_trl.dataset              = proc.dataset;
cfg_trl.trialdef.triallength = 10;
cfg_trl                      = ft_definetrial(cfg_trl);

    dat_all = ft_redefinetrial(cfg_trl, dat_raw);
    
    
%% Apply ICA to raw data
cfg            = [];
cfg.channel    = proc.ica.topolabel;
dat_all_forica = ft_selectdata(cfg, dat_all);


cfg           = [];
cfg.unmixing  = proc.ica.unmixing;
cfg.topolabel = proc.ica.topolabel;
    comp = ft_componentanalysis(cfg, dat_all_forica);

cfg           = [];
cfg.component = proc.ica.rejcomp; % Excluded Components
    dat_all_ica = ft_rejectcomponent(cfg, comp);
    
%% Remove artifactual epochs

cfg = [];
cfg.artfctdef.first = proc.rejections.first.artfctdef.summary;
cfg.artfctdef.final = proc.rejections.final.artfctdef.summary;
dat_cln      = ft_rejectartifact(cfg, dat_all_ica);


%% Replace bad channels in raw data
% Note: zero channels should also be reconstructed! (can happen to 1/2 of
% ref channels if other ref was noisy)
missing       = union(proc.impedence.bads, proc.rejections.badchans);
cln_data      = horzcat(dat_cln.trial{:});
zero_chan_idx = find(all(cln_data == 0, 2));

missing       = [missing, dat_cln.label{zero_chan_idx}];

if ~isempty(missing)
    cfg                                  = [];
    cfg.method                           = 'template';
    cfg.channel                          = {'all'};
    cfg.elecfile                         = 'easycapM10-acti61_elec.sfp';
    cfg.template                         = 'easycapM10-acti61_neighb.mat';
    neighbours                           = ft_prepare_neighbours(cfg);

    cfg = [];
    cfg.method                           = 'spline';
    if size(missing,1) > size(missing,2) % column vector
        cfg.badchannel                   = missing;
    else % row vector
        cfg.badchannel                   = missing';
    end
    cfg.neighbours                       = neighbours;
    cfg.elecfile                         = 'easycapM10-acti61_elec.sfp';

    %dat_raw_ica                           = ft_channelrepair(cfg, dat_all_ica);
    dat_cln                           = ft_channelrepair(cfg, dat_cln);
end

%% Create predictor time-series 
%
% goal: create matrix M (predictors x eeg_times)
%
% 1. per passage psg...
% 2. ...per predictor prd...
% 3. ...p_prd     = interp1(predictor_onsets, predictor_value)...
% 4. ...P_psg     = [p_prd_1... p_prd_n]...
%       where P_psg has columns times, predictor1, predictor2, etc...
% 5. ...P_psg.Samples = floor(P_psg.Times * sampleRate) + eeg_sample_of_passage_onset
% 6. ...P_all = [P_psg_1; P_psg_n]
% 7. convert P_all to fieldtrip "raw" data structure with 1 "continuous"
% trial
% 8-onward: repeat fieldtrip epoching and rejections, then fit mTRF in
% 10sec intervals
%
% matlab function: use interp1 with method = 'previous' 

% prep stimulus predictors
stim_raw          = dat_raw; % 
stim_raw.label    = use_predictors'; % as column vector
stim_raw.trial{1} = zeros(length(stim_raw.label), length(stim_raw.time{1}));
Fs                = stim_raw.fsample;
time_raw          = stim_raw.time{1};
% get passage onsets in samples
evt          = ft_read_event(proc.dataset, 'type', 'Stimulus');
psg_triggers = {evt(:).value};

if length(psg_triggers{1})>2 % handle trigger as "S  1" or "1" format
    psg_triggers = cellfun(@(x) str2num(x(2:4)), psg_triggers);
else 
    psg_triggers = cellfun(@(x) str2num(x), psg_triggers);
end

psg_onsets     = cell2mat({evt(:).sample});
psg_onsets_sec = psg_onsets / Fs;


for i_prd = 1:length(stim_raw.label) % for each predictor
    predictor = readtable(['predictors/' stim_raw.label{i_prd} '.txt'], 'TreatAsEmpty', '--undefined--');
    for i_psg = psg_triggers % for each passage
        passage = predictor(predictor.passage == i_psg,:);
        time    = passage.time;
        if any(isnan(time))                % handle missing data in time axis
            interval   = mode(diff(time)); % assume sampling is regular for missing time data
            if isnan(time(1))              % we need to reconstruct the start time
                first_good_time = find(~isnan(time), 1);
                start_time      = time(first_good_time) - (interval * (first_good_time-1));
            else
                start_time = time(1);
            end
            % create corrected time axis
            time_cor = start_time + (interval * (0:length(time)));
            time_cor = time_cor(1:(end-1))'; %... and make row vector
            
            % check our work
            check    = [time, time_cor, time - time_cor];
            if max(check(:,3)) >= 0.001 
                error('Reconstructed time axis is not accurate!');
            else % else all good!
                time = time_cor;
            end
        end
        value = fillmissing(passage.value, 'constant', 0);
        tmin = time(1);
        tmax = time(end);
%         tmin = time(1);
%         tmax = time(end);

        % interpolate to sampling rate of the data
        new_times   = 0:(1/Fs):tmax; % MAKE SURE NO NANs!
        %new_values  = interp1(time, value, new_times, 'nearest', 0);
        % interp1 for word onset/logfreq should be previous
        new_values  = interp1(time, value, new_times, 'previous', 0);
        
        
        % check!!
        %plot(time, value);
        %hold on
        %plot(new_times, new_values, 'r');
        %hold off
        
        % align with time-axis of the raw data
        first_sample = psg_onsets(psg_triggers==i_psg);
        last_sample  = first_sample + length(new_times)-1;
        new_samples  = first_sample : last_sample;
        
        stim_raw.trial{1}(i_prd, new_samples) = new_values;
    end   
end

% % transform variables based on "just_content" parameter
% for i_cnt = 1:length(just_content)
%    I_lex = find(strcmp('LexFunc', stim_raw.label));
%    I     = find(strcmp(just_content(i_cnt), stim_raw.label));
%    %stim_raw.trial{1}(I,:) = stim_raw.trial{1}(I,:) * stim_raw.trial{1}(I_lex,:);
%    stim_raw.trial{1}(I,:) = stim_raw.trial{1}(I,:) .* stim_raw.trial{1}(I_lex,:);
% end
% 
% % then transform variables based on 'just_onset' parameter
% for i_cnt = 1:length(just_content)
%     I_onset = find(strcmp('WordOnset100ms', stim_raw.label));
%     I     = find(strcmp(just_content(i_cnt), stim_raw.label));
%     stim_raw.trial{1}(I,:) = stim_raw.trial{1}(I,:) .* stim_raw.trial{1}(I_onset,:);
% end 

% test -- transform variables based on "just_content" parameter
for i_cnt = 1:length(just_content)
   I_lex = find(strcmp('LexFunc', stim_raw.label));
   I_onset = find(strcmp('WordOnset100ms', stim_raw.label));
   I     = find(strcmp(just_content(i_cnt), stim_raw.label));
   stim_raw.trial{1}(I,:) = stim_raw.trial{1}(I,:) .* stim_raw.trial{1}(I_lex,:).*stim_raw.trial{1}(I_onset,:);
end

% scale all predictors by SD
for i = 1:length(stim_raw.label)
    my_sd                  = std(stim_raw.trial{1}(i,:));
    stim_raw.trial{1}(i,:) = stim_raw.trial{1}(i,:) / my_sd;
end

% Epoch stim time-series
stim_all = ft_redefinetrial(cfg_trl, stim_raw);

% Apply rejections to stim time-series
cfg = [];
cfg.artfctdef.first = proc.rejections.first.artfctdef.summary;
cfg.artfctdef.final = proc.rejections.final.artfctdef.summary;
stim_cln            = ft_rejectartifact(cfg, stim_all);

% Remove "trials" where predictors are > 50% zeros
num_trials  = length(stim_cln.trial);
keep_trials = zeros(1,num_trials);
for t = 1:num_trials
    prp_zero = sum(stim_cln.trial{t}(:) == 0) / numel(stim_cln.trial{t});
    if prp_zero < 0.5
        keep_trials(t) = 1;
    end
end

stim_cln.trial      = stim_cln.trial(find(keep_trials));
stim_cln.time       = stim_cln.time(find(keep_trials));
stim_cln.sampleinfo = stim_cln.sampleinfo(find(keep_trials),:);

dat_cln.trial       = dat_cln.trial(find(keep_trials));
dat_cln.time        = dat_cln.time(find(keep_trials));
dat_cln.sampleinfo  = dat_cln.sampleinfo(find(keep_trials),:);


% DEBUG: plots

% for i = 1:length(stim_raw.label)
%     subplot(4,2,i)
%     plot(stim_raw.time{1}, stim_raw.trial{1}(i,:));
%     title(stim_raw.label{i});
% end

%% Make evoked data

cfg = [];
cfg.lpfilter    = 'yes';
%cfg.lpfreq      = 40;
cfg.lpfreq      = 12;
dat_evokd = ft_preprocessing(cfg, dat_cln);
%dat_evokd = ft_preprocessing(cfg, dat_raw_ica);

%% Delta 1-4 Hz 
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [1 4];
cfg.hilbert   = 'abs';
dat_delta = ft_preprocessing(cfg,dat_cln);

%% Theta 5-8
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [4 8];
cfg.hilbert   = 'abs';
dat_theta = ft_preprocessing(cfg,dat_cln);

%% Gamma 30-50
cfg = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [30 50];
cfg.hilbert   = 'abs';
dat_gamma = ft_preprocessing(cfg, dat_cln);


%% Crossvalidate to estimate lambda against the EVOKED data
% [R,P,RMSE] = MTRFCROSSVAL(STIM,RESP,FS,MAP,TMIN,TMAX,LAMBDA)

% downsample  for tractability
train_Fs       = 64;

cfg            = [];
cfg.resamplefs = train_Fs;
dat_evokd_rs   = ft_resampledata(cfg, dat_evokd);
stim_cln_rs    = ft_resampledata(cfg, stim_cln);

% rotate so features are columns
% ...and NaN -> zeros
stim = stim_cln_rs.trial;
resp = dat_evokd_rs.trial;
for i = 1:length(stim)
    stim{i} = stim{i}';
    resp{i} = resp{i}';
    stim{i}(isnan(stim{i})) = 0; 
    resp{i}(isnan(resp{i})) = 0; 
end

%try_lambda = logspace(1, 5, 30); %base
%try_lambda = logspace(-1,5,50);
% try_lambda = logspace(-1,6,100);
% try_lambda = logspace(0, 6, 20);
% lambda range
%k = linspace(0,20,11);
k = linspace(-15,15,31);
for i=1:length(k)
    try_lambda(i) = 2^k(i);
end

% old version mtrf toolbox script (2019)
% [R,P,RMSE]  = mTRFcrossval(stim, resp, 64, 1, 0, 100, try_lambda); 
% R           = squeeze(nanmean(nanmean(R, 3)));
% RMSE        = squeeze(nanmean(nanmean(RMSE, 3)));
% [~,I_R]     = max(R);
% [~,I_RMSE]  = min(RMSE);

% for new version mtrf toolbox script (2020)

%%% CONTINUE FROM HERE: correlation values are complex (!!) maybe because
%%% trying to take the SQRT of a negative when evaluating the fits?? see
%%% mTRFevaluate ln. 102 etc...
[stat, T]   = mTRFcrossval(stim, resp, train_Fs, 1, 0, 500, try_lambda);
R           = stat.r;
R           = squeeze(nanmean(nanmean(R, 3)));
RMSE        = stat.err;
RMSE        = squeeze(nanmean(nanmean(RMSE, 3)));
[~,I_R]     = max(R);
[~,I_RMSE]  = min(RMSE);

subplot(1,2,1)
plot(try_lambda, R, '-o', 'linewidth', 2)
hold on; plot(try_lambda(I_R), R(I_R), 'k.', 'markersize', 30); hold off;
set(gca, 'XScale', 'log')
subplot(1,2,2)
plot(try_lambda, RMSE, '-o', 'linewidth', 2, 'color', 'red')
hold on; plot(try_lambda(I_RMSE), RMSE(I_RMSE), 'k.', 'markersize', 30); hold off;
set(gca, 'XScale', 'log')

% Smallest lambda that either minimizes RMSE or maximizes R:
lambda = try_lambda(min(I_R, I_RMSE));

%% Train Models
% [model,t,c] = mTRFtrain(stim,resp,fs,map,tmin,tmax,lambda)
% z-score predictors before fitting
% train per trial, then average
% lambda set via cross-validation above...

% loop over four datasets!
dsets   = {'dat_evokd', 'dat_delta', 'dat_theta', 'dat_gamma'}; 
test_Fs = 128;

for d = 1:length(dsets)
    use_dat = eval(dsets{d});
    
    % downsample  for tractability
    cfg            = [];
    cfg.resamplefs = test_Fs;
    dat_cln_rs     = ft_resampledata(cfg, use_dat);
    if d == 1 % resample the stim only the first time around
        stim_cln_rs = ft_resampledata(cfg, stim_cln);
    end
    
    % M is epochs x predictors x time-lags x channels
    % time dimension = min lag to max lag in resampled intervals...
    % set this based on the time interval for the lags, set in the model
    % call in the loop below
    M = zeros(  length(dat_cln_rs.trial), ...  % epochs
                length(stim_cln_rs.label), ... % pedictors
                385, ...                       % time lags 117 for  time lag-100~800; 385 for -1000~2000
                length(dat_cln_rs.label));     % channels

    for i = 1:length(dat_cln_rs.trial) 
        stim              = stim_cln_rs.trial{i}(:,:)';
        stim              = (stim - nanmean(stim)) ./ nanstd(stim); % scale
        stim(isnan(stim)) = 0; %  NaN -> zero
        
        %resp              = dat_rs.trial{i}(:,:)';
        resp              = dat_cln_rs.trial{i}(:,:)';
        resp              = (resp - nanmean(resp)) ./ nanstd(resp); % scale
        
        % time interval for the model lags set here!
        % old version script
%         [model,t,c] = mTRFtrain(stim, resp, 128, 1, -100, 800, lambda);
%         M(i,:,:,:)  = model;

        % new 2020 version script
        %[model] = mTRFtrain(stim, resp, test_Fs, 1, -100, 800, lambda);
        [model] = mTRFtrain(stim, resp, test_Fs, 1, -1000, 2000, lambda);
        M(i,:,:,:)  = model.w;
    end
    
    M       = squeeze(nanmean(M, 1));
    varname = ['M' dsets{d}(4:6)];
    eval([varname ' = M;']); % creates M_ev, M_de, M_th, M_ga...
end
%chan = find(strcmp('1', dat.label));
%plot(t, squeeze(M(:,:,chan)), 'linewidth', 2); legend(pred.label, 'fontsize', 18, 'boxoff')

%% Save result
v = stim_cln_rs.label;
e = dat_cln_rs.label;
t = model.t; % 2020 version script
save(['models_20200615/'  proc.subject '_20200615.mat'], 'M_ev', 'M_de', 'M_th', 'M_ga', 't', 'v', 'e');
%save(['models/test.mat'], 'M_ev', 'M_de', 'M_th', 'M_ga', 't', 'v', 'c');
%% Figure
figure;
chan = find(strcmp('3', e));
%plot(t, squeeze(M_ev(:,:,chan)), 'linewidth', 2); legend(v, 'fontsize', 18, 'boxoff')
%plot(t, squeeze(M_ev(:,:,chan)), 'linewidth', 2); legend(v, 'fontsize', 18)
y =squeeze(M_ev(:,:,chan));
plot(t,y(9,:), 'linewidth', 2);legend(v{9}, 'fontsize', 18);
xlim([-1000, 2000]);
saveas(gcf, ['figs/' proc.subject '.png']);
%saveas(gcf, ['figs/test.png']);
