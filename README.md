# mTRF_alice

## Goal 
<li> Apply the multivariate temporal response function (mTRF) to Alice datasets (33 subjects)

## Package
<li> <a href= "http://www.fieldtriptoolbox.org/ "> fieldtrip toolbox  </a>for eeg data anaylsis </li>
<li> <a href= "https://www.mathworks.com/matlabcentral/fileexchange/74260-mtrf-toolbox "> mTRF-Toolbox  </a>for building mTRF models </li>
   
## Scripts versions
### do_alice_mtrf.m (version: 2020/07/30)
<li> Release version </li>

### do_alice_mtrf.m (version: 2020/07/15)
<li> Major changes: </li>
<li> baseline correction </li>
<pre><code>baseline = M(:,[91:129],:); %[-0.3 0]
baseline_average = nanmean(baseline,2);
baseline_average1 = repmat(baseline_average, [1 385 1]);
M = M-baseline_average1;</code></pre> 

### do_alice_mtrf.m (version: 2020/06/15)
<li> Major changes: </li>
<li> interpolation </li>
<pre><code>new_values  = interp1(time, value, new_times, 'previous', 0);</code></pre>
<li> transform variables (LogFreq and Closebrackets) based on WordOnset100ms and LexFunc </li>
<pre><code>% test -- transform variables based on "just_content" parameter
for i_cnt = 1:length(just_content)
   I_lex = find(strcmp('LexFunc', stim_raw.label));
   I_onset = find(strcmp('WordOnset100ms', stim_raw.label));
   I     = find(strcmp(just_content(i_cnt), stim_raw.label));
   stim_raw.trial{1}(I,:) = stim_raw.trial{1}(I,:) .* stim_raw.trial{1}(I_lex,:).*stim_raw.trial{1}(I_onset,:);
end</code></pre>

### do_alice_mtrf.m (version: 2020/05/31)
<li> Major changes: </li>
<li> Low-pass filter: 12 Hz </li>
<pre><code>%% Make evoked data
cfg = [];
cfg.lpfilter    = 'yes';
cfg.lpfreq      = 12;
dat_evokd = ft_preprocessing(cfg, dat_cln); </code> </pre>
<li> Lambda range: 2^0, 2^2, ..., 2^20 </li>
<pre><code>k = linspace(0,20,11);
for i=1:length(k)
    try_lambda(i) = 2^k(i);
end</code></pre>
<li> Time lags: -1000~2000ms </li>
<pre><code>[model] = mTRFtrain(stim, resp, test_Fs, 1, -1000, 2000, lambda);</code></pre>

### mtrf_plot_jrb.m  (version: 2020/07/30)
<li> Release version </li>
<li> Plot results across subjects. Each line represents each channel. </li>

## Files 
<li> <a href="https://umich.box.com/s/tw206e6kid6pj6og5vgsrkhdroihvlmb" > Predictors</a> </li>

## Model pilot results
<a href="https://docs.google.com/document/d/19UscK-aBd9DBrC2d08MNdrNf46zX557kOHEIKHxp_uQ/edit?usp=sharing" > Plots (test log)</a>
<li> <a href= "https://docs.google.com/presentation/d/1ksen6Z7AjV4sGlXhzdczhrRMcbfHzUM_n95v5T1ZtuY/edit?usp=sharing"> Model results 2020/07/01 (log frequency vs. frequnecy bins)</a> </li>
<li> <a href= "https://docs.google.com/presentation/d/1AyeqNDTFWX9w-bKPvfFp32O2WdFZSkCLJI9JxGzhxi8/edit?usp=sharing"> Model results 2020/07/07 (baseline correction)</a> </li>

## Analysis and preliminary results 
<a href="https://docs.google.com/document/d/18uYyKKKWDdwy63MBitztqVP91SU9sQGMuj_efAxRaYE/edit?usp=sharing"> Analysis and results </a>
