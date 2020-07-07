# mTRF_alice

## Goal 
<li> Apply mTRF to Alice data (33 subjects)

## Scripts versions
### do_alice_mtrf.m (version: 2020/07/06)
<li> Major changes: </li>
<li> baseline correction </li>
<pre><code>baseline = M(:,[1:129],:); %pre-stimulus [-1000 0]
    other = M(:,[130:end],:); %post-stimulus [0 2000]
    baseline_average = nanmean(M(:,1:129,:)); %1x129x61
    baseline_average1 = repmat(baseline_average, [9, 1, 1]); %9x129x61
    M_corrected = baseline-baseline_average1; 
    M = horzcat(M_corrected,other);%9x385x61</code></pre> 

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

### mtrf_plot_jrb.m (version: 2020/05/31)
Plot results for Log frequency, Intensity, and Close brackets

## Files 
<li>Models based on do_alice_mrtf.m(2020/05/31)</li>
<a href="https://umich.box.com/s/tbxxkr33rnx7hxifh4zrlagwuq2t3r0f" > Models in Mbox</a> 
<li> <a href="https://umich.box.com/s/tw206e6kid6pj6og5vgsrkhdroihvlmb" > Predictors</a> </li>

## Model results
<a href="https://docs.google.com/document/d/19UscK-aBd9DBrC2d08MNdrNf46zX557kOHEIKHxp_uQ/edit?usp=sharing" > Plots</a>
<li> <a href= "https://docs.google.com/presentation/d/1ksen6Z7AjV4sGlXhzdczhrRMcbfHzUM_n95v5T1ZtuY/edit?usp=sharing"> Model results 2020/07/01</a> </li>
