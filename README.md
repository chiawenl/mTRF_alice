# mTRF_alice

## Goal 
<li> Apply mTRF to Alice data (33 subjects)

## Scripts versions
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
<a href="https://umich.box.com/s/tw206e6kid6pj6og5vgsrkhdroihvlmb" > Predictors</a>

## Model results
<a href="https://docs.google.com/document/d/19UscK-aBd9DBrC2d08MNdrNf46zX557kOHEIKHxp_uQ/edit?usp=sharing" > Plots</a>

