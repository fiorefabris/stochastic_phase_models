# adler
Package for generating time series from the theta phase model and some stochasti variations, and data analysis and visualization, following the data analysis protocol detailed in [*Intermittent ERK oscillations downstream of FGF in mouse embryonic stem cells*](https://journals.biologists.com/dev/article/149/4/dev199710/274396/Intermittent-ERK-oscillations-downstream-of-FGF-in).


<ins>Functionalities</ins>

**Get time series module**
    - *compute time series*. 
    - *get simulated dt mp*.
    - *get_simulated_dt_fixed_mp*
    - *get_epsilon_plus_pop_mp*

**Plotting module**
    - *plot time series*
    - *plot fft*
    - *plot FPT*
    - *plot pulses* 
    - *plot quantifiers*
    - *plot activity*
    - *plot consecutiveness*
    - *plot dt* (imports get simulated dt)
    - *plot simulated dt*
    - *plot epsilon plus* (imports all_comninations from time_series)

**Pulse detection module**
    - *compute FFT* 
    - *compute FDT* 
    - *get consecutive* 
    - *compute pulses quantifiers* 
    - *compute pulse detection* 

