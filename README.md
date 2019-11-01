# PeriodicMeanSubtractionForHumNoise
This subtraction filter analyses the shape of the power line interference ("Hum noise") and subtracts it from the time series. Method: moving estimate of Hum noise by periodic median of the high-pass filtered signal. Highly robust to any harmonic contributions and transients in the signal.

%% Syntax:
% Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod)
% Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod, Window)
% Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod, Window, Plot)
% Filtered,Filter] = periodicMedianFilter(RawData, HumPeriod, Window, Plot, PowerLineHum)
%
%% mandatory parameters:
% RawData   = series of data with power line hum
% Humperiod = the (whole) number of samples over which the hum repeats (e.g. 1000/50=20 if sample freq=1000)
%
%% optional parameters
% Window [optional] = the number of periods on the sliding window (default 50; 0=entire range)
% Plot   [optional] default 0=none; 1=make a comparison plot and FFT; 2=also plot the filter
%
