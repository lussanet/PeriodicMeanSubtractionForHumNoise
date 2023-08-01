function [Filtered,Filter,FgHandle] = periodicMedianFilter(RawData, HumPeriod, Window, DoPlot, PowerLineHum)
%% Remove power line interference ("mains hum") -with harmonics- from signal.
% This subtraction filter analyses the shape of the power line interference ("Hum noise") and 
% subtracts it from the time series RawData.
% Method: moving estimate of Hum noise by periodic median of the high-pass filtered signal.
% 
% This method is highly robust to any harmonic contributions, frequency and amplitude flutuations 
% and transients in the signal. It performs really good at signal onset and offset.
% 
% The optimal window is about 1 sec, but if the data have a strong frequency component at 1 Hz (e.g.
% tuck jumps), the filter can fail dramatically. A check and solution (halving the window) is
% implemented to deal with this case.
%
% SYNTAX
% Filtered = periodicMedianFilter(RawData, HumPeriod);
% [Filtered,Filter,FgHandle] = periodicMedianFilter(RawData, HumPeriod, Window, Plot, PowerLineHum);
%
% INPUT
%   mandatory parameters:
%     RawData      (n-D x NSamples double) series of data with power line hum
%                  (If  NSamples x n-D,  this is corrected automatically)
%     Humperiod    (double) the (whole) number of samples over which the hum repeats (e.g. 
%                  1000/50=20 if sample freq=1000)
%   optional parameters
%     Window       [double] the number of periods on the sliding window (default 50; 0=entire range)
%     Plot         [double] default 0. 1=A comparison plot and FFT; 2=also plot the filter
%     PowerLineHum [double] the frequency of the hum (default 50 Hz)
%
% OUTPUT
%     Filtered  (n-D x NSamples double) clean data
%     Filter    (n-D x NSamples double) shape of the subtraction filter
%     FgHandle  (handle) figure handle
% 
% EXAMPLE
%  Frequency = 1500; % Hz
%  Mains     = 60; % Hz
%  TwoSconds = 2*Frequency;
%  DoPlot    = true;
% Filtered = periodicMedianFilter(RawData, round(Frequency/Mains), TwoSeconds, DoPlot, Mains);
% 
% Local functions: tryPeriodicMedianFilter
% 
% (c) 2023 by Movement Science, WWU Muenster
% Author: Marc de Lussanet
% Version 1  (15. 5.2019)
% Version 15 (27.01.2023) comments & header; check dimensions, etc.

%% handle optional parameters
narginchk(2,5);
if nargin<3 || isempty(Window),       Window       = 50;    end
if nargin<4 || isempty(DoPlot),       DoPlot       = false; end
if nargin<5 || isempty(PowerLineHum), PowerLineHum = 50;    end

%% Init
Filtered = RawData;
Filter   = [];

%% check that hum period has no decimals
if abs(HumPeriod-round(HumPeriod)) > 0.00001
    error('HumPeriod (%f) must be integer',HumPeriod);
else
    HumPeriod = round(HumPeriod);
end

%% check dimensions
Dims = size(RawData);
if length(Dims)>2
    error('Too many dimensions of RawData (%d > 2)',length(Dims));
end
NSamples = max(Dims);
if NSamples < HumPeriod
    warning('Hum filter not applied since no of samples (%d) less than period of hum (%d)',NSamples,HumPeriod);
    return;
end
% if the data are of NSample x n-D, invert the data
if Dims(1)>Dims(2)
    InvertData = true;
    RawData = RawData';
else
    InvertData = false;
end

%% Test whether the current window works for the file
% The periodic median filter has almost always excellent performance, except when the signal has a
% main frequency with a period in the range of the window (i.e. typically 1 sec). This can be the
% case e.g., for tuck jumps.
% In that case (i.e. if a warning is produced), try reduce the window. If that does not help, no
% filtering is applied.
[Filtered,Filter,FgHandle,Warn] = tryPeriodicMedianFilter(RawData, HumPeriod, Window, DoPlot, PowerLineHum);
if any(Warn,'all')
    [Filtered,Filter,FgHandle,Warn] = tryPeriodicMedianFilter(RawData, HumPeriod, Window/2, DoPlot, PowerLineHum);
    if any(Warn,'all')
        warning('Reducing the window did not help: no hum filtering is applied')
        Filtered = RawData;
    end
end

% if the data are of NSample x n-D, invert the data to the original dimension
if InvertData == true
    Filtered = Filtered';
end
end


%% =================================================================================================

function [Filtered,Filter,FgHandle,Warn] = tryPeriodicMedianFilter(RawData, HumPeriod, Window, DoPlot, PowerLineHum)
%% Subfunction that does the actual computation

%% handle optional parameters
narginchk(2,5);
if nargin<3 || isempty(Window),       Window       = 50;    end
if nargin<4 || isempty(DoPlot),       DoPlot       = false; end
if nargin<5 || isempty(PowerLineHum), PowerLineHum = 50;    end

%% Constants
HPCutFreq  = 20; % high pass frequency for constructing the subtraction filter
HPOrder    = 4 ; % (effective order is doubled by filtfilt)
MessFreq   = PowerLineHum * HumPeriod;
Sz         = size(RawData);
PlCh       = 1;  % channel that is plotted

%% Init
Filtered= RawData;
Filter  = [];
FgHandle= [];
Warn    = '';

%% error handling
if length(Sz)>2 || Sz(2)<Sz(1)
    error('Data must be shaped as Data(channels : timeseries)');
end
if length(RawData) < 2^nextpow2(PowerLineHum)/2
    disp('data too short for filtering'); return;
end
if any(all(isnan(RawData),2))
    % The data contains only NaN values in at least one dimension: there is nothing to be filtered
    return;
end
Gaps = isnan(RawData);
if any(isnan(RawData),'all')
    % Data contains NANs; these are filled with mean and removed again after filtering
    for Ch=1:Sz(1)
        RawData(Ch,Gaps(Ch,:)) = mean(RawData(Ch,:),'omitnan');
    end
end

%% prepare the signal: remove offset and drift
% high-pass filter removes drifts and variations. 2nd Order (4th order effectively) and cutoff
% not too close to hum frequency, to prevent deformations.
Signal  = filth(HPCutFreq,MessFreq,HPOrder,RawData,'h');

%% if window longer than data, then simply take the entire signal (added 190902)
if HumPeriod*Window > length(Signal)
    Window=0;
end

%% create the filter. Either from entire signal or moving window
DoUseMovingWindow = Window ~= 0;
if DoUseMovingWindow
    for Ch=1:Sz(1)
        %% create a moving window of Window periods for which to compute a gradually changing filter
        % make a repeating array of the signal
        FiltMat  = reshape(repmat(Signal(Ch,:),[1,Window]),1,length(Signal)*Window);
        % make a matrix of this in which each column starts one period later in the signal
        RepeatLen= HumPeriod*Window;
        FiltMat  = reshape([FiltMat FiltMat(1:RepeatLen)], [length(Signal) + HumPeriod,Window]);
        % the filter is the median
        Median = median( FiltMat,2)';
        % the end of the Filter contains increasingly more of the start of the signal: omit this
        Len      = size(FiltMat,1) - RepeatLen;
        Median(Len+1 : end) = [];
        % Repeat the first and the last Periods and add those to the filter
        Vor      = repmat(Median(1:HumPeriod)        ,1,ceil((Window-1)/2));
        Nach     = repmat(Median(end-HumPeriod+1:end),1,2*Window);
        if Ch==1,      Filter=zeros(Sz(1),length(Vor)+length(Median)+length(Nach));      end
        Filter(Ch,:) = [Vor Median Nach]; %#ok<AGROW> 
        Filter(Ch,:) = Filter(Ch,:)-mean(Filter(Ch,:)); %#ok<AGROW> %zero-offset
    end
    % subtract
    Filtered = Filtered-Filter(:, 1:length(RawData));
    FilterRep= squeeze(Filter(PlCh, 1:length(RawData)));
    NoiseAmplitude = max(abs(Filter(:, 1:length(RawData))),[],2);
else
    for Ch=1:Sz(1)
        if Ch==1,      Filter=zeros(Sz(1),HumPeriod);      end
        %% chop the signal in pieces of one period and take the median as filter
        NPer           = floor(length(Signal)/HumPeriod);
        Repeat         = reshape(Signal(Ch,1:HumPeriod*NPer),HumPeriod,NPer);
        Filter(Ch,:)   = median(Repeat,2)'; %#ok<AGROW> 
        Filter(Ch,:)   = Filter(Ch,:)-mean(Filter(Ch,:)); %#ok<AGROW> %zero-offset
        % repeat the filter for all periods
        FilterRep      = reshape(squeeze(repmat(Filter(Ch,:),1,NPer+1)),1,HumPeriod*(NPer+1));
        % ... and subtract it from the signal
        Filtered(Ch,:) = RawData(Ch,:) - FilterRep(1:length(Signal));
    end
    NoiseAmplitude = max(abs(Filter),[],2);
end
% remove the gaps again
Filtered(Gaps) = nan;
% report warnings
Error = abs(RawData-Filtered);
Warn = Error>2*NoiseAmplitude;
if any(Warn,'all')
    warning('The filter is not well suited for quickly repeated jumps (mean error = %f). Try using a shorter window',mean(Error,'all'));
end

%% figures if desired
if DoPlot
    RawCh = RawData(PlCh,:); FiltCh = Filtered(PlCh,:); %(this will plot the first channel only!)
    %% Frequency spectrum
    NFFT     = 2^nextpow2(length(FiltCh)); % Next power of 2 from length of y
    Freqs    = MessFreq/2*linspace(0,1,NFFT/2+1);
    Time     = (0:length(FiltCh)-1) / MessFreq;
    PowerVor = fft(RawCh -mean(RawCh),NFFT)/length(RawCh);
    PowerNach= fft(FiltCh-mean(RawCh),NFFT)/length(FiltCh);
    PowerVor = 2*abs(PowerVor( 1:NFFT/2+1));
    PowerNach= 2*abs(PowerNach(1:NFFT/2+1));
    
    FgHandle=figure('name','periodicMedianFilter');
    subplot(2,1,1);hold on; plot(Time,RawCh); plot(Time,FiltCh);
    title('Dataseries before (blue) and after PMS filter');   xlabel('time (s)');
    subplot(2,1,2);hold on; title('FFT spectrum before (blue) and after PMS filter')
    xlabel('frequency (s^{-1})'); ylabel('power');
    plot(Freqs,PowerVor);
    plot(Freqs,PowerNach);
    axis([0 inf 0 1.05*max(PowerVor(Freqs>PowerLineHum*0.9 & Freqs<PowerLineHum+1))]);
    if DoPlot==2
        figure;hold on; plot(RawData); plot(FilterRep);   title('dataseries / filter')
    end
end
end