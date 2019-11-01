function [Filtered,Filter] = removeHumByPMS(RawData, HumPeriod, Window, Plot, PowerLineHum)

%% This subtraction filter analyses the shape of the power line interference ("Hum noise") and subtracts it from the time series RawData.
% Method: moving estimate of Hum noise by periodic median of the high-pass filtered signal.
% Highly robust to any harmonic contributions and transients in the signal.
%
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
% Marc de Lussanet, Movement Science, WWU Muenster
% Version 1 (15.5.2019)
% Version 2 (21.8.2019) default window width 50 periods (1 sec)
% Version 3 (25.8.2019) improved high-pass filter settings
% Version 4 (02.9.2019) various bugs (see below)

	%% handle optional parameters
	if nargin<2,        error('RawData and Humperiod are required parameters');
	elseif nargin == 2, PowerLineHum = 50; Plot = 0; Window = 50;
	elseif nargin == 3, PowerLineHum = 50; Plot = 0;
	elseif nargin == 4, PowerLineHum = 50; 
	end
	
	%% Constants (25.8.2019)
	CutFreq  = 20; % was PowerLineHum/2;
	LPOrder  = 4 ; % was 2
	MessFreq = PowerLineHum * HumPeriod; 
	
	%% error handling
	if length(RawData) < 2^nextpow2(PowerLineHum)/2
		Filtered = [];Filter = []; 
		disp('data too short for filtering'); return;
	end
	if sum(isnan(RawData))
		warning('data contains NANs; Fix: these will be filled with mean');
		RawData(isnan(RawData)) = nanmean(RawData);
	end
	
	%% prepare the signal: remove offset and drift
	Signal  = RawData;
	% high-pass filter removes drifts and variations. 2nd Order (4th order effectively) and cutoff 
	% not too close to hum frequency, to prevent deformations.
	Signal  = filth(CutFreq,MessFreq,LPOrder,Signal,'h');
	
	%% if window longer than data, then simply take the entire signal (added 190902)
	if HumPeriod*Window > length(Signal)
		Window=0;
	end
	
	%% create the filter. Either from entire signal or moving window
	if ~Window
		%% chop the signal in pieces of one period and take the median as filter
		NPer     = floor(length(Signal)/HumPeriod);
		Repeat   = reshape(Signal(1:HumPeriod*NPer),HumPeriod,NPer);
		Filter   = median(Repeat,2)';
		Filter   = Filter-mean(Filter); %zero-offset
		% repeat the filter for all periods
		FilterRep= reshape(repmat(Filter,1,NPer+1),1,HumPeriod*(NPer+1));
		% ... and subtract it from the signal
		Filtered = RawData - FilterRep(1:length(Signal));
	else
		%% create a moving window of Window periods for which to compute a gradually changing filter
		% make a repeating array of the signal
		FiltMat  = reshape(repmat(Signal',[1,Window]),1,length(Signal)*Window);
		% make a matrix of this in which each column starts one period later in the signal
		RepeatLen= HumPeriod*Window;
		FiltMat  = reshape([FiltMat FiltMat(1:RepeatLen)], [length(Signal) + HumPeriod,Window]); % BUG 190902
		% the filter is the median 
		Filter   = median( FiltMat,2)';
		% the end of the Filter contains increasingly more of the start of the signal: omit this
		Len      = size(FiltMat,1) - RepeatLen;                                                  % BUG 190902
		Filter(Len+1 : end) = [];
		% Repeat the first and the last Periods and add those to the filter
		Vor      = repmat(Filter(1:HumPeriod)        ,1,ceil((Window-1)/2));
		Nach     = repmat(Filter(end-HumPeriod+1:end),1,2*Window);
		Filter   = [Vor Filter Nach];
		Filter   = Filter-mean(Filter); %zero-offset
		% subtract
		Filtered = RawData;
		Filtered = Filtered-Filter(1:length(RawData));
		FilterRep= Filter(1:length(RawData));
	end
	
	%% figures if desired
	if Plot
		%% Frequency spectrum
		NFFT     = 2^nextpow2(length(Filtered)); % Next power of 2 from length of y
		Freqs    = MessFreq/2*linspace(0,1,NFFT/2+1); 
		Time     = (0:length(Filtered)-1) / MessFreq;
		PowerVor = fft(RawData -mean(RawData),NFFT)/length(RawData);
		PowerNach= fft(Filtered-mean(RawData),NFFT)/length(Filtered);
		PowerVor = 2*abs(PowerVor( 1:NFFT/2+1));
		PowerNach= 2*abs(PowerNach(1:NFFT/2+1));
		
		figure;
		subplot(2,1,1);hold on; plot(Time,RawData); plot(Time,Filtered); 
		title('Dataseries before (blue) and after MPS filter');   xlabel('time (s)'); 
		subplot(2,1,2);hold on; title('FFT spectrum before (blue) and after MPS filter') 
		xlabel('frequency (s^{-1})'); ylabel('power');
		plot(Freqs,PowerVor);		
		plot(Freqs,PowerNach);
		axis([0 inf 0 1.05*max(PowerVor(Freqs>PowerLineHum*0.9 & Freqs<PowerLineHum+1))]);
		if Plot==2
			figure;hold on; plot(RawData); plot(FilterRep);   title('dataseries / filter')
		end
	end
end


%% =====================================================
%% =====================================================
%% =====================================================


function [Filtered]=filth(CutFreq,MessFreq,Order,Data,f)
%% Version 2 : 24.10.2017 Marc de Lussanet
%% high pass of low pass filter (f= 'high'  'stop' of 'low')
%% der effektive Order ist verdoppelt durch filtfilt

	if f=='h'
		[B,A]=butter(Order,(2*CutFreq)/MessFreq,'high');
	elseif f=='s'
		[B,A]=butter(Order,(2*CutFreq)/MessFreq,'stop');
	else % if f=='l' %% lowpass
		[B,A]=butter(Order,(2*CutFreq)/MessFreq);
	end

	[n,~]=size(Data);
	Filtered = Data;
	for i=1:n
		Filtered(i,:)=filtfilt(B,A,Data(i,:));
	end
end


