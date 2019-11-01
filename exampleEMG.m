% Script for testing Humfilering on simulated data
%
% Data are white noise; Hum can be any kind of repeated signal.
%
% Marc de Lussanet, Movement Science, University of Muenster
% 27.10.2019

	%% constants
	Fs     =1000;   
	Frq    =50; 
	Period =Fs/Frq;
	Periods=5000;
	Time   =(0:Period*Periods)/Fs;	
	ErrPlot= 1; 
	Nloop  = 10; %100; %500; 	
	Win    = 50;  % Window for MPS filter (in periods of hum; 0=infinite)

	%% flags
	RemoveHum    = 1;
	AddStep      = 1;
	HumAmpl      = 10;
	%LoopHFsignal= [1 0.5 0.2 0.1 0.05 0]; % 1=white noise data; 0=60Hz low pass
	%LoopHFsignal= 0.1; % 1=white noise data; 0=60Hz low pass
	LoopHFsignal = 1; % 1=white noise data; 0=60Hz low pass
		 
	% time windows of interest
	TRemoveHum= 20000: 20000+2000;
	TStep     = 1    : 10000;		% 
	Win0      = 1    : 200; 
	Win2      = TStep(end)   + (-100 : 100); %2450:2550; 
	Win3      = TRemoveHum(1)+ (-100 : 100); %3450:3550; 
	WinStabl  = 25000:   length(Time)-1000;  %WinStabl=500:2000;
	
	%% kind of hum
	HumTypes = {'SinHum','WaveHum','RndHum','PeakHum'};
	HumType  = 2;
	
	% Create different kinds of HUM
	if     HumType == 1 	% 1. sinusoidal hum 
		Hum    = cos(2*pi*(Frq)*Time)*1;
	elseif HumType == 2 	% 2. random-like repeats
		Repeat = [10 5 0 2 -4 -9 -10 -6 8 9 0 0 2 -9 10 -5 -4 -3 8 9]/10; % random-like
	elseif HumType == 3 	% 3. wave with harmonics
		Repeat = [0 2 4 6 8 10 8 6 4 2 0 -2 -4 -6 -8 -10 -8 -6 -4 -2]/10; % wave with harmonics
	elseif HumType == 4 	% 4. 50Hz pulse
		Repeat = [10 -10 -10 -10 -10 -10 -10 -10 -10 -10  -10 -10 -10 -10 -10 -10 -10 -10 -10 -10]/10; % pulse
	end
	if HumType>1
		Hum = repmat(Repeat,1,Periods+1) - mean(Repeat); Hum(length(Time)+1:end) = [];
	end
	Hum = Hum * HumAmpl;
	
	%% TEST : sudden removal of hum
	if RemoveHum,  Hum(TRemoveHum) = 0.8 * Hum(TRemoveHum); end

	CmErrMPSs_HF = [];
	
	for HF = LoopHFsignal
		CmErrMPS0 = [];
		CmErrMPSs = [];
		CmErrMPS2 = [];
		CmErrMPS3 = [];
		for i=1:Nloop	
			%% Create white noise data
			Data = randn(1,length(Time)); % *sqrt(10);
			Data = HF * Data + (1-HF) * filth(60,1000,1,Data,'l'); % niet zo veel hoogfrequent

			%% TEST : add a step
			if AddStep,  Data(TStep) = Data(TStep) + 100; end

			%% Simulated raw data
			DataHum = Hum + Data;

			%% apply MPS filter
			DataMedPhase = removeHumByPMS(DataHum, round(Fs/Frq), Win, 0);

			ErrMPS  = abs(DataMedPhase-Data);
			ErrMPSs = filth(20,1000,2,abs(ErrMPS),'l');
			CmErrMPSs = [CmErrMPSs ErrMPS(WinStabl)]; %#ok<*AGROW>
			CmErrMPS0 = [CmErrMPS0 ErrMPS(Win0)];
			CmErrMPS2 = [CmErrMPS2 ErrMPS(Win2)];
			CmErrMPS3 = [CmErrMPS3 ErrMPS(Win3)];
			if ErrPlot && i<10 && HF == LoopHFsignal(end)
				if i==1,	H = figure; hold on; end
				figure(H)
				plot(Time(20:end), ErrMPSs(20:end), 'r');
			end
		end
		CmErrMPSs_HF = [CmErrMPSs_HF CmErrMPSs];
	end

	if length(LoopHFsignal) > 1
		J      = figure; hold on; 
		Labels = {};
		for i=1:length(LoopHFsignal)
			Labels = [Labels, repmat({sprintf('%.2f',LoopHFsignal(i))},1,length(CmErrMPSs))];
		end
		boxplot(CmErrMPSs_HF,Labels,'PlotStyle','compact','symbol',''); title('PMS'); grid on;
		ylabel('absolute error')
		ylim([0 0.5]);
	else
	
		I      = figure; hold on; 
		PlMPS  = [CmErrMPS0, CmErrMPSs, CmErrMPS2, CmErrMPS3];
		Labels = [
			repmat({'Onset'},    1,length(CmErrMPS0)), ...
			repmat({'Stable'},   1,length(CmErrMPSs)), ...
			repmat({'Transient'},1,length(CmErrMPS2)), ...
			repmat({'Noise Off'},1,length(CmErrMPS3))];
		ylabel('absolute error')
		boxplot(PlMPS,Labels,'PlotStyle','compact','symbol',''); title('PMS'); grid on;
		ylabel('absolute error')
		ylim([0 1.4]);

		if ErrPlot
			figure(H);
			WindowSec = Win * Frq / Fs;
			xlabel('time (s)'); ylabel('abs. error (low-pass filtered)');
			ylim([0 1]);
		end

		fprintf('HumAmpl = %.1f. Abs. error in stable period: %.3f %.3f / %.3f (%.3f) [smoothed abs(median, mean / max); (sd)]\n', HumAmpl, ... 
			median(CmErrMPSs), ...
			mean(CmErrMPSs), ...
			max( CmErrMPSs), ...
			std( CmErrMPSs ));

		figure; Win2sm = TStep(end)-90 : TStep(end)+10;
		subplot(1,2,1); hold on; plot(Time(Win2),Data(Win2),'g','LineWidth',2); 
										 plot(Time(Win2),DataMedPhase(Win2),'r');
		title('data with step'); xlabel('time (s)'); ylabel('data');xlim([9.9 10.1]);
		subplot(1,2,2); hold on; plot(Time(Win2sm),Data(Win2sm),'g','LineWidth',2); 
										 plot(Time(Win2sm),DataMedPhase(Win2sm),'r');
		title('data at step'); xlabel('time (s)'); ylabel('data'); legend('white noise','PMS-filtered','Location','southwest');
		xlim([9.97 10.002]);ylim([90 105]);hold off; 

		S = figure; 
		S(1)=subplot(1,3,1); hold on; plot(Time(Win0),abs(DataMedPhase(Win0)-Data(Win0))); title('onset')
		xlabel('time (s)'); ylabel('error'); ylim([0 2.5]);
		S(2)=subplot(1,3,2); hold on; plot(Time(Win2),abs(DataMedPhase(Win2)-Data(Win2))); title('transient')
		S(3)=subplot(1,3,3); hold on; plot(Time(Win3),abs(DataMedPhase(Win3)-Data(Win3))); title('noise off')
		linkaxes(S,'y');
	end
	
	
	
	
	
	
%%  =======================================================================
%%  =======================================================================
%%  =======================================================================
%%  =======================================================================
%%  =======================================================================





