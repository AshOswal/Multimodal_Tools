
% simulation makes two time series with high and low beta (13-19Hz and 21-30Hz)
close all;
clear all;
figure(1);
estdelays = [];
% initialisation
fsample   = 1000;                          % frequency
amp1      = 2/sqrt(2);                     % amplitude
time      = 0:1/fsample:2;                 % time points
time(end) = [];
delay1    = 0.04;                          % s (60ms)
delay2    = 0.01;                          % s (10ms)
stdev     = 2;
f1        = 13:1:19;                       % low beta                
f2        = 21:1:30;                       % high beta
f         = {f1;f2};
delay     = [delay1;delay2];
ntrials   = 150;
y1        = 0;
y2        = 0;
SNR       = 1;                            % SNR to set noise
noise     = sqrt(amp1/SNR);


for i = 1:size(f,1)
for n = 1:numel(f{i}(:))
    T = 1/(f{i}(n)); % period
    ranph = 0;       % an additional random phase set to 0
    phase1 = (delay(i)/T)*2*pi;
    y1 = repmat(amp1*sin(time*(f{i}(n))*2*pi+phase1+ranph),ntrials,1) + y1;
    y2 = repmat(amp1*sin(time*(f{i}(n))*2*pi+ranph),ntrials,1) +y2; 
end
end

y1 = y1 + noise*randn(size(y1));
y2 = y2 + noise*randn(size(y2));
figure;plot(mean(y1,1));hold on;plot(mean(y2,1),'r');


% get into ft or spm structure
cfg = [];
cfg.fsample = fsample;
cfg.label = {'1';'2'};
cfg.time = repmat({time},1,ntrials);
for n = 1:ntrials
cfg.trial{n} = [y1(n,:);y2(n,:)];
end
D = spm_eeg_ft2spm(cfg,fullfile(pwd,'cohtest.mat'));
timelock = cfg;


%% coherence analysis ft


cfg = [];
cfg.output ='powandcsd';
cfg.channelcmb = timelock.label;
cfg.keeptrials = 'no';
cfg.keeptapers='no';
cfg.taper = 'dpss';
cfg.method = 'mtmfft';
cfg.foilim = [1 45];
cfg.tapsmofrq = 1;

inp = ft_freqanalysis(cfg, timelock);
phs = unwrap(angle(inp.crsspctrm));
figure(1);
plot(inp.freq,phs);hold on;

ind  = find(inp.freq >= min(f1) & inp.freq <= max(f1));
ff   = inp.freq(ind);
dph  = phs(ind);

ind1 = find(inp.freq >= min(f2) & inp.freq <= max(f2)); 
ff1  = inp.freq(ind1);
dph1 = phs(ind1);

plot(ff,dph,'r');hold on;
plot(ff1,dph1,'k')
[a ,~]     = regress(dph',[ff' ones(size(ff))']);
[a1,~]     = regress(dph1',[ff1' ones(size(ff1))']);

% plot regression line
xlim = min(ff-2):max(ff+2);
plot(xlim,a(1)*xlim+a(2),'r-','LineWidth',2);
xlim1 = min(ff1-2):max(ff1+2);
plot(xlim1,a1(1)*xlim1+a1(2),'k-','LineWidth',2);

estdelay = [a(1)/(2*pi);a1(1)/(2*pi)].*1e3;
tit      = sprintf('SNR = %.2f \nlow beta delay = %.1f ms \nhigh beta delay = %.1f ms',SNR,estdelay(1),estdelay(2));  
h = title(tit);
set(h,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Unwrapped Phase (radians)','FontSize',14);
set(gca,'FontSize',12);



%% coherence analysis spm
%{
S = [];
S.D = D;
S.chancomb =  timelock.label';
S.pretrig  = 0;
S.posttrig = 2000;
S.timewin  = 2000*(D.nsamples-1)/D.fsample;
S.timestep = 1/D.fsample;
S.freqwin  = [1 45];
S.freqres = 1;
S.imaginary = 0; %% added field by AO
S.robust   = 'no';
S.robust.savew = false;
S.robust.bycondition = true;
S.robust.ks = 5;
t = {strcat('spm_robust_','ks',num2str(S.robust.ks))};

[D1 ,Dph] = spm_eeg_ft_multitaper_coherence(S);

figure(1);
subplot(2,1,2);
phs  = unwrap(Dph);
plot(D1.frequencies,phs);hold on;
ind  = find(D1.frequencies >= min(f1) & D1.frequencies <= max(f1));
ff   = D1.frequencies(ind);
dph  = phs(ind);

ind1 = find(D1.frequencies >= min(f2) & D1.frequencies <= max(f2)); 
ff1   = D1.frequencies(ind1);
dph1  = phs(ind1);

plot(ff,dph,'r');hold on;
plot(ff1,dph1,'k');hold on
[a ,~]     = regress(dph',[ff' ones(size(ff))']);
[a1,~]     = regress(dph1',[ff1' ones(size(ff1))']);

% plot regression line
xlim = min(ff-2):max(ff+2);
plot(xlim,a(1)*xlim+a(2),'r-','LineWidth',2);
xlim1 = min(ff1-2):max(ff1+2);
plot(xlim1,a1(1)*xlim1+a1(2),'k-','LineWidth',2);


estdelay = [a(1)/(2*pi);a1(1)/(2*pi)].*1e3;
tit      = sprintf('SNR = %.2f \nlow beta delay = %.1f ms \nhigh beta delay = %.1f ms',SNR,estdelay(1),estdelay(2));  
h = title(tit);
set(h,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',14);
ylabel('Unwrapped Phase (radians)','FontSize',14);
set(gca,'FontSize',12);
%}