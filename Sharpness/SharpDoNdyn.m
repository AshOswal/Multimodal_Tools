function SharpDoNdyn
% This function shows an example in calculating sharpness & DoNgap
% Author: Chien-Hung Yeh
filename1 = 'OFFx.txt'; filename2 = 'ONx.txt';
RAW1 = dlmread(filename1,'	',1,1)';
RAW2 = dlmread(filename2,'	',1,1)';
ch = 1; srate = 1000;
RAW1 = RAW1(ch,srate:end-srate);
RAW2 = RAW2(ch,srate:end-srate);
% 50Hz comb filt
fo = 50; q = 35; bw = (fo/(srate/2))/q;
[b,a] = iircomb(srate/fo,bw,'notch');
FIL1 = filtfilt(b,a,RAW1');
FIL2 = filtfilt(b,a,RAW2');
% 1Hz highpass
FIL1 = ft_preproc_highpassfilter(FIL1',srate,1,6,'but','twopass');
FIL2 = ft_preproc_highpassfilter(FIL2',srate,1,6,'but','twopass');
% wavelet decomp
taxis1 = (0:length(FIL1)-1)/srate;
taxis2 = (0:length(FIL2)-1)/srate;
faxis = 1:1:40;
[COMP1,~] = compute_wavelet(FIL1,srate,taxis1,faxis,'CH');
[COMP2,~] = compute_wavelet(FIL2,srate,taxis2,faxis,'CH');
% welch power
[Fd1,fd1] = pwelch(FIL1',srate,srate/2,srate,srate); Fd1 = Fd1'; fd1 = fd1';
% peak frequency bin
BetaInx = find(fd1>=12 & fd1<=30);
[~,I] = max(Fd1(BetaInx));
BetaPeak = fd1(BetaInx(I'));
% evolving wavelet amp
idx = find(faxis==BetaPeak);
evoAmp1 = sum(COMP1(faxis(idx-2:idx+2),:));
evoAmp2 = sum(COMP2(faxis(idx-2:idx+2),:));
evoAmp1 = smooth(evoAmp1,0.2*srate)';
evoAmp2 = smooth(evoAmp2,0.2*srate)';
for kk=1:floor(length(evoAmp1)/(20*srate))
    evoAmp1((kk-1)*20*srate+1:kk*20*srate) = evoAmp1((kk-1)*20*srate+1:kk*20*srate)-mean(evoAmp1((kk-1)*20*srate+1:kk*20*srate));
end
for kk=1:floor(length(evoAmp2)/(20*srate))
    evoAmp2((kk-1)*20*srate+1:kk*20*srate) = evoAmp2((kk-1)*20*srate+1:kk*20*srate)-mean(evoAmp2((kk-1)*20*srate+1:kk*20*srate));
end
% burst threshold
THRESH = (quantile(evoAmp1,0.75)+quantile(evoAmp2,0.75))/2;
% EMD
rng('default'); IMFs1 = feemd(FIL1,0.1,100,10); IMFb1 = sum(IMFs1(3:4,:));
rng('default'); IMFs2 = feemd(FIL2,0.1,100,10); IMFb2 = sum(IMFs2(3:4,:));
% peaks above threshold
menv1 = max_win(evoAmp1,tukeywin(100,0.9));
mt1 = find(eq(menv1',evoAmp1));
mt1 = mt1(evoAmp1(mt1)>THRESH);
menv2 = max_win(evoAmp2,tukeywin(100,0.9));
mt2 = find(eq(menv2',evoAmp2));
mt2 = mt2(evoAmp2(mt2)>THRESH);
i1 = sort(evoAmp1(mt1)); idx1 = find(evoAmp1==i1(end-4));
i2 = sort(evoAmp2(mt2)); idx2 = find(evoAmp2==i2(end-4));
% within burst sharpness & DoNgap
[SHARP1,gapDoN1] = compute_withinBurst(IMFb1,srate,idx1);
[SHARP2,gapDoN2] = compute_withinBurst(IMFb2,srate,idx2);
fprintf('under OFF med...\n')
fprintf(['sharpness(1st) = ' num2str(SHARP1(1)) '\n'])
fprintf(['sharpness(2nd) = ' num2str(SHARP1(2)) '\n'])
fprintf(['DoNgap(1st) = ' num2str(gapDoN1(1)) '\n'])
fprintf(['DoNgap(2nd) = ' num2str(gapDoN1(2)) '\n'])
fprintf('under ON med...\n')
fprintf(['sharpness(1st) = ' num2str(SHARP2(1)) '\n'])
fprintf(['sharpness(2nd) = ' num2str(SHARP2(2)) '\n'])
fprintf(['DoNgap(1st) = ' num2str(gapDoN2(1)) '\n'])
fprintf(['DoNgap(2nd) = ' num2str(gapDoN2(2)) '\n'])
end
function [wav,wav_complex] = compute_wavelet(bipolar,srate,taxis,faxis,STR_ch)
[DataFT,Cfg] = preSpec(bipolar,taxis,faxis,srate,STR_ch);
DataSpec = ft_freqanalysis(Cfg,DataFT);
wav_complex = squeeze(DataSpec.fourierspctrm); % wavelet complex
wav = abs(wav_complex); % wavelets envelope
for ii=1:size(wav,1)
    getNan = isnan(wav(ii,:));
    wav(ii,getNan) = 0;
end
repwav = repmat(mean(wav,2),1,size(wav,2));
wav = ((wav-repwav)./repwav)*100; % t-score
wav = imgaussfilt(wav,2); % gaussian filter
end
function [Data,Cfg] = preSpec(data,taxis,faxis,samplerate,channel)
Cfg.toi = taxis;
Cfg.foi = faxis;
Cfg.method = 'wavelet';
Cfg.output = 'fourier';
Cfg.channelcmb = {'all' 'all'};
Cfg.keeptrials ='yes';
Cfg.width = 10;
Cfg.gwidth = 5;
Data.fsample = samplerate;
Data.label = {channel};
Data.time = {taxis};
Data.trial = {data};
end
function [SHARP,gapDoN] = compute_withinBurst(X,srate,idx)
pt1 = idx-0.2*srate:idx; X1 = X(pt1);
pt2 = idx:idx+0.2*srate; X2 = X(pt2);
menv11 = max_win(X1,tukeywin(50,0.8)); mt11 = find(eq(menv11',X1));
menv12 = max_win(-X1,tukeywin(50,0.8)); mt12 = find(eq(menv12',-X1));
gap1=[]; mt1 = union(mt11,mt12);
for ii=2:length(mt1)-1
    gap1 = [gap1 mt1(ii)-srate*5/1000:mt1(ii)+srate*5/1000];
end
gap1(gap1<=0 | gap1>length(X1))=[];
menv21 = max_win(X2,tukeywin(50,0.8)); mt21 = find(eq(menv21',X2));
menv22 = max_win(-X2,tukeywin(50,0.8)); mt22 = find(eq(menv22',-X2));
gap2=[]; mt2 = union(mt21,mt22);
for ii=2:length(mt2)-1
    gap2 = [gap2 mt2(ii)-srate*5/1000:mt2(ii)+srate*5/1000];
end
gap2(gap2<=0 | gap2>length(X2))=[];
% gap-DoN
tmp1 = gapcdn(X1',1/srate,gap1); gapDN1 = tmp1(end);
tmp2 = gapcdn(X2',1/srate,gap2); gapDN2 = tmp2(end);
gapDoN = [gapDN1 gapDN2];
% sharpness
[peaks1,troughs1] = findpt(X1,X1); Es1 = union(peaks1,troughs1);
[peaks2,troughs2] = findpt(X2,X2); Es2 = union(peaks2,troughs2);
widthS = srate*5/1000;
ampPC = 0; amps1 = X1; amps2 = X2;
Es1(Es1<=5 | Es1>=length(X1)-5)=[];
Es2(Es2<=5 | Es2>=length(X2)-5)=[];
SHARP1 = quantile(Esharp(X1,Es1,widthS,ampPC,amps1),0.75);
SHARP2 = quantile(Esharp(X2,Es2,widthS,ampPC,amps2),0.75);
SHARP = [SHARP1 SHARP2];
end