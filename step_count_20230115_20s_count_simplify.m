function [] = step_count_20230115_20s_count()
% 改自 step_count_20221027
% acc window 5.5s => 3s
% specify_global_constants_2016_0729 => specify_global_constants_2022_1013
% filter_signal_33 => filter_signal_2022_1013
% stream_acc2Eul_0709 => stream_acc2Eul_2022_0805
% copy_ghost_bins_0707 => copy_ghost_bins_2022_1013
% move A_min to the front of autocorr 2023_0115
% *************************************************************************
global DEBUG_MODE;
DEBUG_MODE = 0;
icase = 3; %*************** HERE 1=> cf fitbit and truth; 2=> cf to truth; 3: run  one exp only.
sintest = 0;
tbegin = 0.02;
tend = 0;
fs0 = 50; % data sampe rate
% *************************************************************************

[GlbCnt] = specify_global_constants_2022_1013(fs0);

if (icase == 3)
    file_index_min = 102; % 6 4;  % 102; % **************** HERE
    file_index_max = file_index_min; % default 200
elseif (icase == 1 || icase==2)     
    file_index_min = 1; file_index_max = 200; % default 200  
end
tstart = tic; % DEBUG ONLY

global_step = 0;     
for ff = file_index_min:file_index_max   
   [ax, ay, az, fileName] = input_accelerometer_data4_new_format(GlbCnt, icase, ff); %SINSIN  BUG in fs; 0503/2016
   if (tend ~= 0)
       [ax, ay, az] = clip_signal(ax, ay, az, GlbCnt.Fs, tbegin, tend, sintest);
   end
   if (isempty(ax))
       continue;
   end
   fprintf('ff=%d filename=%s\n',ff,fileName);                 
   n = length(ax);
   
   [file_step] = effEuler_step_count_main_0709(GlbCnt, n, ax, ay, az);
   
   global_step = global_step + file_step;
   
   fprintf(' ########################  End of Program file\n Id=%d Total steps  = %d\n', ff, global_step);  
end % ff read data loop
time_elapsed = toc(tstart);
fprintf(' Elapse Time = %d\n',time_elapsed);
return; % online_step_count_test_184_new_format 

function [ax, ay, az] = clip_signal(ax, ay, az, Fs, tbegin, tend, sintest)
global DEBUG_MODE;
if (sintest)
    freq = 1.6; % Hz
    walk_time = tend - tbegin;
    nonwalk_time = 2; % sec
    amp = 1000; % mg
    truth = walk_time * freq; % step
    fprintf("truth step : %d\n", truth);
    
    tw = 0:1/Fs:walk_time;
    xw = amp * sin(2*pi*freq*tw);
    
    tn = 0:1/Fs:(nonwalk_time/2) - 1/Fs;
    xn = zeros(1, length(tn));
    
    x = [xn xw xn];
    
    ax = x + 1000;
    ay = zeros(1, length(ax));
    az = ay;
else
    ns = tbegin * Fs + 1;
    ne = tend * Fs + 1;
    ax = [zeros(1,Fs)+1000 ax(ns:ne) zeros(1,Fs)+1000];
    ay = [zeros(1,Fs) ay(ns:ne) zeros(1,Fs)];
    az = [zeros(1,Fs) az(ns:ne) zeros(1,Fs)];
end

if (DEBUG_MODE)
    a = sqrt(ax.^2 + ay.^2 + az.^2);
    figure(11);clf;
    plot(a);
end
return

function [Fv,Y] = fft_freq_domain(Fs,x)
L = length(x);
Fn = Fs/2;
FTa = fft(x)/L;
Fv = linspace(0,1,fix(L/2)+1)*Fn;
Iv = 1:length(Fv);
Y = abs(FTa(Iv))*2;
return

function [acf] = autocorr_matlab_edit_brute_force(y, numLags)
y = y - mean(y); % 2022_1026
acf = nan(numLags+1,1);
    for j = 0:numLags
        cross   = y(1:end-j).*y(j+1:end);
        iNonNaN = ~isnan(cross);
        if any(iNonNaN)
            T        = sum(iNonNaN)+sum(~isnan(y(1:j)));
            acf(j+1) = nansum(cross)/T;
        end
    end
acf = acf./acf(1);
acf = acf';
return

function [acf] = autocorr_brute_force_CY(x, tau, n, dinc)
x = x - mean(x);
tau = floor(tau / dinc);

j = 0;
for (i = 0 : tau)
    sum = 0;
    for (k = 0 : dinc : n - j - 1)
        sum = sum + x(k + 1) * x(j + k + 1);
    end
    acf(i + 1) = sum;
%     assert(abs(acf(i + 1)) < 2^31);
    j = j + dinc;
end
% acf = acf./acf(1);
return

function [global_step, tt] = effEuler_step_count_main_0709(GlbCnt, tbreak, ax, ay, az)
global compBgnPt compEndPt  Sstamp

% allocate memory
accWrk = zeros(1,GlbCnt.Phys_win_size + GlbCnt.Ovlp);
acc = zeros(1,GlbCnt.Phys_win_size +  GlbCnt.Ovlp);

if (GlbCnt.FilterType == 1)
    nzeros = GlbCnt.Fs; % length of filter zeros padding
    acc_iir = zeros(1, nzeros + GlbCnt.Phys_win_size + GlbCnt.Ovlp);
end

[Eul] = allocate_Euler_struct_2(GlbCnt.MaxCompBin2);

global_step = 0;
if (tbreak < GlbCnt.Di) % BUG, add this condition
    return;
end
global_step = 0;

tt = 0; ic = 0; ib = 0; n = 0; bin = 0; cwin = 1;
Sstamp(1).bgn = 0; Sstamp(1).end = 0; Sstamp(1).step = 0; % BBB
maxCompBin = GlbCnt.MaxCompBin1;% maxCompBin = 20;%for1
firstCal = 1; % III

while (1) % each online window
    tt = tt + 1; % global 
    ib = ib + 1; % relative to bin
    ic = ic + 1; % point relative to comp domain
    
    accWrk(ic) = sqrt(ax(tt).^2 + ay(tt).^2 + az(tt).^2); % cstyle
    
    if (ib == GlbCnt.Di || tt==tbreak) % extract bin
        bin = bin + 1;        
        ib = 0;
    end % bin
    
    if (bin == maxCompBin || tt==tbreak)
        
        if (GlbCnt.FilterType == 1)
            acc_iir(1:nzeros) = zeros(1, nzeros) + mean(accWrk);
            acc_iir(1 + nzeros:GlbCnt.Phys_win_size + GlbCnt.Ovlp + nzeros) = accWrk;
            acc_iir = filter_signal_2022_1013(acc_iir, GlbCnt);
            acc = acc_iir(1 + nzeros:GlbCnt.Phys_win_size + GlbCnt.Ovlp + nzeros);
        else
            acc = accWrk;
            acc = filter_signal_2022_1013(acc, GlbCnt);
        end
%         [Fv,Y] = fft_freq_domain(GlbCnt.Fs,acc);
%         figure(1);clf
%         plot(Fv, Y);
%         [FvWrk,YWrk] = fft_freq_domain(GlbCnt.Fs,accWrk);
%         figure(2);clf
%         plot(FvWrk, YWrk);
%         ylim([0 max(Y)]);
%         [1];
        [Eul] = stream_acc2Eul_2023_0115(GlbCnt, bin, acc, Eul, firstCal);
        
        if (firstCal == 1) % III
            firstCal = 0;
        end
        
        Eul.nbin = bin;
        phyBgnBin = 1 + GlbCnt.OvlpDi;% phyBgnBin = 1 + 3;%for1
        if (GlbCnt.Di == 25)
            phyEndBin = Eul.nbin - GlbCnt.OvlpDi;
            compBgnPt = (cwin-1)*GlbCnt.Phys_win_size + 1 - 2*GlbCnt.Ovlp; % in global coord
        elseif (GlbCnt.Di == 50)
            phyEndBin = Eul.nbin - 2;
            compBgnPt = (cwin-1)*GlbCnt.Phys_win_size + 1 - (GlbCnt.Ovlp + GlbCnt.Ncomp);
        end
        compEndPt = (cwin)*GlbCnt.Phys_win_size;
             
        bgnBin = 1; 
        endBin = Eul.nbin;
         
        if (cwin==1)
            phyBgnBin = 1;
            compBgnPt = 1;     
        end
        if (tt==tbreak)
            phyEndBin =  Eul.nbin;
            
            if (cwin > 1)
                compEndPt = compBgnPt+GlbCnt.Ovlp+ic-1;
            else
                compEndPt = ic;
            end
            
        end
        
        [winStep] = eff_euler_online_step_count_1(bin, Eul, GlbCnt, phyBgnBin, phyEndBin, bgnBin, endBin);
        
        fprintf('cwin = %d  nbin = %d  maxCompBin=%d  winStep=%d\n', cwin, endBin, maxCompBin, winStep);
        global_step = global_step + winStep;
        assert(global_step>=0); % NOTES it is possible that winStep < 0; 因為除了end bin 之外,我們刪除的是前一個window的step
        if (tt == tbreak)  % DDD
            break;
        end
% ************************************************************ 2022_1013_DY
        [accWrk, Eul] = copy_ghost_bins_2022_1127(GlbCnt.Ovlp, GlbCnt.OvlpDi, GlbCnt.OvlpBin, ic, accWrk, Eul, maxCompBin);% [acc ,Eul] = copy_ghost_bins_2022_0805(100, 3, n, acc, Eul, maxCompBin);
        bin = GlbCnt.OvlpBin; % bin = 5 sec
        n = GlbCnt.Ovlp; % length(acc) - BinSize
        maxCompBin = GlbCnt.MaxCompBin2; % maxCompBin = 25;%for1
        cwin = cwin + 1;
        ic = GlbCnt.Ovlp; % ic = 100;%for1
% ******************************************************** 2022_1013_DY_end
    end
end % while (1) t loop nwrk ib
fprintf('global_step = %d\n',global_step);
return % effEuler_step_count_main_0709 

function [Eul] = stream_acc2Eul_2023_0115(GlbCnt, bin, acc, Eul, firstCal)

if (GlbCnt.Di == 25)
    if (firstCal) % III
        bgnBin_acc = 1;
    else
        bgnBin_acc = 2*GlbCnt.OvlpDi + 1;
    end
elseif (GlbCnt.Di == 50)
    if (firstCal) % III
        bgnBin_acc = 1;
    else
        bgnBin_acc = 2*GlbCnt.OvlpDi;
    end
end

if (firstCal)
    acc_bgn = 1;
    bgnBin_acf = 1;
else
    acc_bgn = GlbCnt.Ovlp + 1;
    bgnBin_acf = GlbCnt.OvlpDi + 1;
end

[Eul] = acc_max_and_min(GlbCnt, bin, bgnBin_acc, acc, Eul, acc_bgn);
[Eul] = cal_acc_mean_1(GlbCnt, bin, Eul); % calculate acc_mean
[Eul] = process_similarity_and_pitch(GlbCnt, bin, acc, Eul, bgnBin_acf);
return; % stream_acc2Eul_2023_0115

function [Eul] = acc_max_and_min(GlbCnt, nbin, bgnBin, acc, Eul, acc_bgn)
% 2023_0115_CY
n = length(acc);
for (bin = bgnBin:nbin) % cal max_acc & min_acc in each bin
    acc_end = min(acc_bgn + GlbCnt.Di - 1, n);
    assert(acc_end <= n);
    Eul.max_acc(bin) = max(acc(acc_bgn:acc_end));
    Eul.min_acc(bin) = min(acc(acc_bgn:acc_end));
    acc_bgn = acc_bgn + GlbCnt.Di;
end % bin 
return % acc_max_and_min

function [Eul] = cal_acc_mean_1(GlbCnt, nbin, Eul)
% ng2 = max(floor(GlbCnt.Win_acf/GlbCnt.Di) - 2,0); % number of euler in 2.5 (3(win_acf)-0.5(di)) sec; 2016-0614 fix bug
ng2 = max(floor(GlbCnt.Win_acfDi) - 2,0); %   number of euler in 2.5 (3(win_acf)-0.5(di)) sec; 2016-0614 fix bug
for (jj = 1:nbin) % cal acc_mean   
    % fprintf('2016-0707 I am here BUG for nbin !!! !!! 48->36(cal vs exact\n')
    jj2 = min(jj + ng2,nbin); % is this a bug 0411-2016 ??
    Eul.acc_mean(jj) = max(Eul.max_acc(jj:jj2)) - min(Eul.min_acc(jj:jj2));
end
return; % cal_acc_mean_1

function [Eul] = process_similarity_and_pitch(GlbCnt, nbin, acc, Eul, bgnBin_acf)
% 2023_0115_CY
max_lag = GlbCnt.Max_lag;
acc_bgn = 1;
acc_end = GlbCnt.Ncomp;
for (bin = bgnBin_acf:nbin)
    if (Eul.acc_mean(bin) < GlbCnt.Amin)
%         fprintf("123456789\n");
        acf_val_cur = 0;
        p1 = 4 * GlbCnt.Fs;
        p2 = p1;
    else
        data = acc(acc_bgn:acc_end);
        max_lag1 = min(length(data)-1, max_lag);
        [acf_val_cur, p1, p2] = fast_cal_similarity_0717(GlbCnt.Dinc, GlbCnt, max_lag1, data); % CCC yhw0714
    end
    Eul.acf_val(bin) = acf_val_cur; % value  yhw: can be next peak ??
    Eul.pitch1(bin) = p1;
    Eul.pitch2(bin) = p2; 
    acc_bgn = acc_bgn + GlbCnt.Di;
    acc_end = min(acc_end + GlbCnt.Di, length(acc));
end
return % process_similarity_and_pitch

function [acc, Eul] = copy_ghost_bins_2022_1127(Ovlp, OvlpDi, OvlpBin, nmax, acc, Eul, maxCompBin)
% 沼澤之灰色地帶長度 = Ovlp;  copy acc
acc(1:Ovlp) = acc(nmax - Ovlp + 1:nmax);
acc(Ovlp + 1:nmax) = zeros(1, nmax - Ovlp);
if (OvlpBin == 10)
    Eul.max_acc(1:OvlpBin) = Eul.max_acc(maxCompBin-OvlpBin + 1:maxCompBin);
    Eul.min_acc(1:OvlpBin) = Eul.min_acc(maxCompBin-OvlpBin + 1:maxCompBin);
    Eul.acf_val(1:OvlpDi) = Eul.acf_val(maxCompBin-OvlpBin + 1:maxCompBin - OvlpDi);
    Eul.pitch1(1:OvlpDi) = Eul.pitch1(maxCompBin-OvlpBin + 1:maxCompBin - OvlpDi);
    Eul.pitch2(1:OvlpDi) = Eul.pitch2(maxCompBin-OvlpBin + 1:maxCompBin - OvlpDi);
elseif (OvlpBin == 5)
    Eul.max_acc(1:OvlpBin) = Eul.max_acc(maxCompBin-OvlpBin + 1:maxCompBin);
    Eul.min_acc(1:OvlpBin) = Eul.min_acc(maxCompBin-OvlpBin + 1:maxCompBin);
    Eul.acf_val(1:OvlpDi) = Eul.acf_val(maxCompBin-OvlpBin + 1:maxCompBin - OvlpDi+1);
    Eul.pitch1(1:OvlpDi) = Eul.pitch1(maxCompBin-OvlpBin + 1:maxCompBin - OvlpDi+1);
    Eul.pitch2(1:OvlpDi) = Eul.pitch2(maxCompBin-OvlpBin + 1:maxCompBin - OvlpDi+1);
end
return; % copy_ghost_bins_2022_1013

function [winStep] = eff_euler_online_step_count_1(nbin, Eul0, GlbCnt, phyBgnBin, phyEndBin, bgnBin, endBin)
% Note if decrease win_acck => accuracy get worse, bit fp get better a
% little e.g. exp0331 rr=20

[Eul0] = modify_sim_and_pitch_12(GlbCnt, Eul0);  

% transfer form physical coord to computational boundary;
sbgn = bgnBin;
send = endBin;

pitch1 = Eul0.pitch1(sbgn:send);
pitch2 = Eul0.pitch2(sbgn:send);
acf_val = Eul0.acf_val(sbgn:send);
acc_mean = Eul0.acc_mean(sbgn:send);
wd = Eul0.wd(sbgn:send);
wdWrk1 = Eul0.pitch1(sbgn:send);
LeftBnd = Eul0.LeftBnd(sbgn:send);
RightBnd = Eul0.RightBnd(sbgn:send);

Eul0.nbin = send-sbgn+1;
assert(Eul0.nbin  == length(wd));
      
[wd, pitch1, pitch2, acf_val, numBnd, LeftBnd, RightBnd] = set_motion_status_6_2(GlbCnt, nbin, acc_mean, pitch1, pitch2, acf_val, wd, wdWrk1, LeftBnd, RightBnd);
  
[winStep] = process_each_window_Euler_only_0625(GlbCnt, numBnd, LeftBnd, RightBnd, pitch1, pitch2, phyBgnBin, phyEndBin);

assert(isempty(winStep)==0);

return; % eff_euler_online_step_count_1 pbgnPt pEndPt  

function [Eul] = allocate_Euler_struct_2(maxBin)
Eul.maxBin = maxBin; % maximum possible bin #
Eul.nbin = 0;
Eul.max_acc = zeros(1,Eul.maxBin); 
Eul.min_acc = zeros(1,Eul.maxBin);  
Eul.acc_mean = zeros(1,Eul.maxBin); 
Eul.acf_val = zeros(1,Eul.maxBin); % 2 byte ?

Eul.pitch1 = zeros(1,Eul.maxBin); % 1 or 2byte ?
Eul.pitch2 = zeros(1,Eul.maxBin);

Eul.wd = zeros(1,Eul.maxBin); % one byte
Eul.wdWrk1 = zeros(1,Eul.maxBin); 

Eul.LeftBnd = zeros(1,Eul.maxBin);
Eul.RightBnd = zeros(1,Eul.maxBin);
return; % allocate_Euler_struct_2

function [acc_filtered] = filter_signal_2022_1013(acc0, GlbCnt)

if (GlbCnt.FilterType)
    if (GlbCnt.FilterType == 1)
        acc_filtered = filter(GlbCnt.CB, GlbCnt.CA, acc0);
    elseif (GlbCnt.FilterType == 2)
        maxSift = 1;maxPhase0 = 4;noiseLevel = 100;
        imf = upemd(acc0(:)', maxSift, maxPhase0, noiseLevel, 50);
        ni1 = 2;ni2 = 6;
        ni2 = min(ni2,size(imf,1)); % 0306-2015 FIX bug
        if (ni1==ni2) % 0403-2015 fix it
            imf_high = (imf(ni1:ni2,:));
        else
            imf_high = sum(imf(ni1:ni2,:));
        end 
        acc_filtered = imf_high(:)';
    end
else
    acc_filtered = acc0 - mean(acc0);
end
return; % filter_signal_2022_1013

function [winStep] = process_each_window_Euler_only_0625(GlbCnt, numSeg, LeftBnd, RightBnd, pitch1, pitch2, phyBgnBin, phyEndBin)
global Glb_cur_tbgn  Glb_cur_tend  Glb_prev_tbgn  Glb_prev_tend Sstamp
global compBgnPt
% a window may consist of many segments, count step (=segStep) in each segment and then sum into steps in a window (winStep)

winStep = 0; 
if (numSeg == 0) % AAA
    Glb_prev_tbgn = compBgnPt + (phyBgnBin-1)*(GlbCnt.Di);
    Glb_prev_tend = compBgnPt + (phyEndBin)*(GlbCnt.Di);
    return;
end % AAA

for (s=1:numSeg)
    jbgn = LeftBnd(s); jend = RightBnd(s);
    %  serch for peak in between [dibgn diend) or [jbgn jend)
    
    binBgn = max(phyBgnBin,jbgn);
    binEnd = min(phyEndBin,jend); % SINSIN phyBgn YHW0707 seg2 = 36, but exact =38

    [segStep] = segment_step_by_integration_0609_2016(GlbCnt, binBgn, binEnd, pitch1, pitch2); % physical domain 的 wd(bgn) wd(end)=?
    
    winStep = winStep + segStep; 

    Glb_cur_tbgn = compBgnPt + (binBgn-1)*(GlbCnt.Di);
    Glb_cur_tend = compBgnPt + binEnd*(GlbCnt.Di) -1;
   
    Sstamp(2).step = segStep; % assign 
    Sstamp(2).bgn = Glb_cur_tbgn;
    Sstamp(2).end = Glb_cur_tend;
    [del_step, Sstamp] = dump_small_step_0713(GlbCnt,Sstamp);
    %************簡化**************
    del_step=0;
    %*******************************
    winStep = winStep - del_step;
    
    Glb_prev_tbgn = Glb_cur_tbgn; 
    Glb_prev_tend = Glb_cur_tend;     
end % loop segment
assert(isempty(winStep)==0);
return % process_each_window_Euler_only_0625     

function [wd, pitch1, pitch2, acf_val] = get_initial_walk_status_3(GlbCnt, nbin, acc_mean, pitch1, pitch2, acf_val, wd)

for (i=1:nbin)
    ip1 = 0;
    if (i >=nbin-2)
        ip1 = max(nbin-3,1);
    end
    i1 = max(1,i-1);
    if (10*acc_mean(i) < 5* acc_mean(i1))
        ip1 = i1;
    end
    i2 = max(1,i-2);
    if (10*acc_mean(i) < 3*acc_mean(i2) && acc_mean(i1) < acc_mean(i2)) % NEW
        ip1 = i2;
    end
    if (ip1 > 0)
        ip2 = max(1,i-(GlbCnt.Win_acfDi)+1); % A1
        [vmax, idx] = max(acf_val(ip2:ip1));
        ip2 = idx + ip2 - 1;
        for (j=i:-1:ip2) % fix long scale bias due to ABD
            pitch1(j) = pitch1(ip2);
            pitch2(j) = pitch2(ip2);
            acf_val(j) = acf_val(ip2);
        end
    end % if
end % i
% ************************* PART 2: main parts ***************************  
[wd, pitch1, pitch2, acf_val] = get_initial_walk_status_part_II(GlbCnt, 1, nbin, 1, 1, acc_mean, pitch1, pitch2, acf_val, wd); % 0428

wdWrk1 = wd;
span = 2;
for (i=1:nbin)  % KEPP_SMOOTH_PITCH  
    pp = pitch2(i);
    i1 = max(i-span,1); i2 = min(i+span,nbin);
    c1 = 0; c2 = 0;
    thresh1 = 2;
    for (j=i1:i)
        if (10*abs(pitch2(j) - pp) < thresh1*pp)
            c1 = c1 + 1;
        end % if
    end % i
    for (j=i:i2)
        if (10*abs(pitch2(j) - pp) < 2*thresh1*pp)
            c2 = c2 + 1;
        end % if
    end % i
    cond1 =  (c1 + c2 >= 3);
    if (cond1 == 0 && wd(i) > 1)
        wdWrk1(i) = 0;    
    end 
end % i
wd = wdWrk1;
return; % get_initial_walk_status_3

function [wd, pitch1, pitch2, acf_val] = get_initial_walk_status_part_II(GlbCnt, j1, j2, icase, floating_value, acc_mean, pitch1, pitch2, acf_val, wd)
% ************************* PART 2: main parts *************************** 
for (jj=j1:j2)
    st_pitch = 2*( pitch2(jj) < GlbCnt.Phigh && pitch2(jj) > GlbCnt.Plow );   
    pp = pitch2(jj); pp = max(pp,1); pp = min(pp,ceil(1.2*GlbCnt.Fs)); % FFF notes 2016-0729; SINSIN 1.2 sec , into global parameter !! 

    if (icase==1) % deafult
        gmin = GlbCnt.GminTable(pp)/GlbCnt.FacAmp + GlbCnt.Amin;
    else %  % exception treasmill, for exp0115,r=36; violate g=g(f) value; : add few  secs at tail
        gmin = (GlbCnt.GminTable(pp)/GlbCnt.FacAmp)/2 + GlbCnt.Amin; % 
    end
    %***********簡化*************
    gmin = GlbCnt.Amin;
    %****************************
    st_acc = 2*(acc_mean(jj) > gmin); % 1228/2014
    
    if (2*acf_val(jj) > GlbCnt.Fac_acf) % change to 0.55 ~ 0.6 ?
        st_simi = 2;
    elseif (10*acf_val(jj) < GlbCnt.Fac_acf*4)
        st_simi = 0;
    else
      if (GlbCnt.FacAmp*acc_mean(jj) > 2000) %  if acceleration is large, relax acf criteria, it is reasonable  physically
           st_simi = 2;
      else
           % st_simi = floating_value; % in between 0427-2016 
           st_simi = 1; % in between 0427-2016
      end
    end
    
    if ((st_acc == 0) || (st_pitch == 0)  || (st_simi == 0) )
        wd(jj) = 0; % non-walk
    elseif (st_acc == 2 && st_pitch==2 && st_simi==2)
        wd(jj) = 2; % walk
    else
        wd(jj) = floating_value; % to be determined
    end
end % jj
return; % get_initial_walk_status_part_II

function [wd, pitch1, pitch2, acf_val, numBnd, LeftBnd, RightBnd] = set_motion_status_6_2(GlbCnt, nbin, acc_mean, pitch1, pitch2, acf_val, wd, wdWrk1, LeftBnd, RightBnd)

% set motion status (flags) based on acc_mean, similarity coefficient and
% pitch; 
% wd == walk status (should use ws instead of wd, SIN)
[pitch1, pitch2] =  modify_irratic_pitch_7(GlbCnt, pitch1, pitch2); % yhw0405 move here
[wd, pitch1, pitch2, acf_val] = get_initial_walk_status_3(GlbCnt, nbin, acc_mean, pitch1, pitch2, acf_val, wd);   
[numBnd, LeftBnd, RightBnd, wd] = trim_walk_status_8(GlbCnt.Sec5, GlbCnt,nbin, acc_mean, wd, LeftBnd, RightBnd, wdWrk1, acf_val, pitch1, pitch2);

return; %set_motion_status_6_2 acc wdWrk1

function [numBnd, LeftBnd, RightBnd, wd] = trim_walk_status_8(winSize, GlbCnt, nbin, acc_mean, wd, LeftBnd, RightBnd, wdWrk1, acf_val, pitch1, pitch2)
% 加入凹? version
% part I: determine floating ( wd=1) status to be 0(non-walk) or 2 (walk)
cc=zeros(1,3); 
for (i=1:winSize:nbin) 
    jbgn = i; 
    jend = min(jbgn + winSize-1,nbin);
    cc(1:3) = 0; % cc(1): wd=0; cc(2):wd=1; cc(3): wd=2;
    for (j=jbgn:jend)
        cc(wd(j)+1) = cc(wd(j)+1) + 1; % idx = wd(j); cc(idx+1) = cc(idx+1) + 1;
    end % j
    
    if (cc(2) == 0)
        continue;
    end
    
    majority = 2*(cc(1) <= cc(3));
    for (j=jbgn:jend)
        if (wd(j) == 1)
            wd(j) = majority;
        end
    end % j
end % i
%*************簡化**************
%{
for (iloop=1:2) % kill 孤島
    wdWrk1(1:nbin) = wd; % store
    if (iloop==1)
        if (GlbCnt.Di == 25)
            span = 10; prob1 = 41;
        elseif (GlbCnt.Di == 50)
            span = 5; prob1 = 21;
        end
    else
        if (GlbCnt.Di == 25)
            span = 2;  prob1 = 51;
        elseif (GlbCnt.Di == 50)
            span = 1;  prob1 = 41;
        end
    end
    
    for (i=1:nbin)
%         if (LeftBnd(i) > 0 || RightBnd(i) > 0)
%             fprintf('How Come ?\n'); 
%             pause;
%         end
        %0610 rm redundent code dis = RightBnd(i) - LeftBnd(i) + 1;
        % if (dis >= trunkSize || wd(i)==0) %1st: neglect big trunk; 2nd: save run time 
        if (wd(i)==0) %1st: neglect big trunk; 2nd: save run time 
            continue;
        end
        j1= max(i-span,1); j2= min(i+span,nbin);
        c1 = 0; c2 = 0; %neglect BND now
        for (j=j1:i)
            c1 = c1 + (wd(j) > 0);
        end % j
        for (j=i:j2)
            c2 = c2 + (wd(j) > 0);
        end % j
        if (100*c1 <= prob1*span && 100*c2 <= prob1*span) % pr1 == probability
            wdWrk1(i) = 0; 
        end
   end % i
   wd = wdWrk1(1:nbin); % restore
end % iloop

% part III: % add new walk duration (in which amplitude is much smaller than normal % walking, ..)
[numBnd, LeftBnd, RightBnd] = find_walk_boundary(nbin, LeftBnd, RightBnd, wd);
 
for (i=1:numBnd) % add inertia time
    i1 = min(RightBnd(i) + 1,nbin); % SINSIN 6 is default maxium ghost length
    if (GlbCnt.Di == 25)
        i2 = min(RightBnd(i) + 8,nbin); % SINSIN 6 is default maxium ghost length
    else
        i2 = min(RightBnd(i) + 4,nbin);%for1
    end
    maxNum = i2-i1+1;
    % I am here SINSIN 0430/2016
    wdWrk1(i1:i2) = wd(i1:i2);

    [wdWrk1, pitch1, pitch2, acf_val] = get_initial_walk_status_part_II(GlbCnt, i1, i2, 2, 0, acc_mean, pitch1, pitch2, acf_val, wdWrk1); % 0428
    num = 0;
    jlast = -1;
    for (j=i1:i2) % travel ghost points
        if (wdWrk1(j) ==2) 
            num = num + 1;
            jlast = j;
        end
    end % j
    if (GlbCnt.Di == 25)
        acceptableNum = 2;
    elseif (GlbCnt.Di == 50)
        acceptableNum = 1;
    end
%     if (num >= maxNum-2)
    if (num >= maxNum - acceptableNum) %for1
        for (j=i1:jlast) % travel ghost points
            wd(j) = 2;
        end
    end
end % i
     
[numBnd, LeftBnd, RightBnd] = find_walk_boundary(nbin, LeftBnd, RightBnd, wd); %0328F8
% [wd] = kill_short_segment4(GlbCnt.Sec1,2*GlbCnt.Sec0p5,nbin,wd, LeftBnd, RightBnd);
[wd] = kill_short_segment4(GlbCnt.Sec1, GlbCnt.Sec1, nbin, wd, LeftBnd, RightBnd); % for1
[wd] = kill_short_segment4(GlbCnt.Sec4, GlbCnt.Sec3, nbin, wd, LeftBnd, RightBnd);
%}  
[numBnd, LeftBnd, RightBnd] = find_walk_boundary(nbin, LeftBnd, RightBnd, wd); %0328F8

return; %trim_walk_status_8 acc

function [wd] = kill_short_segment4(thresh1, thresh2, nbin, wd, LeftBnd, RightBnd) % new
% consider solid boundary (cannot change status !! ??); find_boundary of different wd
nr = 1;
LeftBnd(1) = 1;
for (i=2:nbin) % HAs Bug ?? questioned 0401-2016 , pls refer to find_walk_bnd 
    if (wd(i) ~= wd(i-1) )
        RightBnd(nr) = i-1;
        nr = nr + 1;
        LeftBnd(nr) = i;
    end
end
RightBnd(nr) = nbin;

% for (r=1:nr) %  DEBUG ONLY
%     for (i=LeftBnd(r):RightBnd(r)-1) assert(wd(i) == wd(i+1)); end
% end  %  DEBUG ONLY

prevStatus = wd(LeftBnd(1));
for (r=2:nr-1)
    curStatus = wd(LeftBnd(r));
    if (prevStatus ~= curStatus && (RightBnd(r+1) - LeftBnd(r+1)) >= thresh2 && (RightBnd(r-1)-LeftBnd(r-1)) >= thresh2 && (RightBnd(r) - LeftBnd(r)) <= thresh1)
        for (i=LeftBnd(r):RightBnd(r))
            wd(i) =prevStatus;
        end 
    end
    prevStatus = curStatus; % it is correct
end % r

return; % kill_short_segment4

function [numBnd, LeftBnd, RightBnd] = find_walk_boundary(nbin, LeftBnd, RightBnd, wd)
%  C code should be very simple
c1 = 0; numBnd = 0;
if (wd(1) > 0) 
    LeftBnd(1) = 1;  c1 = 1;
end
for (i=2:nbin)
    if (wd(i) > 0 && wd(i-1) ==0 )
        c1 = c1 + 1;
        LeftBnd(c1) = i;
    end
end
for (i=1:nbin-1)
    if (wd(i) > 0 && wd(i+1) ==0)
        numBnd = numBnd + 1;
        RightBnd(numBnd) = i;
    end
end
if (wd(nbin) > 0)
    numBnd = numBnd + 1;  
    RightBnd(numBnd) = nbin;
end
assert(c1==numBnd && nbin > 0);
return; % find_walk_boundary

function [pitch1, pitch2] = modify_irratic_pitch_7(GlbCnt, pitch1, pitch2)

numBin = length(pitch1);
for (i=1:numBin) % it  will get smoother p2ave for exp1230 r=2; % yhw0405DD move here
    if (pitch1(i) >= 4*GlbCnt.Fs)
        maxBin = min(numBin,i+GlbCnt.Sec4);
        for (j=i+1:maxBin)
            if (pitch1(j) < 4*GlbCnt.Fs)
                pitch1(i) = pitch1(j);
                pitch2(i) = pitch2(j);
                break;
            end
        end % j
    end % pitch = 0
end % i

return % modify_irratic_pitch_7

function  [sumStep] = segment_step_by_integration_0609_2016(GlbCnt, binBgn, binEnd, pitch1, pitch2)

fac_freq = GlbCnt.Fac_Step;
sumStep = 0;
for (i=binBgn:binEnd)
    
    freq = floor(GlbCnt.Fac_Step*GlbCnt.Fs/pitch1(i));       

    if (pitch1(i) ~= pitch2(i)) 
        freq = 2*freq;  
    end
    
    sumStep = sumStep + freq;
end % i   
sumStep = round(sumStep*GlbCnt.Di/(fac_freq*GlbCnt.Fs));

return; % segment_step_by_integration

function [Eul0] = modify_sim_and_pitch_12(GlbCnt, Eul0)
% 0416-2016
fac_acf = GlbCnt.Fac_acf;
nbin = Eul0.nbin;
% *********************** part 1, normal freq 1Hz ~ inf 3Hz ***********************
for (jj=1:nbin) 
    jj0 = max((jj) - (GlbCnt.Sec0p5),1);     
    amp0 = max(Eul0.max_acc(jj0:jj)) - min(Eul0.min_acc(jj0:jj));  
    jj1 =  min(jj +1,nbin);    
    jj2 =  min(jj+ 1 + (GlbCnt.Sec1),nbin); 
    amp1 = max(Eul0.max_acc(jj1:jj2)) - min(Eul0.min_acc(jj1:jj2));
    if (GlbCnt.FacAmp*amp0 < 300) %  possibly non-walk; try to 降低 non-walk status 的 false positive   
        thresh2 = 5;
    else
        thresh2 = 2;  
    end
    
    if (10*amp0 >= thresh2*amp1) % not a ABI
        if (GlbCnt.FacAmp*Eul0.acc_mean(jj) > 1800) % 0424-2016 move here % 高速運動=> 增加權重
            if (Eul0.pitch2(jj) < floor(GlbCnt.Fs/3))
                Eul0.acf_val(jj) = min(Eul0.acf_val(jj) + (6*fac_acf/64),10000); % add by 0.05; yhw0409_2016fix bug
            else
                Eul0.acf_val(jj) = min(Eul0.acf_val(jj) + (fac_acf/16),10000); % add by 0.05; yhw0409_2016fix bug
            end
        end
    else  % active ACF length is too short around right end
        if (GlbCnt.FacAmp*amp0 < 300) % 
            Eul0.acf_val(jj) = 0; % value  yhw: can be next peak ??
            Eul0.pitch1(jj) = 4*GlbCnt.Fs; % assign an impossible long pitch; correct later
            Eul0.pitch2(jj) = 4*GlbCnt.Fs; % need to penalizethe pitch ?
        % since corr is an average concept, take it previous value to make it stable occurs (ABI) @ 
        % transisition from  (A)non-walk -> walk or walk; del FP; (B) walk ->手亂揮, the previous one is reliable
        else
            jn1 = max(jj-1,1);
            Eul0.acf_val(jj) = Eul0.acf_val(jn1); %  
            Eul0.pitch1(jj) = Eul0.pitch1(jn1); %   
            Eul0.pitch2(jj) =  Eul0.pitch2(jn1); %   
        end
    end
end % jj %  count count count
return; % modify_sim_and_pitch_12
 
function [acf_val, pitch1, pitch2] = fast_cal_similarity_0717(dinc, GlbCnt, max_lag0, acc0)
% map to coarser grid;
%dinc: delte increment for acf calculation, if dinc > 1=> speed up; suggest dinc=1 or 2 
%choose dinc such that coarse grid sample rate >=16Hz , but coarse grid sf >=32 for running with  pitch <=0.25 sec 
if (GlbCnt.NoDownSample)
    dinc = 1;
end
assert(dinc >=1);
max_lag = floor(max_lag0/dinc);
acc = acc0(1:dinc:end);

% figure(13);clf;
% plot(acc);

[acf_val, pitch1, pitch2] = cal_similarity0602_2016(GlbCnt, max_lag, acc);
% [1];
if (acf_val < 2)
   return; 
end

if (GlbCnt.NoDownSample)
    return;
end

if (pitch1 <= GlbCnt.Fs/4 && (pitch1-pitch2) < 3) % single wave high frequency motion, go back to original fine grid coord
    [acf_val, pitch1, pitch2] = cal_similarity0602_2016(GlbCnt, max_lag0, acc0); %SINSIN Should be cal_similarity0602_2016 ??
    return; 
end
% *********************************************************
idx = pitch1 + 1;
if (pitch2 ~= pitch1)
    idx_slave = pitch2 + 1;
else
    idx_slave = -1;
end
% map to fine grid (riginal coordinate) check nbrs
ave0 = mean(acc0); %0606-2016
idx_nbr(2) = dinc * idx - 1 ; % map to fine grid
idx_nbr(1) = idx_nbr(2) - 1;
idx_nbr(3) = idx_nbr(2) + 1;
    
x1 = acc0(1:end - idx_nbr(1) + 1);  x2 = acc0(idx_nbr(1):end);
acf_nbr(1) = sum((x1-ave0).*(x2-ave0));
    
x1 = acc0(1:end - idx_nbr(2) + 1);  x2 = acc0(idx_nbr(2):end);
acf_nbr(2) = sum((x1-ave0).*(x2-ave0)); % mid
    
x1 = acc0(1:end - idx_nbr(3) + 1);  x2 = acc0(idx_nbr(3):end);
acf_nbr(3) = sum((x1-ave0).*(x2-ave0));
    
acf_max = acf_nbr(1); idx_max = 1;
for (i=2:3) % 1: left nbr; 2: center; 3: right nbr
    if (acf_nbr(i) > acf_max)
        idx_max = i; acf_max = acf_nbr(i);
    end
end
assert(idx_max > 0);
pitch1 = idx_nbr(2) + idx_max - 2;
if (idx_slave == -1)
    pitch2 = pitch1;
else
    pitch2 = dinc * idx_slave - 1; % use coarse grid; mau lose accuracy (resolution)
end
pitch1 = pitch1 - 1;  %SINSIN FIX BUG YHW0507
pitch2 = pitch2 - 1;

return; % fast_cal_similarity_0717  acf_val jdx

function [acf_val, pitch1, pitch2] = cal_similarity0602_2016(GlbCnt, max_lag, acc)

nacf = max_lag + 1;
% acf = autocorr(acc,max_lag);
% acf = autocorr_matlab_edit_brute_force(acc, max_lag);
acf = autocorr_brute_force_CY(acc, max_lag, length(acc), 1);
acf00 = acf(1);

% figure(14);clf;
% plot(acf);

% ave = mean(acc); %0606-2016
% acf00 = sum((acc-ave).*(acc-ave));
% acf = acf00*acf;
    
[acf_val_max, idx] = step_local_max(0, 1, nacf, acf);

if (idx > 0) % need to modify! say, 32Hz
    if (10*idx < 35*GlbCnt.Fs && acf(max(idx-3,1)) > acf(idx)) % SINSIN flase local max; FIX BUG 0423-2016
        [acf_val_max, idx] = step_local_max(0, min(idx+1,max_lag), nacf, acf);
    end
end

if (idx <= 1 || length(acf) < 10)
    acf_val=0; pitch1= 4*GlbCnt.Fs;  pitch2=pitch1; 
    return;
end
    
pitch1 = idx;  pitch2 = idx;
[acf_val] = normalized_acf_0602_2016(GlbCnt, acc, acf_val_max, acf00, idx);

jdx1 = max(GlbCnt.Num_wave_in_bin+2,floor(idx*3/8)); % has double frequency ?
jdx2 = floor(idx*5/8); jdx2 = min(jdx2+1,max_lag);
[m_slave, idx_slave] = step_local_max(0, jdx1, jdx2, acf);

if (idx_slave > jdx1) % if idx_slv=1=> must be bnd point and no local max
    pitch2 = idx_slave; % double frequency found
end
  
pitch1 = pitch1 - 1;  %SINSIN FIX BUG YHW0507
pitch2 = pitch2 - 1;

return; % cal_similarity0602_2016  

function  [acf_val] = normalized_acf_0602_2016(GlbCnt, x, acf_val_max, acf00, idx1)
% we use the pitch corresponding to the max value of the acf, then recaluclate the similarity again (when val1 > 0.2; 
% otherwise it seems no hope to get high similarity coefficient)
fac_acf = GlbCnt.Fac_acf;
acf_val_max = floor(acf_val_max);

if (10*acf_val_max < (2*acf00)) 
    acf_val = (fac_acf*acf_val_max/acf00); % normalized
    return;
end

n = length(x);
A = x(1:n-idx1+1);
B = x(idx1:n);
nA = norm(A);
nB = norm(B);

%  p1 = sum(A.*B)/(nA*nB); r = min(nA,nB)/max(nA,nB); then def acf_val  =
%  p1*sqrt(ratio);
v0 = floor(fac_acf*sum(A.*B)/(nA*nB)); % more correct value than acf  
if (nA > nB)
    max2 = nA; min2 = nB;
else
    max2 = nB; min2 = nA;
end

if (max2 >= 3*min2) % YHW0330 MODIFY to be x/256; r(==max2/min2) > 3
    acf_val = floor(147*v0)/256;
elseif (max2 >= 2*min2) % 2 <= r < 3, approximate the sqrt by a line in each region
    acf_val = floor((-floor(33*max2/min2) + 247)*v0/256);
elseif (max2 >= min2) % 1 <= r < 2
    acf_val = floor((-floor(75*max2)/min2 + 331)*v0/256);
else
    acf_val = floor(v0);
end  
%*******簡化*********
acf_val = floor(v0);
%********************
return; % normalized_acf_0602_2016  acf_val0

function [init_max, idx] = step_local_max(init_max, i1, i2, x)
%SINSIN yhw0418-2016 BUG fixed 
% return with local max > input is not a min
idx = -1;
for (i=i1+1:i2-1)
    if (x(i) >= x(i+1) && x(i) >= x(i-1) && (2*x(i) > x(i+1) + x(i-1)) ) % one equal sign local max
        if (x(i) > init_max)
            init_max = x(i);
            idx = i;
        end
    end
end % i
return; % step_local_max

function [del_step, Sstamp] = dump_small_step_0713(GlbCnt, Sstamp) 
% BBB
del_step = 0;  
if (Sstamp(2).bgn - Sstamp(1).end > GlbCnt.Fs)  % 分開的兩個區間: 處理上一個連續區間
    assert(Sstamp(2).bgn >= Sstamp(1).end);
    is_small_step = (Sstamp(1).end-Sstamp(1).bgn) < GlbCnt.Fs*GlbCnt.MinWalkTime ||  Sstamp(1).step < GlbCnt.MinWalkStep;
    if (is_small_step) 
       % fprintf('AB %d %d\n',condA, condB);
      del_step = Sstamp(1).step;
    end
    Sstamp(1).bgn = Sstamp(2).bgn;
    Sstamp(1).end = Sstamp(2).end;
    Sstamp(1).step = Sstamp(2).step;
else % 因為online boundary)造成連接的 n 個區間
    Sstamp(1).end = Sstamp(2).end;
    Sstamp(1).step = Sstamp(1).step + Sstamp(2).step;
end
return; % dump_small_step_7 Cont_step

function [GlbCnt] = specify_global_constants_2022_1013(fs0)

GlbCnt.Fs = fs0;
GlbCnt.Dinc = 2;
GlbCnt.FacAmp = 1; %  factor of acceleration resolution reduction
% ********************************************************************** CY
GlbCnt.BinSizeCase = 1; % (1) 0.5 sec (2) 1 sec
GlbCnt.FilterType = 1; % (0) no filter(minus the average) (1) IIR (2) UPEMD
GlbCnt.IIR_passband = 3; % (1) low pass (2) high pass (3) bandpass
GlbCnt.NoDownSample = 1; % (1) turn off down sample autocorr 
% ****************************************************************** CY_end

if (GlbCnt.BinSizeCase == 1)
    GlbCnt.Di = 25;
elseif (GlbCnt.BinSizeCase == 2)
    GlbCnt.Di = 50;
end

if (GlbCnt.FilterType == 1)
    if (GlbCnt.IIR_passband == 1)
        [GlbCnt.CB,GlbCnt.CA] = butter(2, [24/50*2], 'low');
    elseif (GlbCnt.IIR_passband == 2)
        [GlbCnt.CB,GlbCnt.CA] = butter(2, [0.75/50*2], 'high');
    elseif (GlbCnt.IIR_passband == 3)
%         [GlbCnt.CB,GlbCnt.CA] = butter(2,[0.75/GlbCnt.Fs*2 24/GlbCnt.Fs*2], 'bandpass'); % acc base
        [GlbCnt.CB,GlbCnt.CA] = butter(2,[0.75/GlbCnt.Fs*2 24/GlbCnt.Fs*2], 'bandpass');
    end
end

GlbCnt.Win_acf = floor(3*GlbCnt.Fs);
GlbCnt.Win_acfDi = floor(GlbCnt.Win_acf/GlbCnt.Di); %  AAA

GlbCnt.Max_lag = floor(3*GlbCnt.Fs/2);  

% step_counting_stage    
GlbCnt.Amin = floor(200/GlbCnt.FacAmp);%   慣用腳的最小加速度
GlbCnt.Phigh = floor(141*GlbCnt.Fs/64); % acf based only  
% 本程式可容許之最高freq (or pitch) i.e. 人的步行極限; 若pitch 設低一點, 可減少誤動作,
%但若設的太高,則不能 catch 項下樓梯的step. generally speaking down stair <= 6Hz; splint <= 5Hz; derived from acf 
GlbCnt.Plow = ceil(16*GlbCnt.Fs/100); 
GlbCnt.MinWalkTime = 4; % in sec 2016-0614 chg 
GlbCnt.MinWalkStep = 8;

GlbCnt.Fac_acf = 10000; % default 1 or 1000; yhw0323-2016 chg to 10^4

GlbCnt.Fac_Step = 1000; % Euler_Step  0609

GlbCnt.Ovlp = floor(5*GlbCnt.Fs/2); % SIN SIN should use GlbCnt.win_acf - GlbCnt.Di later;
GlbCnt.Phys_win_size = floor(20*GlbCnt.Fs); %0204-2016

if (GlbCnt.BinSizeCase == 1)
    GlbCnt.OvlpDi = floor(GlbCnt.Ovlp/GlbCnt.Di);
    GlbCnt.Ovlp = GlbCnt.OvlpDi * GlbCnt.Di;
    GlbCnt.OvlpBin = 2 * GlbCnt.OvlpDi; % OvlpBin = 10 (5 sec)
elseif (GlbCnt.BinSizeCase == 2)
    GlbCnt.OvlpDi = floor(GlbCnt.Ovlp/GlbCnt.Di) + 1; % 0729 BBB
    GlbCnt.Ovlp = (GlbCnt.OvlpDi - 1) * GlbCnt.Di; %BBB
    GlbCnt.OvlpBin = 2*GlbCnt.OvlpDi - 1; % OvlpBin = 5 (5 sec)
end

GlbCnt.Num_wave_in_bin = floor(GlbCnt.Fs/8); % in sec 每個 bin 最多有一個 wave  !0.13sec  
%0630-2016 comment out 

Fsdi = (GlbCnt.Fs/GlbCnt.Di);
GlbCnt.Sec1 = floor(Fsdi);
GlbCnt.Sec0p5 =  floor(GlbCnt.Sec1/2); 
GlbCnt.Sec2 = 2*GlbCnt.Sec1;
GlbCnt.Sec3 = 3*GlbCnt.Sec1;
GlbCnt.Sec4 = 4*GlbCnt.Sec1;
GlbCnt.Sec5 = 5*GlbCnt.Sec1;

assert(floor(GlbCnt.Sec0p5)==ceil(GlbCnt.Sec0p5));
assert( GlbCnt.Fs >= 24 &&  GlbCnt.Fs <= 100);
assert(GlbCnt.Max_lag <= GlbCnt.Win_acf);
assert(GlbCnt.Phys_win_size > GlbCnt.MinWalkTime); % DEBUG ONLY
% DDD assert(ceil(GlbCnt.Fs/GlbCnt.Di) == floor(GlbCnt.Fs/GlbCnt.Di)); %0602-2015 added

% ************************************************************ 2022_0805_DY
GlbCnt.Max_acc_length = GlbCnt.Win_acf;
% ******************************************************** 2022_0805_DY_end

if (GlbCnt.BinSizeCase == 1)
    GlbCnt.MaxCompBin1 = GlbCnt.Phys_win_size/GlbCnt.Di; % Euer Bin 的長度 for cwin=1    
    GlbCnt.MaxCompBin2 = ((GlbCnt.Phys_win_size + 2*GlbCnt.Ovlp )/GlbCnt.Di); % for cwin > 1
elseif (GlbCnt.BinSizeCase == 2)
    GlbCnt.MaxCompBin1 = GlbCnt.Phys_win_size/GlbCnt.Di;
    GlbCnt.MaxCompBin2 = ((GlbCnt.Phys_win_size + 2*GlbCnt.Ovlp )/GlbCnt.Di) + 1;
end

GlbCnt.Ncomp = max(GlbCnt.Max_acc_length,GlbCnt.Ovlp); % 

%SINSIN only true for FacAmp = 1
GlbCnt.GminTable = [1440 1440 1440 1440 1440 1440 1440 1440 1440 1440 1440 1440 1366 1234 1120 1020 895 726 590 480 389 314 253 202 160 142 126 112 99 88 79 70 62 55 48 43 37 33 28 25 21 18 15 12 9 7 5 3 1 0 0 0 0 0 0 0 0 0 0 0]; 

pmax = ceil(1.2*GlbCnt.Fs) + 2; % 2 : safety factor
for (pp = 1:pmax)
    % code  pp = pitch2(jj); pp = max(pp,1); pp = min(pp,ceil(1.2*GlbCnt.Fs));
    gmin = min_acc_freq_table1(GlbCnt,pp,GlbCnt.Amin) - GlbCnt.Amin;
    GlbCnt.GminTable(pp) = gmin;
end 

if (1) % % **************    Generate C Code Script     ******************* 
    fprintf('C code Defintion Scripy Begins\n\n\n\n');
    fprintf('#define Fs %d\n',GlbCnt.Fs);
     
    fprintf('#define Dinc %d\n',GlbCnt.Dinc); 

    fprintf('#define Fac_acf %d\n',GlbCnt.Fac_acf);
    fprintf('#define FacAmp %d\n',GlbCnt.FacAmp); 
    fprintf('#define Ovlp %d\n',GlbCnt.Ovlp);
    fprintf('#define  Phys_win_size %d\n',GlbCnt.Phys_win_size);

    fprintf('#define Di %d\n',GlbCnt.Di);
    fprintf('#define Win_acf %d\n',GlbCnt.Win_acf);

    fprintf('#define OvlpDi %d\n',GlbCnt.OvlpDi);
    fprintf('#define Win_acfDi %d\n',GlbCnt.Win_acfDi);

    fprintf('#define Max_acc_length %d\n',GlbCnt.Max_acc_length); 
    fprintf('#define MaxCompBin1 %d\n',GlbCnt.MaxCompBin1); 
    fprintf('#define MaxCompBin2 %d\n',GlbCnt.MaxCompBin2); 
    fprintf('#define Ncomp %d\n',GlbCnt.Ncomp); 
   
    fprintf('#define Max_lag %d\n',GlbCnt.Max_lag);
    fprintf('#define Phigh %d\n',GlbCnt.Phigh);
    fprintf('#define Plow %d\n',GlbCnt.Plow);
    fprintf('#define Num_wave_in_bin %d\n',GlbCnt.Num_wave_in_bin);
   
    fprintf('#define Amin %d\n',GlbCnt.Amin);
 
    fprintf('#define MinWalkTime %d\n',GlbCnt.MinWalkTime);
    fprintf('#define MinWalkStep %d\n',GlbCnt.MinWalkStep);
    fprintf('#define Sec1 %d\n',GlbCnt.Sec1);
    fprintf('#define Sec0p5 %d\n',GlbCnt.Sec0p5);
    fprintf('#define Sec3 %d\n',GlbCnt.Sec3);
    fprintf('#define Sec4 %d\n',GlbCnt.Sec4); 
    fprintf('#define Sec5 %d\n',GlbCnt.Sec5);
  
    assert(length(GlbCnt.GminTable) < 2^15);  
    ntab = length(GlbCnt.GminTable);
    fprintf('short int GminTable[%d] = {0,',ntab+1);
    
    for (k=1:length(GlbCnt.GminTable))
        if (k < length(GlbCnt.GminTable))
            fprintf('%d, ',GlbCnt.GminTable(k));
        else
            fprintf('%d};\n',GlbCnt.GminTable(k));
        end
    end
    
  % ********************** Print filter coefficient **********************
    if (GlbCnt.FilterType == 1)
        fac_filter = 32768;  
        fprintf('int CA[6] = {0, '); 
        for (k=1:length(GlbCnt.CA))
            if (k < length(GlbCnt.CA))
                fprintf('%d, ',round(fac_filter*GlbCnt.CA(k)));
            else
                fprintf('%d};',round(fac_filter*GlbCnt.CA(k)));
            end
        end % k
        fprintf(' //coefficient have been normalized by fac_filter true for  fac_filter=%d, fs=%d, passband = 0.75~24Hz\n',fac_filter,GlbCnt.Fs);
        fprintf('int CB[6] = {0, ');
        for (k=1:length(GlbCnt.CB))
            if (k < length(GlbCnt.CB))
                fprintf('%d, ',round(fac_filter*GlbCnt.CB(k)));
            else
                fprintf('%d};\n',round(fac_filter*GlbCnt.CB(k)));
            end
        end % k
    end
    fprintf('\n\n\n');
    fprintf(' C code Defintion ends\n');
end % **************    End Generate C Code Script    **************
return;  % specify_global_constants_2022_1013

function gmin = min_acc_freq_table1(GlbCnt,p,gLowLim)
% (freq	acc value) = (1.0  600), (1.5 1000), (2.0 1400), (2.5  1800), (3.0 2200)
% this table is used to generate gminTable
%SINSIN Integer version ?? yhw0410-2016 question
f = GlbCnt.Fs/p;
f = min(4,f);
if (f < 1)
    c2 = 0;
elseif (f < 2) % f=[1,2]=> c=[0.1 0.2]
    c2 = f;
elseif (f < 3) % f=[2,3] => c=[0.2 0.6]
    c2 = 4*f-6;
else
    c2 = 6;
end  
gmin = floor(gLowLim + c2*(80)*(f-1));
return; % min_acc_freq_table
% *************************************************************************
function [ax, ay, az, fileName] = input_accelerometer_data4_new_format(GlbCnt, icase, k)

ax = []; ay=ax; az=ax; fileName=[];
dir1 = '\experiment2\';   

if (icase == 1) % compare ours to fitbit and ground truth
    file_index_list = [1 2 3 4 5 6 7  101 102]; 
elseif (icase == 2) %  compare ours to  ground truth
    file_index_list =  [1 2 5 6 7 8 9  101 102]; 
elseif (icase == 3) % run single case
    file_index_list =  [1 2 3 4 5 6 7 8 9 10 11 12 101 102]; 
end

fileName{1} = 'exp-walk-1230-2014';
fileName{2} = 'exp-wrist-0115-2015';
fileName{3} = 'exp_positive_0312';
fileName{4} = 'exp_positive_0318';
fileName{5} = 'exp_0320';
fileName{6} = 'exp_0325';
fileName{7} = 'exp_0331';
fileName{8} = 'exp1201-2014';
fileName{9} = 'paper-data';
%fileName{10} = 'exp-waist-0115-2015'; wrong data: a lot of zero readings
fileName{11} = 'total-ieee';
fileName{101} = 'exp-non-walk-1230-2014';
fileName{102} = 'exp_negative_0313';

fileName{200} = []; 

if (isempty(intersect(file_index_list,k)) == 1)
    return; 
end

filenameIn = strcat(fileName{k},'.mat');
pp = [cd dir1 filenameIn];

if (exist(pp) == 0)
   return; 
end
load(pp);

ax = ((ax/GlbCnt.FacAmp)); %  0611-2016
ay = ((ay/GlbCnt.FacAmp));
az = ((az/GlbCnt.FacAmp));
return; % input_accelerometer_data4_new_format