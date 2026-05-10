% ========================================================================
% Kraken Time Front Code
% Profile: North Pacific | SD = 1300 m | Fc = 200 Hz
% ========================================================================

clear
close all
fclose all;
dbstop if error

% ------------------------------------------------------------------------
% INPUTS
% ------------------------------------------------------------------------
filename = 'NorthPacific';
savematname = 'NP_test.mat';
OutDir = 'C:\Users\abiga\Desktop\at\North Pacific\Kraken';
SSPfile = 'C:\Users\abiga\Desktop\at\North Pacific\Kraken\NorthPacific.mat';

SD = 1300;          % source depth (m)
Range = 100;        % range (km)
nmodes = 200;       % number of modes
NR = 10 * Range;    % number of ranges
NRD = 5000;         % number of receiver depths
fc = 200;           % center frequency (Hz)
Q = 4;              % Q = fc / bandwidth
bw = fc / Q;        % bandwidth (Hz)

% ------------------------------------------------------------------------
% Load SSP and setup
% ------------------------------------------------------------------------
cd(OutDir);
load(SSPfile);
SSP.C = double(SSP.c);
c0 = mean(SSP.C(:));

figure
plot(SSP.C, SSP.z);
xlabel('Sound Speed (m/s)');
ylabel('Depth (m)');
title('Sound Speed Profile');
axis ij;

% ------------------------------------------------------------------------
% Estimate T from modal dispersion (min low f, max high f)
% ------------------------------------------------------------------------
f_low  = fc - bw/2;
f_high = fc + bw/2;

fprintf('\nEstimating group-speed limits between %.1f and %.1f Hz...\n', f_low, f_high);

% --- LOW frequency run ---
CreateENVfileKRAKEN(SSP, f_low, SD, Range, OutDir, filename);
CreateFLPfileKRAKEN(f_low,SD,Range, NR, NRD, 5000, OutDir, filename)
kraken(filename);
g_low_all = read_kraken_mod_binary([filename '.mod'], f_low, c0);

% --- HIGH frequency run ---
CreateENVfileKRAKEN(SSP, f_high, SD, Range, OutDir, filename);
CreateFLPfileKRAKEN(f_high,SD,Range, NR, NRD, 5000, OutDir, filename)
kraken(filename);
g_high_all = read_kraken_mod_binary([filename '.mod'], f_high, c0);

g_low_min  = min(g_low_all);
g_high_max = max(g_high_all);
% Filter out unphysical speeds
g_low_all  = g_low_all(g_low_all > 1300 & g_low_all < 1700);
g_high_all = g_high_all(g_high_all > 1300 & g_high_all < 1700);

if isempty(g_low_all) || isempty(g_high_all)
    warning('No valid group speeds found, reverting to c0.');
    g_low_min  = c0;
    g_high_max = c0;
else
    g_low_min  = min(g_low_all);
    g_high_max = max(g_high_all);
end



tdelay_slowest = (Range * 1000) / g_low_min;
tdelay_fastest = (Range * 1000) / g_high_max;
dt_spread = abs(tdelay_slowest - tdelay_fastest);
T = dt_spread + 1;  % 1 s buffer

fprintf('Min group speed (low f): %.2f m/s\n', g_low_min);
fprintf('Max group speed (high f): %.2f m/s\n', g_high_max);
fprintf('Arrival spread ≈ %.2f s → Setting T = %.2f s\n\n', dt_spread, T);

% ------------------------------------------------------------------------
% FFT setup
% ------------------------------------------------------------------------
fs = 4 * fc;
dt = 1 / fs;
N = round(fs * T);
df = fs / N;

frq = [df:df:bw];
frq = [-fliplr(frq) 0 frq] + fc;
nf = length(frq);
nyqst = ceil((nf + 1) / 2);

wind = sinc((frq - fc) / bw);
psif = complex(zeros(NRD, nf));

% ------------------------------------------------------------------------
% Run Kraken over all frequencies
% ------------------------------------------------------------------------
for i = 1:nf
    freq = frq(i);
    CreateENVfileKRAKEN(SSP, freq, SD, Range, OutDir, filename);
    CreateFLPfileKRAKEN(freq,SD,Range, NR, NRD,5000, OutDir, filename)
    kraken(filename);
    
    [~, ~, ~, ~, Pos, pressure] = read_shd([filename '.shd.mat']);
    omega = 2 * pi * freq;
    k0 = omega / c0;
    scale = exp(1i * k0) / (4 * pi);
    pfield = squeeze(pressure);
    psif(:, i) = scale * pfield(:, end);

    fclose all;
    fprintf('Freq = %.2f Hz done\n', freq);
end

figure;
imagesc(abs(pfield));
xlabel('Range index'); ylabel('Depth index');
title('Kraken Pressure Magnitude (Raw)');

% ------------------------------------------------------------------------
% Time-domain synthesis
% ------------------------------------------------------------------------
tdelay = (Range * 1000) / c0;
zmin = min(Pos.r.z); 
zmax = max(Pos.r.z);
taxis = tdelay + (0:N-1)/fs+10;
tmin = min(taxis); 
tmax = max(taxis);

ptz = zeros(NRD, N);
for iz = 1:NRD
    data = wind .* conj(psif(iz,:)) .* exp(1i * 2 * pi * frq * tdelay);
    data = [data(nyqst:nf), zeros(1, N - nf), data(1:nyqst - 1)];
    ptz(iz,:) = ifft(data);
end

zg = Pos.r.z;
for d = [150, 500, 1000, 3000]
    [~, iz] = min(abs(zg - d));
    plot(taxis, real(ptz(iz,:))); hold on;
end
xlabel('Time (s)');
ylabel('Pressure');
legend('150 m', '500 m', '1000 m', '3000 m');
title('Time-domain arrivals at multiple depths');

% ------------------------------------------------------------------------
% Plot output
% ------------------------------------------------------------------------
threshold = -70;
figure(2);
cti = 6; nsbi = 6;
ncol = -nsbi * threshold / cti;
ntv = -threshold / ncol;
scale = [threshold+cti:(0-threshold)/(ncol):0] - ntv/2;
cmap = clrscl('wcbmry', ncol);
colormap(cmap);

pkdm = max(abs(ptz(:)));
data = 20 * log10(abs(ptz) / pkdm);
data(data == -inf) = threshold;
clf
imagesc(scale, [0 1]', [scale' scale']', [threshold 0.0])
set(gca, 'position', [0.1 0.1 0.2 0.02]);
hold on
contour(scale, [0 1]', [scale' scale']', [threshold:cti:0], 'k');
hold off
set(gca, 'xtick', [threshold+cti 0]);
set(gca, 'ytick', []);
axis([threshold+cti 0 0 1]);
xlabel('Power (dB)');

taxis = taxis(end:-1:1);
axes;
set(gca, 'position', [0.1 0.2 0.8 0.65]);
imagesc(taxis, zg, flipud(data), [threshold 0.0]);
xlabel('Travel time (s)');
ylabel('Depth (m)');
set(gca, 'YDir', 'reverse');
grid on; zoom on;
set(gca, 'xlim', [tmin tmax]);
set(gca, 'ylim', [zmin zmax]);
drawnow;

save(fullfile(OutDir, savematname), 'data', 'taxis', 'zg');

% ========================================================================
% Local function: read_kraken_mod_binary
% ========================================================================
function g_all = read_kraken_mod_binary(modfile, freq, c0)
    % Robust binary .mod reader for Kraken (Fortran unformatted file)
    fid = fopen(modfile, 'rb');
    if fid < 0
        warning('Cannot open %s, using c0.', modfile);
        g_all = c0;
        return;
    end
    A = fread(fid, 'float32');
    fclose(fid);

    % Filter plausible values for real wavenumber (1/m)
    % In 100 km range at ~1500 m/s, k ≈ 2πf/c ≈ 0.84 1/m
    kvals = A(A > 0.5 & A < 1.0);       % keep only physical range
    kvals = kvals(1:min(400, numel(kvals)));   % keep first ~400 modes

    if isempty(kvals)
        % Try again as 64-bit
        fid = fopen(modfile, 'rb');
        A = fread(fid, 'double');
        fclose(fid);
        kvals = A(A > 0.5 & A < 1.0);
        kvals = kvals(1:min(400, numel(kvals)));
    end

    if isempty(kvals)
        warning('No valid k found in %s, using c0.', modfile);
        g_all = c0;
        return;
    end

    g_all = 2 * pi * freq ./ kvals;  % modal group speeds (m/s)
end
