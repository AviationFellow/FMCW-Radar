% Reference: MathWorks example "Design and Simulate an FMCW Long-Range Radar (LRR)"

c = 3e8; % Speed of light (m/s)
fc = 77e9; % Operating frequency (Hz)
lambda = c / fc;% Wavelength (m)

% Target parameters
targetRange = 300;% True target range (m)
targetVelocity = -35;% m/s, (approaching radar)
targetAngle = 0;% True target angle (measured from boresight)

% Radar design specifications
rangeMax = 500; % Maximum detection range (m)
rangeRes = 1; % Desired range resolution (m)
vMax = 50; % Maximum velocvity of target(m/s)
tRoundTrip = 2 * rangeMax / c;
tm = 5.5 * tRoundTrip;  % Sweep time (s)

% FMCW parameters:
bw = c / (2 * rangeRes);  % Sweep bandwidth (Hz)
sweepSlope = bw / tm; % Sweep slope (Hz/s)

% Estimate maximum beat frequency and Doppler shift:
fbeatMax = sweepSlope*(2*rangeMax/c);
fdopMax = 2 * abs(targetVelocity)/lambda; % round-trip Doppler (Hz)
fIFMax = fbeatMax + fdopMax;
fs = max(2 * fIFMax, bw);  % Sample rate (Hz)

% Create FMCW waveform object
waveform = phased.FMCWWaveform('SweepTime', tm,'SweepBandwidth', bw,'SampleRate', fs);
Ns = waveform.SampleRate * waveform.SweepTime; % Samples per sweep

%Tx, Rx Setup
tx = phased.Transmitter('PeakPower', 1e3, 'Gain', 30);
rx = phased.ReceiverPreamp('Gain', 30, 'NoiseFigure', 4.5, 'SampleRate', fs);

rxArray = phased.ULA('NumElements', 2, 'ElementSpacing', lambda/2); %two-element uniform linear array for angle estimation

% Determine initial target position based on range and angle.
% Here, x is the range direction and y is the cross-range.
targetPos = [targetRange * cosd(targetAngle); targetRange * sind(targetAngle); 0];
% Assume the target is approaching: velocity components accordingly.
targetVel = [-targetVelocity * cosd(targetAngle); targetVelocity * sind(targetAngle);0];

% Create a platform object for the target.
targetPlatform = phased.Platform('InitialPosition', targetPos, 'Velocity', targetVel);

% Used free-space propagation model (Friss)
channel = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc,'SampleRate', fs, 'TwoWayPropagation', true);

% Simulation Loop: Generate Multiple Sweeps for Doppler Processing
Nsweeps = 64;  % sweeps for Doppler estimation
rxData = complex(zeros(Ns, Nsweeps, 2)); % Store dechirped signal for two antennas
txSignal = waveform();  % Reference transmitted FMCW chirp

for m = 1:Nsweeps
    % Update target position (simulate motion)
    [tgtPos, tgtVel] = targetPlatform(waveform.SweepTime);
    radarPos = [0; 0; 0];        % Radar assumed at origin
    radarVel = [0; 0; 0];        % Radar stationary
    
    % Generate and transmit FMCW signal
    sig = txSignal;
    sig = tx(sig);
    
    % Propagate signal to target and back (reflection assumed with RCS=1)
    sig = channel(sig, radarPos, tgtPos, radarVel, tgtVel);
    
    [az, el, ~] = cart2sph(tgtPos(1), tgtPos(2), tgtPos(3));
    angles = [rad2deg(az); rad2deg(el)]; 

    % Use 'angles' as the direction of arrival for collectPlaneWave
    sig = collectPlaneWave(rxArray, sig, angles, fc);

    sig = rx(sig);
    
    % Dechirp: mix the received signal with the transmitted FMCW signal
    dechirped = dechirp(sig, txSignal);
    
    % Store dechirped data for each antenna element
    rxData(:, m, :) = dechirped;
end

% Signal Processing for Range and Doppler Estimation
% Process first antenna data for range estimation
data1 = squeeze(rxData(:, :, 1));  % Data from antenna 1
rangeFFT = fftshift(fft(data1, [], 1), 1);
f_range = (-Ns/2 : Ns/2-1)' * (fs / Ns);  % Frequency vector for range FFT
% Convert beat frequency to range using: R = (f_beat * c) / (2 * sweepSlope)
rangeAxis = (f_range * c) / (2 * sweepSlope);

% Average the magnitude over sweeps to obtain a range profile.
magRange = mean(abs(rangeFFT), 2);

% For Doppler (velocity) estimation, select the range bin with the maximum average magnitude.
[~, maxIdx] = max(magRange);
dopplerData = fftshift(fft(rangeFFT(maxIdx, :), [], 2), 2);
fd = (-Nsweeps/2 : Nsweeps/2-1) * (1 / (Nsweeps * tm));  % Doppler frequency vector
velocityAxis = fd * lambda / 2;  % v = (fd * lambda) / 2
[~, maxDopIdx] = max(abs(dopplerData));
estimatedVelocity = velocityAxis(maxDopIdx);

% SNR Estimation
% Compute received signal power and estimate noise power.
tau = tm;
receivedPower = mean(abs(data1(:)).^2);
k = 1.38e-23; T = 290;           % Boltzmann constant and noise temperature
noisePower = k * T * (1 / tau) * 10^(4.5/10);  % Approximate noise power (using bandwidth ~1/tau)
SNR_sim = 10 * log10(receivedPower / noisePower);

% Angle Estimation (Using Two-Element Array)
% Compute phase difference between the two antennas.
data2 = squeeze(rxData(:, :, 2));  % Data from antenna 2
phaseDiff = angle(data1) - angle(data2);
avgPhaseDiff = mean(phaseDiff(:));
d = rxArray.ElementSpacing;
% Estimate angle using: theta = asin((phaseDiff * lambda)/(2*pi*d))
estimatedAngle = asind((avgPhaseDiff * lambda) / (2*pi*d));

% Plotting Results
% Create a figure with three subplots for range profile, Doppler profile, and frequency spectra.
figure;
subplot(3,1,1);
plot(rangeAxis, 20*log10(abs(magRange)), 'LineWidth', 1.5);
xlabel('Range (m)'); ylabel('Magnitude (dB)');
title('Range Profile'); grid on;
subplot(3,1,2);
plot(velocityAxis, abs(fftshift(fft(rangeFFT(maxIdx, :)))), 'LineWidth', 1.5);
xlabel('Velocity (m/s)'); ylabel('Magnitude');
title('Doppler (Velocity) Profile'); grid on;
% Compute FFT of the transmitted and one received dechirped signal (first sweep, first antenna)
TxSpectrum = fftshift(fft(txSignal));
RxSpectrum = fftshift(fft(squeeze(rxData(:,1,1))));
f_axis = f_range / 1e6;  % Convert frequency axis to MHz
subplot(3,1,3);
plot(f_axis, 20*log10(abs(TxSpectrum)), 'b', 'LineWidth', 1.5);
hold on;
plot(f_axis, 20*log10(abs(RxSpectrum)), 'r', 'LineWidth', 1.5);
xlabel('Frequency (MHz)'); ylabel('Magnitude (dB)');
title('Frequency Spectrum: Transmitted (blue) vs. Received (red)'); grid on;
legend('Tx Spectrum', 'Rx Spectrum');

% Display Estimated Parameters
[maxV, indx] = max(20*log10(abs(magRange)));
fprintf('Estimated Range: %.2f m (True Range: %.2f m)\n', rangeAxis(indx), targetRange);
fprintf('Estimated Velocity: %.2f m/s (True Velocity: %.2f m/s)\n', estimatedVelocity, targetVelocity);
fprintf('Estimated Angle: %.2f deg (True Angle: %.2f deg)\n', estimatedAngle, targetAngle);
fprintf('Estimated SNR: %.2f dB\n', SNR_sim);

[instFreq, tInst] = instfreq(txSignal, fs);
figure;
plot(tInst, instFreq/1e6, 'LineWidth', 1.5); % Convert frequency to MHz
grid on;
xlabel('Time (s)');
ylabel('Instantaneous Frequency (MHz)');
title('Instantaneous Frequency vs Time of Transmitted Signal');
