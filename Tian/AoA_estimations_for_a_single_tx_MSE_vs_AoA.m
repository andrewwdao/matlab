clear all; 
clc;
debug_flag = 0;
plot_spectrum_flag = 0;

if debug_flag ~= 1
    plot_spectrum_flag = 0;
    close all;
end
plot_colours = {'r', 'b', 'k', 'm'};    % colours for the 4 estimators

%% Specify angles and the gain
num_txs = 1;
distance = 1;  % not used since we specify the loss manually
E_amplitude_gain = 1;
SNR_dB = -20:20:20;
true_AoA = -85:5:85;

if plot_spectrum_flag == 1   % if plotting, test for a few cases
    SNR_dB = 20;  
    true_AoA = [-89 -60 -30 0 30 60 89];
end

%% Tx Signal Model and Rx Antenna Array
% Tx Signal model
c = physconst('LightSpeed');
fc = 2.4e9; % Operating frequency (Hz), assume all tx's operating at the same centre freq
Delta_f = 1000;
sub_carrier = (1:num_txs)' * Delta_f;
Fs = 2 * max(sub_carrier);  % sample frequency
lambda = c / fc; % Wavelength
T = 1/Fs;% 8/Delta_f;
t = 0:1/Fs:(T-1/Fs);  % Time vector for the signal
P_t = ones(num_txs,1);    % Transmit signal power

element_num = 4;  % Number of elements in the ULA
ula_spacing = 0.5 * lambda;  % Element spacing (ULA)

% --- Generate transmitted signal
% Transmitted signal (simple sinusoidal signal)
s_t = sqrt(P_t) .* exp(1j * 2 * pi * sub_carrier * t);  % Complex sinusoid
E = E_amplitude_gain^2 * P_t(1) * T * Fs;% average received signal energy, E[amplitude_gain^2]*sum(s_t)

% steering vectors, each column is a steering vector of a particular angle of interest
angles = -90:1:90;   % anslges of interest
steering_vectors = zeros(element_num, length(angles));
for i = 1:length(angles)
    steering_vectors(:, i) = exp(-1j * 2 * pi * ula_spacing * (0:element_num-1)' * sind(angles(i)) / lambda);
end

%% Simulation, four estimators
num_itr = 1; % number of iterations
if plot_spectrum_flag == 1
    num_itr = 1;
end

MSE_ML_st_known = zeros(length(SNR_dB), length(true_AoA));
MSE_conventional = zeros(length(SNR_dB), length(true_AoA));
MSE_MVDR = zeros(length(SNR_dB), length(true_AoA));
MSE_MUSIC = zeros(length(SNR_dB), length(true_AoA));

tic
for idx_true_AoA = 1: length(true_AoA)
    if debug_flag == 1
        disp(['------------------------------------------  True AoA = ', num2str(true_AoA(idx_true_AoA)), ' -----------------------------------------------------'])
    end
    for idx_snr = 1: length(SNR_dB)
        sigma_n2 = E / 10^(SNR_dB(idx_snr)/10); % noise variance
        amplitude_gain = E_amplitude_gain * exp(1i*unifrnd(0,2*pi));   % channel apmlitude gain, assume it's a random phase shift
        
        % Received Signal
        channel = Channel(lambda, distance, sigma_n2, element_num, true_AoA(idx_true_AoA), ula_spacing);
        
        % Friis loss
        %y_t_friis = channel.FriisModel(s_t);
        % LoS with specified loss
        y_t_los = channel.LoS_specified_Loss(s_t, amplitude_gain);
        % After channel steering vector
        y_t_ula = channel.ULA(y_t_los);
        
        if plot_spectrum_flag == 1
            figure;
            hold on
        end
        for idx_iter = 1: num_itr
            % -------------------------------------------- AWGN -------------------------------------
            y_t_ula_awgn = channel.AWGN(y_t_ula);
            
            % ------------------------- ML assume the trasnmitted signal and the channel gain are known --------------
            % ------ angle spectrum
            ML_spectrum = zeros(size(angles));
            for i = 1: length(angles)
                ML_spectrum(i) = real(amplitude_gain * s_t * y_t_ula_awgn' * steering_vectors(:,i));
            end
            % ------- estimation
            [spectrum_value_max_ML, idx_max] = max(ML_spectrum);
            est_ML = angles(idx_max);
            MSE_ML_st_known(idx_snr, idx_true_AoA) = MSE_ML_st_known(idx_snr, idx_true_AoA) + (est_ML - true_AoA(idx_true_AoA))^2;
            
            if debug_flag == 1
                disp(['ML with known s(t) and channel gain, Estimated AoA = ', num2str(est_ML)]);
            end
            if plot_spectrum_flag == 1 % plotting the spectrum
                plot(angles, ML_spectrum/spectrum_value_max_ML , plot_colours{1}, 'DisplayName', 'ML - known s(t)');   % normalise the spectrum 
                scatter(est_ML, 1, 64, plot_colours{1}, 'DisplayName', 'ML estimation');
            end
            
            
            % ----------------------------- Conventional beamscan -----------------
            % --- Calculate the cov matrix
            R_conventional = y_t_ula_awgn * y_t_ula_awgn';
            % ------ angle spectrum
            conventional_spectrum = zeros(size(angles));
            for i = 1:length(angles)
                %steering_vector = exp(-1j * 2 * pi * ula_spacing * (0:element_num-1)' * sind(angles(i)) / lambda);
                conventional_spectrum(i) = steering_vectors(:,i)' * R_conventional * steering_vectors(:,i);
            end
            conventional_spectrum = abs(conventional_spectrum);
            % ------- estimation
            [spectrum_value_max_conventional, idx_max] = max(conventional_spectrum);
            est_conventional = angles(idx_max);
            MSE_conventional(idx_snr, idx_true_AoA) = MSE_conventional(idx_snr, idx_true_AoA) + (est_conventional - true_AoA(idx_true_AoA))^2;
            
            if debug_flag == 1
                disp(['Conventional beamscan, Estimated AoA = ', num2str(est_conventional)]);
            end
            if plot_spectrum_flag == 1 % plotting the spectrum
                plot(angles, conventional_spectrum/spectrum_value_max_conventional , plot_colours{2}, 'DisplayName', 'Conventional beamscan');   % normalise the spectrum 
                scatter(est_conventional, 1, 64, plot_colours{2}, 'DisplayName', 'Conventional beamscan estimation');
            end
            
            
%             % ----------------------------- MVDR Capon --------------------------------
%             % --- Calculate the cov matrix
%             R_mvdr = y_t_ula_awgn * y_t_ula_awgn' + eps(1)*eye(4); % add a small positive constant
%             % Compute the MVDR spectrum
%             mvdr_spectrum = zeros(size(angles));
%             for i = 1:length(angles)
%                 % steering_vector = exp(-1j * 2 * pi * ula_spacing * (0:element_num-1)' * sind(angles(i)) / lambda);
%                 mvdr_spectrum(i) = 1 / (steering_vectors(:,i)' * inv(R_mvdr) * steering_vectors(:,i));
%             end
%             mvdr_spectrum = abs(mvdr_spectrum);
%             % ------- estimation
%             [spectrum_value_max_mvdr, idx_max] = max(mvdr_spectrum);
%             est_mvdr = angles(idx_max);
%             MSE_MVDR(idx_snr, idx_true_AoA) = MSE_MVDR(idx_snr, idx_true_AoA) + (est_mvdr - true_AoA(idx_true_AoA))^2;
%             
%             if debug_flag == 1
%                 disp(['MVDR, Estimated AoA = ', num2str(est_mvdr)]);
%             end
%             if plot_spectrum_flag == 1 % plotting the spectrum
%                 plot(angles, mvdr_spectrum/spectrum_value_max_mvdr , plot_colours{3}, 'DisplayName', 'MVDR');   % normalise the spectrum 
%                 scatter(est_mvdr, 1, 64, plot_colours{3}, 'DisplayName', 'MVDR estimation');
%             end
%             
%             
%             % ---------------------------- MUSIC Algorithm -------------------------------
%             % --- Calculate the corr matrix
%             R_music = y_t_ula_awgn * y_t_ula_awgn' / size(y_t_ula_awgn, 2);
%             % Perform eigenvalue decomposition
%             [eigenvectors, eigenvalues] = eig(R_music);
%             % Sort eigenvalues and eigenvectors
%             [eigenvalues, idx_max] = sort(diag(eigenvalues), 'descend');
%             eigenvectors = eigenvectors(:, idx_max);
%             % Determine the noise subspace
%             num_signals = num_txs;
%             noise_subspace = eigenvectors(:, num_signals+1:end);
%             % Compute the MUSIC spectrum
%             music_spectrum = zeros(size(angles));
%             for i = 1:length(angles)
%                 %steering_vector = exp(-1j * 2 * pi * ula_spacing * (0:element_num-1)' * sind(angles(i)) / lambda);
%                 music_spectrum(i) = 1 / (steering_vectors(:,i)' * (noise_subspace * noise_subspace') * steering_vectors(:,i) + eps(1)); % add a small positive constant to prevent division by zero. 9.44 in [1]
%             end
%             music_spectrum = abs(music_spectrum);
%             % ------- estimation
%             [spectrum_value_max_music, idx_max] = max(music_spectrum);
%             est_music = angles(idx_max);
%             MSE_MUSIC(idx_snr, idx_true_AoA) = MSE_MUSIC(idx_snr, idx_true_AoA) + (est_music - true_AoA(idx_true_AoA))^2;
%             
%             if debug_flag == 1
%                 disp(['MUSIC, Estimated AoA = ', num2str(est_music)]);
%             end
%             if plot_spectrum_flag == 1 % plotting the spectrum
%                 plot(angles, music_spectrum/spectrum_value_max_music, plot_colours{4}, 'DisplayName', 'MUSIC');   % normalise the spectrum 
%                 scatter(est_music, 1, 64, plot_colours{4}, 'DisplayName', 'MUSIC estimation');
%             end
        end
    end
    if plot_spectrum_flag == 1
        hold off
        grid on
        title(['SNR = ', num2str(SNR_dB(idx_snr)), 'dB, AoA = ', num2str(true_AoA(idx_true_AoA))]);
        xlabel('Angle of interest, degrees')
        ylabel('Angle spectrum')
        legend
    end
end
runtime = toc

MSE_ML_st_known = MSE_ML_st_known / num_itr;
MSE_conventional = MSE_conventional / num_itr;
MSE_MVDR = MSE_MVDR / num_itr;
MSE_MUSIC = MSE_MUSIC / num_itr;

%% calculate CRLB
saclar_part = 8 * pi^2 * ula_spacing^2 / lambda^2 * element_num * (element_num-1) * (2*element_num-1) / 6;
angle_term = cosd(true_AoA).^2;
CRLB = zeros(length(SNR_dB), length(true_AoA));
for idx_snr = 1: length(SNR_dB)
    CRLB(idx_snr, :) = 1 ./ (10^(SNR_dB(idx_snr)/10) * saclar_part * angle_term);
end

%% plot MSE vs. angles
% Plot angle spectrums
figure;
for idx_snr = 1: length(SNR_dB)
    subplot(length(SNR_dB),1,idx_snr)
    semilogy(true_AoA, MSE_ML_st_known(idx_snr, :), plot_colours{1}, 'DisplayName', 'ML - known s(t)');
    hold on
    semilogy(true_AoA, MSE_conventional(idx_snr, :), plot_colours{2}, 'DisplayName', 'Conventional');
    semilogy(true_AoA, MSE_MVDR(idx_snr, :), plot_colours{3}, 'DisplayName', 'MVDR');
    semilogy(true_AoA, MSE_MUSIC(idx_snr, :), plot_colours{4}, 'DisplayName', 'MUSIC');
    % semilogy(true_AoA, CRLB(idx_snr, :), 'DisplayName', 'CRLB');
    grid on
    xlabel('AoA - ground truth')
    ylabel('MSE')
    title(['SNR = ', num2str(SNR_dB(idx_snr)), 'dB'])
end
legend
