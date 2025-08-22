clear; clc;

SNR_dB = 0:2:50;
SNR_lin = 10.^(SNR_dB/10);
num_symbols = 1e3; % symbs per packet
num_iterations = 1e4;

BER_1Tx1Rx = zeros(size(SNR_dB));
BER_1Tx2Rx = zeros(size(SNR_dB));
BER_1Tx4Rx = zeros(size(SNR_dB));

w = waitbar(0, 'Running Simulation...');

L = 3; % Number of paths

for idx = 1:length(SNR_dB)
    total_errors_1Tx1Rx = 0;
    total_errors_1Tx2Rx = 0;
    total_errors_1Tx4Rx = 0;
    total_symbols = 0;
    
    noise_var = 1 / (2 * SNR_lin(idx));
    
    for iter = 1:num_iterations
        tx_symbols = randi([0, 1], num_symbols, 1) * 2 - 1; % BPSK

        H1 = multipathRayleighFadingEnvelope(num_symbols, 100, L);
        H21 = multipathRayleighFadingEnvelope(num_symbols, 100, L);
        H22 = multipathRayleighFadingEnvelope(num_symbols, 100, L);
        H = [multipathRayleighFadingEnvelope(num_symbols, 100, L); multipathRayleighFadingEnvelope(num_symbols, 100, 3); multipathRayleighFadingEnvelope(num_symbols, 100, 3); multipathRayleighFadingEnvelope(num_symbols, 100, 3)];
        
        for n = 1:num_symbols
            % No Diversity (1Tx, 1Rx)
            h1 = H1(n);
            noise = sqrt(noise_var) * (randn + 1j*randn);
            y1 = h1 * tx_symbols(n) + noise;
            rx_decoded_1Tx1Rx = real(y1 / h1) > 0;
            total_errors_1Tx1Rx = total_errors_1Tx1Rx + (rx_decoded_1Tx1Rx ~= (tx_symbols(n) > 0));
            
            % MRRC (1Tx, 2Rx)
            h21 = H21(n);
            h22 = H22(n);
            y21 = h21 * tx_symbols(n) + sqrt(noise_var) * (randn + 1j*randn);
            y22 = h22 * tx_symbols(n) + sqrt(noise_var) * (randn + 1j*randn);
            rx_combined_1Tx2Rx = (conj(h21) * y21 + conj(h22) * y22) / (abs(h21)^2 + abs(h22)^2);
            rx_decoded_1Tx2Rx = real(rx_combined_1Tx2Rx) > 0;
            total_errors_1Tx2Rx = total_errors_1Tx2Rx + (rx_decoded_1Tx2Rx ~= (tx_symbols(n) > 0));
            
            % MRRC (1Tx, 4Rx)
            h = H(:, n).';
            y = h .* tx_symbols(n) + sqrt(noise_var) * (randn(1,4) + 1j*randn(1,4));
            rx_combined_1Tx4Rx = sum(conj(h) .* y) / sum(abs(h).^2);
            rx_decoded_1Tx4Rx = real(rx_combined_1Tx4Rx) > 0;
            total_errors_1Tx4Rx = total_errors_1Tx4Rx + (rx_decoded_1Tx4Rx ~= (tx_symbols(n) > 0));
        end
        
        total_symbols = total_symbols + num_symbols;
    end
    
    % Compute BER
    BER_1Tx1Rx(idx) = total_errors_1Tx1Rx / total_symbols;
    BER_1Tx2Rx(idx) = total_errors_1Tx2Rx / total_symbols;
    BER_1Tx4Rx(idx) = total_errors_1Tx4Rx / total_symbols;
    
    waitbar(idx / length(SNR_dB), w);
end

close(w);

% simulation Parameters
N = 1e3;                      % num symbs
SNR_dB = 0:2:50;              % SNR range in dB
SNR_lin = 10.^(SNR_dB/10);    % linear SNR values
num_iterations = 1e4;        

BER_2Tx1Rx = zeros(length(SNR_dB), num_iterations);
BER_2Tx2Rx = zeros(length(SNR_dB), num_iterations);

h = waitbar(0, 'Running Simulation...');

for iter = 1:num_iterations
    waitbar(iter / num_iterations, h);
    tx_symbols = 2*randi([0,1], 1, N) - 1; % BPSK
    for idx = 1:length(SNR_dB)
        num_pairs = floor(N/2);
        rx_combined_2Tx1Rx = zeros(1, N);
        rx_combined_2Tx2Rx = zeros(1, N);

        H21 = multipathRayleighFadingEnvelope(num_pairs, 100, L); 
        H22 = multipathRayleighFadingEnvelope(num_pairs, 100, L);
        H31 = multipathRayleighFadingEnvelope(num_pairs, 100, L);
        H32 = multipathRayleighFadingEnvelope(num_pairs, 100, L);
        
        for p = 1:num_pairs
            n = 2*p-1; %odd index

            h21 = H21(p); h22 = H22(p);  % 2Tx/1Rx
            h31 = H31(p); h32 = H32(p);  % 2Tx/2Rx
            
            noise_var = 1/(2*SNR_lin(idx));
            
            % Alamouti (2Tx/1Rx)
            s1 = tx_symbols(n)/sqrt(2); % scaled the power of the tx antennas to equal the total power of one
            s2 = tx_symbols(n+1)/sqrt(2);
            y11 = h21*s1 + h22*s2 + sqrt(noise_var) * (randn + 1j*randn);
            y12 = -h21*conj(s2) + h22*conj(s1) + sqrt(noise_var) * (randn + 1j*randn);
            s1_est_2Tx1Rx = conj(h21) * y11 + h22 * conj(y12);
            s2_est_2Tx1Rx = conj(h22) * y11 - h21 * conj(y12);
            h_norm_2Tx1Rx = abs(h21)^2 + abs(h22)^2;
            rx_combined_2Tx1Rx(n) = s1_est_2Tx1Rx / h_norm_2Tx1Rx;
            rx_combined_2Tx1Rx(n+1) = s2_est_2Tx1Rx / h_norm_2Tx1Rx;
            
            % Alamouti (2Tx/2Rx)
            y21 = h31*s1 + h32*s2 + sqrt(noise_var) * (randn + 1j*randn);
            y22 = -h31*conj(s2) + h32*conj(s1) + sqrt(noise_var) * (randn + 1j*randn);
            s1_est_2Tx2Rx = conj(h21) * y11 + h22 * conj(y12) + conj(h31) * y21 + h32 * conj(y22);
            s2_est_2Tx2Rx = conj(h22) * y11 - h21 * conj(y12) + conj(h32) * y21 - h31 * conj(y22);
            h_norm_2Tx2Rx = abs(h21)^2 + abs(h22)^2 + abs(h31)^2 + abs(h32)^2;
            rx_combined_2Tx2Rx(n) = s1_est_2Tx2Rx / h_norm_2Tx2Rx;
            rx_combined_2Tx2Rx(n+1) = s2_est_2Tx2Rx / h_norm_2Tx2Rx;
        end
        
        tx_bits = (tx_symbols + 1) / 2; 
        BER_2Tx1Rx(idx, iter) = sum(real(rx_combined_2Tx1Rx) > 0 ~= tx_bits) / N;
        BER_2Tx2Rx(idx, iter) = sum(real(rx_combined_2Tx2Rx) > 0 ~= tx_bits) / N;
    end
end

close(h);

% compute average BER over iterations
BER_2Tx1Rx_avg = mean(BER_2Tx1Rx, 2);
BER_2Tx2Rx_avg = mean(BER_2Tx2Rx, 2);

% Plot all of rhe results
figure;
semilogy(SNR_dB, BER_1Tx1Rx, 'r-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB, BER_1Tx2Rx, 'b-s', 'LineWidth', 2);
semilogy(SNR_dB, BER_1Tx4Rx, 'g-d', 'LineWidth', 2);
semilogy(SNR_dB, BER_2Tx1Rx_avg, 'x-', 'LineWidth', 2);
semilogy(SNR_dB, BER_2Tx2Rx_avg, '*-', 'LineWidth', 2);
legend('1Tx-1Rx (No Diversity)', '1Tx-2Rx (MRRC)', '1Tx-4Rx (MRRC)', 'Alamouti (2Tx/1Rx)', 'Alamouti (2Tx/2Rx)');
xlabel('SNR (dB)'); ylabel('BER');
grid on;
title('BER Performance of Diversity Schemes Simulation');

function h = multipathRayleighFadingEnvelope(N, fm, L)
    H = zeros(L, N);
    theta = 2*pi*rand(1,L); % introduce random phase delays represent time delays
    for l=1:L
        H(l,:) = exp(-1j*theta(l)) * rayleighFadingEnvelope(N, fm); % put it across the entire envelope
    end
    h = sum(H);
end

function h = rayleighFadingEnvelope(N, fm)
    delta_f = 2*fm / (N-1);
    T = 1/delta_f;
    f = -fm;
    freqs = linspace(-fm, fm, N);
    
    i = normrnd(0, 1, 1, N/2) + 1j*normrnd(0, 1, 1, N/2);
    q = normrnd(0, 1, 1, N/2) + 1j*normrnd(0, 1, 1, N/2);
    
    i = cat(2, conj(flip(i)), i);
    q = cat(2, conj(flip(q)), q);
    
    S = 1.5 ./ (pi*fm*sqrt(1 - ((freqs) / fm).^2));
    S(1) = S(2); S(N) = S(N-1); % prevent singularities
    
    i = i .* sqrt(S);
    q = q .* sqrt(S);
    
    x = ifft(i, 'symmetric');
    y = ifft(q, 'symmetric');
    h = sqrt(x.^2 + y.^2);
    h = 2 * (h - min(h)) / (max(h) - min(h)) - 1;
end