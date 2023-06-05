function z = expDdecay(data, tau_f, tau_s)
% Uses Meucci's double exponential decay method 
% (based on Holt's linear method) to smooth data
% to reduce noise. Follows on to normalise the data for comparative
% purposes. 
% NOTE: can be used in an external 'for' loop of the columns of a timetable
% to smooth all columns.

% INPUTS:
% data = unscaled state varaible data eg. Volitility index, inflation etc. 
% (type: double)
% tau_f = fast half-life at which decay volatilities
% (type: double)
% tau_s = slow half-life at which decay correlations
% (type: double)

    T = length(data);
    timestamp = 1:T;

    % Smoothing State Variable
    state_sig = data';
    
    % row of 0s
    z = zeros(1, T);

    for t = 1:T
        p_es_fast = exp(-log(2) / tau_f * (t - timestamp(1:t)));
        gamma_f = sum(p_es_fast);
        z(t) = sum(p_es_fast .* state_sig(1:t)) / gamma_f;
    end

    % Scoring State Variable
    mu_est = zeros(1, T);
    mu2_est = zeros(1, T);
    sd_est = zeros(1, T);

    for t = 1:T
        p_es_slow = exp(-log(2) / tau_s * (t - timestamp(1:t)));
        gamma_s = sum(p_es_slow);

        mu_est(t) = sum(p_es_slow .* z(1:t)) / gamma_s;
        mu2_est(t) = sum(p_es_slow .* z(1:t).^2) / gamma_s;
        sd_est(t) = sqrt(mu2_est(t) - mu_est(t)^2);
    end

    z = ((z - mu_est) ./ sd_est)';
    
end
