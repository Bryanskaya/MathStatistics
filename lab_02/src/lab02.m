function lab_02()
    x = 
	[-4.58,-5.10,-4.24,-4.82,-6.05,-4.05,-4.48,-4.65,-3.67,-4.01,-3.22,
	-5.79,-4.20,-5.36,-5.16,-4.31,-3.91,-4.04,-5.65,-5.78,-4.03,-5.15,
	-4.89,-4.60,-5.71,-5.67,-4.40,-5.10,-5.29,-3.90,-3.77,-5.47,-5.57,
	-5.91,-3.62,-2.45,-3.13,-5.60,-4.35,-3.36,-3.87,-4.78,-5.72,-4.66,
	-4.34,-4.60,-5.24,-4.43,-5.15,-5.45,-3.22,-4.61,-5.65,-5.63,-4.73,
	-3.82,-4.34,-4.98,-6.43,-4.25,-4.66,-5.49,-4.98,-6.10,-4.44,-5.25,
	-4.89,-2.97,-7.26,-4.14,-4.45,-5.93,-3.49,-3.96,-5.66,-5.04,-5.40,
	-4.89,-4.61,-4.01,-5.44,-2.76,-4.97,-3.87,-4.33,-5.28,-5.63,-3.94,
	-4.56,-4.67,-5.35,-4.89,-3.79,-6.25,-4.38,-4.28,-4.48,-4.16,-5.67,
	-5.34,-4.78,-4.33,-4.89,-3.16,-4.88,-5.40,-4.64,-4.84,-6.72,-5.25,
	-3.27,-6.44,-3.49,-5.92,-5.21,-3.88,-5.08,-3.50,-6.80,-5.87];
    
    gamma = 0.9;
    
    % Задание 1а - вычисление точечных оценок µˆ и S2 для MX и DX
    mu = find_mu(x);
    fprintf('mu = %.2f\n', mu);
    
    S2 = find_S2(x);
    fprintf('S^2 = %.2f\n', S2);
    
    % Задание 1б - вычисление нижней и верхней границ доверительного
    % интервала для МХ
    mu_low = find_mu_low(x, gamma);
    S2_high = find_mu_high(x, gamma);
    fprintf("%.2f < MX < %.2f\n", mu_low, S2_high);
    
    % Задание 1в - вычисление нижней и верхней границ доверительного
    % интервала для DХ
    S2_low = find_S2_low(x, gamma);
    S2_high = find_S2_high(x, gamma);
    fprintf("%.2f < DX < %.2f\n", S2_low, S2_high);
    
    % Задание 3а
    task_mu(x, gamma);
    task_S2(x, gamma); 
end 

function res = find_mu(x)
    res = sum(x) / length(x);
end

function res = find_S2(x)
    res = sum((x - find_mu(x)).^2) / (length(x) - 1);
end

function res = find_mu_low(x, gamma)
    n = length(x);
    mu = find_mu(x);
    S2 = find_S2(x);
    res = mu - sqrt(S2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function res = find_mu_high(x, gamma)
    n = length(x);
    mu = find_mu(x);
    S2 = find_S2(x);
    res = mu + sqrt(S2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function res = find_S2_low(x, gamma)
    n = length(x);
    S2 = find_S2(x);
    res = ((n - 1) * S2) / chi2inv((1 + gamma) / 2, n - 1);
end

function res = find_S2_high(x, gamma)
    n = length(x);
    S2 = find_S2(x);
    res = ((n - 1) * S2) / chi2inv((1 - gamma) / 2, n - 1);
end

function task_mu(x, gamma)
    figure;
    grid on;
    hold on;
    
    n = length(x);
    mu = find_mu(x);
    n_range = 1:1:n;
    
    plot([n_range(1), n_range(n)], [mu, mu], 'b');
    
    mu_array = zeros(n);
    mu_low_array = zeros(n);
    mu_high_array = zeros(n);
    
    for i=1:n
        mu_array(i) = find_mu(x(1:i));
        mu_low_array(i) = find_mu_low(x(1:i), gamma);
        mu_high_array(i) = find_mu_high(x(1:i), gamma);
    end
    
    plot(n_range, mu_array, 'g--');
    plot(n_range, mu_low_array, 'r-');
    plot(n_range, mu_high_array, 'm:', 'Linewidth',2);
    
    xlabel('n');
	ylabel('y');
    legend('$\hat{\mu}^2(\vec x_N)$', ...
        '$\hat{\mu}^2(\vec x_n)$', ...
        '$\underline{\mu}^2(\vec x_n)$', ...
        '$\overline{\mu}^2(\vec x_n)$', ...
        'Interpreter', 'latex');
end

function task_S2(x, gamma)
    figure;
    grid on;
    hold on;
    
    n = length(x);
    S2 = find_S2(x);
    n_range = 1:1:n;
    
    plot([n_range(1), n_range(n)], [S2, S2], 'b');
    
    S2_array = zeros(n);
    S2_low_array = zeros(n);
    S2_high_array = zeros(n);
    
    for i=1:n
        S2_array(i) = find_S2(x(1:i));
        S2_low_array(i) = find_S2_low(x(1:i), gamma);
        S2_high_array(i) = find_S2_high(x(1:i), gamma);
    end
    
    plot(n_range, S2_array, 'g--');
    plot(n_range, S2_low_array, 'r-');
    plot(n_range, S2_high_array, 'm:', 'Linewidth',2);
    
    legend('$\hat S^2(\vec x_N)$', ...
        '$\hat S^2(\vec x_n)$', ...
        '$\underline{\sigma}^2(\vec x_n)$', ...
        '$\overline{\sigma}^2(\vec x_n)$', ...
        'Interpreter', 'latex');
end