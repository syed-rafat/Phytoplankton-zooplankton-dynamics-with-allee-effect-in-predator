%% Classic rosenzweig-macarthur

K_values = [1.3, 2.4, 2.6];

for i = 1: length(K_values)
    noiseless_timeseries_soln = [];
    use classic predator prey.ini
    simtime 0 2000 2000
    solver -n
    K=K_values(i);
    % Run until positive solution is found
    noiseless_timeseries = time('-s');
    noiseless_timeseries_soln = noiseless_timeseries;
    noiseless_timeseries_sliced = noiseless_timeseries_soln(1:2000, :);
    
    
    plot(noiseless_timeseries_sliced(:,1), noiseless_timeseries_sliced(:,2), 'Color', 'a7c957', 'LineWidth', 2)
    hold on
    plot(noiseless_timeseries_sliced(:,1), noiseless_timeseries_sliced(:,3), 'Color', '#023047', 'LineWidth', 2)
    % h = legend({'P','Z'}, 'FontSize', 24, 'FontWeight','bold');
    xlabel('time');
    ylabel('population');
    yl = ylim;
    text(600, 0.9*yl(2), ['K= ', num2str(K)], 'FontSize', 24)
    hold off;
    
    saveas(gcf, ['figures\classic\K_',  num2str(K),'.png']) % gcf gets the currernt figure


    noisy_timeseries_soln = [];
    use classic predator prey with noise.ini
    simtime 0 2000 2000
    solver -n
    K=K_values(i);
    % Run until positive solution is found
    noisy_timeseries = time('-s');
    noisy_timeseries_soln = noisy_timeseries;
    noisy_timeseries_sliced = noisy_timeseries_soln(1:2000, :);
    
    
    plot(noisy_timeseries_sliced(:,1), noisy_timeseries_sliced(:,2), 'Color', '#a7c957', 'LineWidth', 2)
    hold on
    plot(noisy_timeseries_sliced(:,1), noisy_timeseries_sliced(:,3), 'Color', '#023047', 'LineWidth', 2)
    % h = legend({'P','Z'}, 'FontSize', 24, 'FontWeight','bold');
    xlabel('time');
    ylabel('population');
    yl = ylim;
    text(600, 0.9*yl(2), ['K= ', num2str(K)], 'FontSize', 24)
    hold off;
    
    saveas(gcf, ['figures\classic\K_',  num2str(K), 'sigma=0.05','.png'])

end

%% Model with Alle effect
%% Deterministic Allee effect

rng(25) % To reproduce the exact same results.

% Deterministic model
% use monod allee effect (no noise).ini
% simtime 0 5000 5000 %starts at 0 t, stops at 5000 days and only output 1 value per day
% K=2.7;
% noiseless = time();
% 
% noiseless(1:1000, :) = [];
% noiseless_phytoplankton_std = std(noiseless(:, 2));
% noiseless_zooplankton_std = std(noiseless(:, 3));
% 
% noiseless_phytoplankton_mean = mean(noiseless(:, 2));
% noiseless_zooplankton_mean = mean(noiseless(:, 3));
% % Calculate coefficient of variation
% cv_intrinsic_phyto = (noiseless_phytoplankton_std/noiseless_phytoplankton_mean) * 100;
% cv_intrinsic_zoo = (noiseless_zooplankton_std/noiseless_zooplankton_mean) * 100;


K1 = [1.4, 2.5, 2.8, 7.4];
K2 = [1.8, 2.5, 2.8, 6.5];
K3 = [2.2, 2.5, 2.8, 6.1];
K4 = [2.9, 3, 5, 5.5];
delta_values = [0.005, 0.01, 0.013, 0.015];

for j = 1: length(delta_values)
    switch j
        case 1
            K_values = K1;
        case 2
            K_values = K2;
        case 3
            K_values = K3;
        case 4
            K_values = K4;
    end
     
    for i = 1: length(K_values)
        noiseless_timeseries_soln = [];
        use monod allee effect (no noise).ini
        simtime 0 5000 5000
        solver -n
        K=K_values(i);
        delta=delta_values(j)
        % Run until positive solution is found
        noiseless_timeseries = time('-s');
        noiseless_timeseries_soln = noiseless_timeseries;
        
        % white_timeseries_soln(1:1000, :) = [];
        
        % white_timeseries_sliced = white_timeseries_soln(1001:2000, :);
        % adjusted_timevalues = white_timeseries_sliced(:,1) - 1000;
        % white_timeseries_sliced(:, 1) = adjusted_timevalues;
        noiseless_timeseries_sliced = noiseless_timeseries_soln(1:1000, :);
        
        
        
        plot(noiseless_timeseries_sliced(:,1), noiseless_timeseries_sliced(:,2), 'Color', '#a7c957', 'LineWidth', 2)
        hold on
        plot(noiseless_timeseries_sliced(:,1), noiseless_timeseries_sliced(:,3), 'Color', '#023047', 'LineWidth', 2)
        % h = legend({'P','Z'}, 'FontSize', 24, 'FontWeight','bold');
        xlabel('time');
        ylabel('population');
        yl = ylim;
        text(600, 0.9*yl(2), ['K= ', num2str(K)], 'FontSize', 24)
        hold off;
        
        saveas(gcf, ['figures\No noise\allee delta=', num2str(delta),' K=',  num2str(K),'.png'])

    end
end


%% Stochastic model (white noise)

% K_values = [1.8, 2.4, 2.7];
K_values = [1.7];
sigma_values = [0.05];

for noise_intensity = 1: length(sigma_values)
for i = 1: length(K_values)
white_timeseries_soln = [];
while true
    use monod allee effect with additive noise.ini
    simtime 0 5000 5000
    solver -n % 
    K=K_values(i);
    delta=0.01;
    sigma=sigma_values(noise_intensity);
    % Run until positive solution is found
    white_timeseries = time('-s');
    if ~isnan(white_timeseries)
        if  all(white_timeseries(:, 3) > 0.00000000001) % sometimes zooplankton soln has very small number
            white_timeseries_soln = white_timeseries;
        break
        end
    end
end

% white_timeseries_soln(1:1000, :) = [];

% white_timeseries_sliced = white_timeseries_soln(1001:2000, :);
% adjusted_timevalues = white_timeseries_sliced(:,1) - 1000;
% white_timeseries_sliced(:, 1) = adjusted_timevalues;
white_timeseries_sliced = white_timeseries_soln(1:1000, :);



plot(white_timeseries_sliced(:,1), white_timeseries_sliced(:,2), 'Color', '#a7c957', 'LineWidth', 2)
hold on
plot(white_timeseries_sliced(:,1), white_timeseries_sliced(:,3), 'Color', '#023047', 'LineWidth', 2)
xlabel('time');
ylabel('population');
yl = ylim;
text(600, 0.9*yl(2), ['K= ', num2str(K)], 'FontSize', 24)
hold off;

saveas(gcf, ['figures\white\whiteallee k_', num2str(K), ' sigma_', num2str(sigma), '.png'])


% noisy_phytoplankton_std = std(white_timeseries_soln(:, 2));
% noisy_zooplankton_std = std(white_timeseries_soln(:, 3));
% 
% noisy_phytoplankton_mean = mean(white_timeseries_soln(:, 2));
% noisy_zooplankton_mean = mean(white_timeseries_soln(:, 3));
% % Calculate coefficient of variation
% cv_noisy_phyto = (noisy_phytoplankton_std/noisy_phytoplankton_mean) * 100;
% cv_noisy_zoo = (noisy_zooplankton_std/noisy_zooplankton_mean) * 100;


end
end
%% Stochastic model(red noise)

% Generate red noise
% time_vector = 1:1:4000;
% lambda = 8;
% SD = 0.05;
% beta = SD * sqrt(2/lambda - 1/lambda^2); % Calculate beta for desired SD
% red = rednoise(time_vector, 0, lambda, beta);



% autocorrelation plot for the red noise
% figure;
% stem(autocorr(red));


% K_values = [2.3, 2.4, 2.5, 2.7, 3, 4];
sigma_values = [0.05];
K_values = [1.8, 2.5, 2.7];
for noise_intensity = 1: length(sigma_values)
for i = 1: length(K_values)
red_timeseries_soln = [];
while true % we need re initialize the model each time, otherwise while loop gets stuck
    use monod allee effect with red noise.ini
    simtime 0 5000 5000
    solver -n % solver euler is default for SDE, uses euler-maruyama
    K=K_values(i);
    delta=0.01;
    sigma=sigma_values(noise_intensity);
    % Run until positive solution is found
    red_timeseries = time('-s');
    if ~isnan(red_timeseries)
        if  all(red_timeseries(:, 3) > 0.00001)
            red_timeseries_soln = red_timeseries;
        break
        end
    end
end

% white_timeseries_soln(1:1000, :) = [];

% white_timeseries_sliced = white_timeseries_soln(1001:2000, :);
% adjusted_timevalues = white_timeseries_sliced(:,1) - 1000;
% white_timeseries_sliced(:, 1) = adjusted_timevalues;
red_timeseries_sliced = red_timeseries_soln(1:1000, :);



plot(red_timeseries_sliced(:,1), red_timeseries_sliced(:,2), 'Color', '#a7c957', 'LineWidth', 2)
hold on
plot(red_timeseries_sliced(:,1), red_timeseries_sliced(:,3), 'Color', '#023047', 'LineWidth', 2)
xlabel('time');
ylabel('population');
yl = ylim;
text(600, 0.9*yl(2), ['K= ', num2str(K)], 'FontSize', 24)
hold off;

saveas(gcf, ['figures\red\redallee sigma_', num2str(sigma), ' k_', num2str(K), '.png'])


% noisy_phytoplankton_std = std(white_timeseries_soln(:, 2));
% noisy_zooplankton_std = std(white_timeseries_soln(:, 3));
% 
% noisy_phytoplankton_mean = mean(white_timeseries_soln(:, 2));
% noisy_zooplankton_mean = mean(white_timeseries_soln(:, 3));
% % Calculate coefficient of variation
% cv_noisy_phyto = (noisy_phytoplankton_std/noisy_phytoplankton_mean) * 100;
% cv_noisy_zoo = (noisy_zooplankton_std/noisy_zooplankton_mean) * 100;


end
end