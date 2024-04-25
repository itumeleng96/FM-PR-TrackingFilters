% Define the range of N (number of cells)
N = 1:40;

% Define the range of P_FA (false alarm probability)
P_FA = [1e-8,1e-6,1e-4,0.01, 0.05, 0.1, 0.2]; % Example values of P_FA

% Initialize a matrix to store alpha values for different P_FA
alpha_values = zeros(length(N), length(P_FA));

% Calculate alpha for each combination of N and P_FA
for i = 1:length(P_FA)
    alpha_values(:, i) = N .* (P_FA(i).^(-1./N) - 1);
end

% Plot alpha against N for different values of P_FA
figure;
plot(N, alpha_values);
legend(cellstr(num2str(P_FA', 'P_{FA} = %0.1e')));
xlabel('Number of cells (N)');
ylabel('Threshold factor (\alpha)');
title('Threshold factor (\alpha) vs. Number of cells (N) for different P_{FA}');
grid on;
ylim([0 100]);