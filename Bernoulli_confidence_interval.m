% Let X1, ..., Xn be independent Bernoulli random variables with parameter
% p, Let p^_n (p hat n) = (1/n)*sum(Xi) for i=1 to n, and define the
% intervals [p^_n - sqrt(log(2/delta)/2n), p^_n + sqrt(log(2/delta)/2n)]
% with delta > 0. For every particular value of n the corresponding
% interval is guaranteed to contain the true parameter p with probability
% at least 1 − delta. However, this does not mean the probability all these
% intervals simultaneously contain the true parameter is larger than
% 1 − delta. We experimentally verify this.

% Clear workspace
clear 

% General parameters
nr_runs = 100;
show_plots = 0;

% Procedure parameters
n = 10000;
p = 1/2;
delta = 0.1;

% Experiment
intervals_contain_p = zeros(nr_runs, 1);
intervals_modified_contain_p = zeros(nr_runs, 1);
intervals_itlog_contain_p = zeros(nr_runs, 1);
for run = 1:nr_runs
    % Generate 10000 bernoulli samples
    x = zeros(1, n);
    for i = 1:n
        x(i) = random('Binomial', 1, p);
    end
    
    % Show sample distribution
    if (show_plots)
        figure('Name', sprintf('Bernoulli samples for n = %d', n))
        h1 = histogram(x, 2);
    end
    
    % Generate intervals for n_i = 3...n
    intervals = zeros((n - 2), 3);
    intervals_modified = zeros((n - 2), 3);
    intervals_itlog = zeros((n - 2), 3);
    for n_i = 3:n
        % Consider only the observations in range 1...n_i
        % since we restrict ourselves to n_i = 3...n
        x_for_interval = x(1:n_i);
        p_hat_n = (1/n_i)*sum(x_for_interval);
    
        % Compute and store interval endpoints
        lower_endpoint = p_hat_n - sqrt(log(2/delta)/(2*n_i));
        upper_endpoint = p_hat_n + sqrt(log(2/delta)/(2*n_i));
        intervals((n_i - 2), :) = [n_i, lower_endpoint, upper_endpoint];

        % b. Using different endpoints
        lower_endpoint = p_hat_n - ...
            sqrt((log(n_i*(n_i + 1)) + log(2/delta))/(2*n_i));
        upper_endpoint = p_hat_n + ...
            sqrt((log(n_i*(n_i + 1)) + log(2/delta))/(2*n_i));
        intervals_modified((n_i - 2), :) = ...
            [n_i, lower_endpoint, upper_endpoint];

        % c. Interval inspired by the law of the iterated logarithm
        lower_endpoint = p_hat_n - ...
            sqrt((log(log(n_i)) + log(2/delta))/(2*n_i));
        upper_endpoint = p_hat_n + ...
            sqrt((log(log(n_i)) + log(2/delta))/(2*n_i));
        intervals_itlog((n_i - 2), :) = ...
            [n_i, lower_endpoint, upper_endpoint];
    end
    
    % Check if p is contained in the intervals
    p_contained = zeros((n - 2), 1);
    p_contained_in_modified_intervals = zeros((n - 2), 1);
    p_contained_in_itlog_intervals = zeros((n - 2), 1);
    for n_i = 3:n
        p_contained(n_i - 2) = ...
            (intervals((n_i - 2), 2) <= p) && ...
            (p <= intervals((n_i - 2), 3));

        p_contained_in_modified_intervals(n_i - 2) = ...
            (intervals_modified((n_i - 2), 2) <= p) && ...
            (p <= intervals_modified((n_i - 2), 3));

        p_contained_in_itlog_intervals(n_i - 2) = ...
            (intervals_itlog((n_i - 2), 2) <= p) && ...
            (p <= intervals_itlog((n_i - 2), 3));
    end
    intervals_contain_p(run, 1) = all(p_contained);
    intervals_modified_contain_p(run, 1) = ...
        all(p_contained_in_modified_intervals);
    intervals_itlog_contain_p(run, 1) = ...
        all(p_contained_in_itlog_intervals);
end

% Get proportion of p contained in intervals and compare
proportion_p_in_intervals = sum(intervals_contain_p) / nr_runs;
probability_guarantee = 1 - delta;

% Indication; very likely, the probability all the intervals will contain
% 1/2 is actually smaller than 1 − δ
fprintf(['Results over %d runs: \n' ...
    'Proportion p contained in intervals: %f \n' ...
    'Probability Guarantee: %f \n'], nr_runs, ...
    proportion_p_in_intervals, probability_guarantee)

% With modified intervals =>
proportion_p_in_modified_intervals = ...
    sum(intervals_modified_contain_p) / nr_runs;
fprintf(['Results over %d runs: \n' ...
    'Proportion p contained in modified intervals: %f \n' ...
    'Probability Guarantee: %f \n'], nr_runs, ...
    proportion_p_in_modified_intervals, probability_guarantee)

% law of the iterated logarithm intervals =>
proportion_p_in_itlog_intervals = ...
    sum(intervals_itlog_contain_p) / nr_runs;
fprintf(['Results over %d runs: \n' ...
    'Proportion p contained in iterated logarithm intervals: %f \n' ...
    'Probability Guarantee: %f \n'], nr_runs, ...
    proportion_p_in_itlog_intervals, probability_guarantee)

