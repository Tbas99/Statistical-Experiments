% We are given m coins. m - 1 coins are fair, but there is one coin which
% is fake and flips heads with probability 1/2 + ε, where ε > 0. However, 
% we don't know the value of ε. We would like to identify the fake coin 
% quickly.
% Solution: Consider the following strategy. Flip both coins n times and
% compute the proportion of heads of each coin (say p1 and p2 for coins 1
% and 2, respectively). Now deem the coin for which the proportion of heads
% is larger to be the fake coin. What is the probability we'll make a
% mistake? We implement the above procedure with ε = 0.1, m = 2 and the 
% following values of n = 10, 100, 500, 1000. For each choice of parameters
% m and n repeat the procedure N = 10000 times and compute the proportion
% of runs where the procedure failed to identify the correct coin. Finally,
% we compare the proportion found with the theoretical (hoeffding
% inequallity) bound.

% Clear workspace
clear 

% General parameters
nr_runs = 10000;
show_plots = 0;

% Procedure parameters
epsilon = 0.1;
m = 2;

% Experiment
for n = [10, 100, 500, 1000]
    % Define error bounds (equivalent for m = 2)
    error_bound = (m - 1)*exp(-n*(epsilon)^2);

    results = zeros((m - 1), nr_runs);
    for run = 1:nr_runs
        % Flip the coins n times and record observations
        fake_coin_tosses = zeros(1, n);
        real_coin_tosses = zeros((m - 1), n);
        for i = 1:n
            % Flip fake coin 
            fake_coin_tosses(i) = flip_coin(epsilon);
    
            % Flip all other coins
            for j = 1:(m - 1)
                real_coin_tosses(j, i) = flip_coin(0);
            end
        end
        
        % Visualize tosses
        if (show_plots)
            figure('Name', sprintf('Flips for fake coin for n = %d', n))
            h1 = histogram(fake_coin_tosses, 2);
            h1.FaceAlpha = 0.3;
            h1.FaceColor = "red";
    
            % Visualize first real coin only
            figure('Name', sprintf('Flips for real coin for n = %d', n))
            h2 = histogram(real_coin_tosses(1,:), 2);
            h2.FaceAlpha = 0.3;
            h2.FaceColor = "green";
        end
    
        % Compute proportion head throws
        fake_coin_head_proportion = sum(fake_coin_tosses == 1) / n;
        real_coin_head_proportion = zeros((m - 1), 1);
        for j = 1:(m - 1)
            real_coin_head_proportion(j, 1) = sum(real_coin_tosses(j, :) == 1) / n;
        end
    
    
        % Return results
        % 0 = Correctly picked the fake coin
        % 1 = Incorrectly picked the fake coin
        for j = 1:(m - 1)
            results(j, run) = fake_coin_head_proportion < real_coin_head_proportion(j, 1);
        end
    end

    % Compute failure proportion
    failure_proportion = zeros((m - 1), 1);
    for j = 1:(m - 1)
        failure_proportion(j, 1) = sum(results(j, :) == 1) / nr_runs;
    end

    % Compare failure proportion with bound
    sprintf(['Failure proportion for n = %d -> %f \n' ...
        'Hoeffding bound for n = %d -> %f'], n, failure_proportion, n, error_bound)
end



% Returns head = 1 or tail = 0
% Parameter bias: Positive towards head, Negative towards tail
function result = flip_coin(bias)
    result = random('Binomial', 1, (1/2) + bias);
end

