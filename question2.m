%% Question 2
% This section of the script generates constants to be used later in the script, including the number of
% particles, mass and probability of scattering. This section will also
% calculate an initial velocity by selecting from the Maxwell-Bolzman distribution
% and display a histogram in a subplot. Note: 
%
% The Boltzman distribution does not have enough particles to properly show it's form. Adding more
% particles makes the live plot too cluttered to see what is happening, so
% for the sake of clarity, 10 particles are being used for 50 timesteps. This can be
% modified freely.

num_particles = 10;     % best displayes what's happening without a mess
colours = ["r", "g", "b", "c", "y", "k", "m"];

kb = 1.38064852;
T = 300;
m = 9.10938356E-31;

mean_time_collision = 0.2;   % measured in picoseconds
mean_free_path = 0;          % initialized to zero in case no particles scatter immediately
timesteps = 1E-3;            % also in picoseconds, so this will be 1 femto second timesteps
pscat = 1 - exp(-timesteps/mean_time_collision);

% particle positions
particles = rand(num_particles, 2);
particles(:, 2) = particles(:, 2)*200;  % x-coordinates
particles(:, 1) = particles(:, 1)*100;  % y-coordinates

% select velocities from Maxwell-Boltzman distribution
% column 3 is x-velocity, 4 is x-velocity
particles(:, 3:4) = randn(num_particles, 2) *  sqrt(kb * T / m) / 1E15;  % scaled to femtoseconds
% subplot(2, 1, 1)  % subplots don't look good in publish
figure(1)
title(sprintf("Initial Avg Velocity: %s", mean(sqrt(particles(:, 3).^2 + particles(:, 4).^2))))
histogram(rssq(particles(:, 3:4), 50));

position_of_last_scatter = particles(:, 3:4);

%%
% The main time loop. Works very similar to quesiton 1 except it will also
% test if a particle scatters. If it does, it selects a new velocity at
% random from the Maxwell-Boltzman distribution. This section also checks the
% distance a particle has been able to travel between scatters and
% calculates the mean free path, which is shown alongside the average
% temperature in the plot title. Over time, the average temperature
% stabalizes at the centroid of the distribution because every time a
% particle is re-thermalized, it has a higher probability of taking on the
% mean velocity.

% main time loop
for i = 0:50      % each step is a femtosecond
    previous_particles = particles;
    
    % update positions
    particles(:, 1) = particles(:, 1) + particles(:, 4);
    particles(:, 2) = particles(:, 2) + particles(:, 3);
    
    % check if any particles passed the boundary and deal with them
    x_boundary_changes_right = particles(:, 2) > 200;
    if any(x_boundary_changes_right)
        particles(:, 2) = particles(:, 2) .* ~x_boundary_changes_right;
    end

    x_boundary_changes_left = particles(:, 2) < 0;
    if any(x_boundary_changes_left)
        particles(:, 2) = particles(:, 2) + 200 * x_boundary_changes_left - abs(particles(:, 2) .* x_boundary_changes_left);
    end
    
    y_boundary_changes_upper = particles(:, 1) > 100;
    if any(y_boundary_changes_upper)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* y_boundary_changes_upper);
        overshoot = (particles(:, 1) - 100) .* y_boundary_changes_upper;
        particles(:, 1) = particles(:, 1) - 2 * overshoot;
    end
    
    y_boundary_changes_lower = particles(:, 1) < 0;
    if any(y_boundary_changes_lower)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* y_boundary_changes_lower);
        overshoot = abs(particles(:, 1)) .* y_boundary_changes_lower;
        particles(:, 1) = particles(:, 1) + 2 * overshoot;
    end
    
    % question 2 scattering logic
    scattered = rand(num_particles, 1) < pscat;
    if any(scattered)
        scattered(:, 2) = scattered;
        free_path = abs(position_of_last_scatter - scattered(:, 2) .* particles(:, 3:4));
        mean_free_path = mean(free_path, "all");
        angles = randn(num_particles, 1) .* 2 * pi;
        rethermalized_velocities(:, 1:2) = randn(num_particles, 2) * sqrt(kb * T / m) / 1E15;
        rethermalized_velocities = rethermalized_velocities .* scattered;
        particles(:, 3:4) = particles(:, 3:4) .* ~scattered + rethermalized_velocities;
        position_of_last_scatter = particles(:, 3:4);
    end
    
    % plot position updates of particles, 
    % ignoring those that passed the x-boundary
    x_boundary_affected = x_boundary_changes_right | x_boundary_changes_left;
    temp_avg = mean(((sqrt(particles(:, 3).^2 + particles(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    % subplot(2, 1, 2)
    hold on
    figure(2)
    title(sprintf("Avg Temperature: %s, Mean Free Path is: %s", temp_avg, mean_free_path))
    for i = 1:length(particles)
        if ~x_boundary_affected(i)
            plot([previous_particles(i, 2) particles(i, 2)], [previous_particles(i, 1) particles(i, 1)], colours(mod(i, length(colours)) + 1))
        end
    end
    axis([0 200 0 100])
    pause(0.1)
end

