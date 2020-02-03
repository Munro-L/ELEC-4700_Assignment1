%% Question 3
% Generate constants to be used later in the script, identicle to Question
% 1 in this case. Staying with 10 oarticles for the sake of report clarity,
% but has been tested with up to 1000.
%
% Note: My weekend was unexpectedly busy, so I have a rough design for all
% component of Q3, however it does no work to the standard I would like. I plan
% to take the resubmission later for the 15% mark recovery.

num_particles = 10;
colours = ["r", "g", "b", "c", "y", "k", "m"];

kb = 1.38064852;
T = 300;
m = 9.10938356E-31;
vth = sqrt(kb * T / m) / 1E15;  % scaled to femtoseconds

mean_time_collision = 0.2;   % measured in picoseconds
timesteps = 1E-3;            % also in picoseconds, so this will be 1 femto second timesteps

%%
% Creates the positions of the particles and assigned a velocity to each,
% chosen from a normal distrobution about the thermal velocity. Question 3
% was not specific about whether to use the velocities of Question 1 or 2,
% so I am choosing to use the velocities of Question 1.

% particle positions
particles = rand(num_particles, 2);
particles(:, 2) = particles(:, 2)*200;  % x-coordinates
particles(:, 1) = particles(:, 1)*100;  % y-coordinates

% normal distribution shifted by themal velocity
% column 3 is x-velocity, 4 is x-velocity
angles = randn(num_particles, 1) .* 2 * pi;
particles(:, 3) = randn(num_particles, 1) + vth*cos(angles);
particles(:, 4) = randn(num_particles, 1) + vth*sin(angles);

%%
% The main time loop. 50 iterations, each step representing a
% femtosecond. This will use the same boundary conditions as Question 1 and
% 2, however with the addition of a bottleneck in the center. The
% bottleneck
% extends from x-coordinates 80 to 120, and y-coordinates 100 to 60 and 0
% to 40. The variable "specular" controls whether boundaries are to be
% treated as specular or diffusive.

% main time loop
specular = 1;    % set to 0 for diffusive
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

    %%
    % The following block is the logic for handling the bottleneck
    % boundaries
    
    % bottleneck boundaries
    bottom_left_changes = particles(:, 1) < 40 & particles(:, 2) > 80 & particles(:, 2) < 120 & particles(:, 3) > 0;
    if any(bottom_left_changes)
        particles(:, 3) = particles(:, 3) - (2 * particles(:, 3) .* bottom_left_changes);
    end
    
    bottom_right_changes = particles(:, 1) < 40 & particles(:, 2) > 80 & particles(:, 2) < 120 & particles(:, 3) < 0;
    if any(bottom_right_changes)
        particles(:, 3) = particles(:, 3) - (2 * particles(:, 3) .* bottom_right_changes);
    end
    
    bottom_top_changes = particles(:, 1) < 40 & particles(:, 2) > 80 & particles(:, 2) < 120 & particles(:, 4) < 0;
    if any(bottom_top_changes)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* bottom_top_changes);
    end
    
    top_left_changes = particles(:, 1) > 60 & particles(:, 2) > 80 & particles(:, 2) < 120 & particles(:, 3) > 0;
    if any(top_left_changes)
        particles(:, 3) = particles(:, 3) - (2 * particles(:, 3) .* top_left_changes);
    end
    
    top_right_changes = particles(:, 1) > 60 & particles(:, 2) > 80 & particles(:, 2) < 120 & particles(:, 3) < 0;
    if any(top_right_changes)
        particles(:, 3) = particles(:, 3) - (2 * particles(:, 3) .* top_right_changes);
    end
    
    top_bottom_changes = particles(:, 1) < 40 & particles(:, 2) > 80 & particles(:, 2) < 120 & particles(:, 4) < 0;
    if any(top_bottom_changes)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* top_bottom_changes);
    end
    
    
    %%
    % This section plots particle position updates
    
    % plot position updates of particles, 
    % ignoring those that passed the x-boundary
    x_boundary_affected = x_boundary_changes_right | x_boundary_changes_left;
    temp_avg = mean(((sqrt(particles(:, 3).^2 + particles(:, 4).^2) .* 1E15).^2) .* m ./ kb);
    title(sprintf("Avg Temperature: %s", temp_avg))
    for i = 1:length(particles)
        if ~x_boundary_affected(i)
            plot([previous_particles(i, 2) particles(i, 2)], [previous_particles(i, 1) particles(i, 1)], colours(mod(i, length(colours)) + 1))
        end
    end
    axis([0 200 0 100])
    hold on
    pause(0.1)
end

%%
% This section constructs a density map based on the final positions of all
% particles.

density_data = ceil((particles(:, 1:2)./10)).*10;      % fancy way of rounding to the nearest 10
grid = zeros(100, 200);
for i = 1:length(density_data)
    grid(density_data(i, 1), density_data(i, 2)) = grid(density_data(i, 1), density_data(i, 2)) + 1;
end
imagesc(grid)


%% 
% Temperature map can be done in a similar way as density map, except
% instead of considering particle positions, we look at the velocities and
% find the temperature of each particle to show in the density map with their position.

heat_data = (particles(:, 3:4) .* 1E15 .* m ./ kb).^2;
position_data = ceil((particles(:, 1:2)./10)).*10;
grid = zeros(100, 200);
for i = 1:length(position_data)
    grid(density_data(i, 1), density_data(i, 2)) = grid(density_data(i, 1), density_data(i, 2)) + sqrt(heat_data(i, 1)^2 + heat_data(i, 2)^2);
end
imagesc(grid)



