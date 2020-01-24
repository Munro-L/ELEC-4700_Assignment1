num_particles = 7;
colours = ["r", "g", "b", "c", "y", "k", "m"];
particles = rand(num_particles, 4);
particles(:, 2) = particles(:, 2)*200;  % x-coordinates
particles(:, 1) = particles(:, 1)*100;  % y-coordinates
particles(:, 3) = particles(:, 3)*10;   % placeholder x-velocity
particles(:, 4) = particles(:, 4)*10;   % placeholder y-velocity

%scatter(particles(:, 2), particles(:, 1))
%axis([0 200 0 100])

for i = 0:100
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
    
    % refresh the plot
    %scatter(particles(:, 2), particles(:, 1))
    for i = 1:length(particles)
        plot([previous_particles(i, 2) particles(i, 2)], [previous_particles(i, 1) particles(i, 1)], colours(i))
    end
    axis([0 200 0 100])
    hold on
    pause(1)
    %refreshdata
    %drawnow
end

