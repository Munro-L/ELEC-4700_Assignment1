particles = rand(1000, 4);
particles(:, 2) = particles(:, 2)*200;  % x-coordinates
particles(:, 1) = particles(:, 1)*100;  % y-coordinates
particles(:, 3) = particles(:, 3)*10;   % placeholder x-velocity
particles(:, 4) = particles(:, 4)*10;   % placeholder y-velocity

scatter(particles(:, 2), particles(:, 1))
xlim([0 200])
ylim([0 100])

for i = 0:100
    % update positions
    particles(:, 1) = particles(:, 1) + particles(:, 4);
    particles(:, 2) = particles(:, 2) + particles(:, 3);
    
    
    % check if any particles passed the boundary and deal with them
    x_boundary_changes = particles(:, 2) > 200;
    if any(x_boundary_changes)
        particles(:, 2) = particles(:, 2) .* ~x_boundary_changes;
    end
    
    y_boundary_changes = particles(:, 1) > 100;
    if any(y_boundary_changes)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* y_boundary_changes);
        particles(:, 1) = y_boundary_changes * 100;
    end
    
    y_boundary_changes_lower = particles(:, 1) < 0;
    if any(y_boundary_changes_lower)
        particles(:, 4) = particles(:, 4) - (2 * particles(:, 4) .* y_boundary_changes_lower);
        particles(:, 1) = particles(:, 1) .* ~y_boundary_changes_lower;
    end
    
    
    scatter(particles(:, 2), particles(:, 1))
    xlim([0 200])
    ylim([0 100])
    refreshdata
    drawnow
end

