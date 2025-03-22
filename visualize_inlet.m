function visualize_inlet()
    % Run the calculation first to get the angles
    [betas, deflection_angles, total_deflection_angles, ~, ~, ~] = ramp_angle_calc();
    
    % Set up figure
    figure('Position', [0, 0, 1200, 700]);
    hold on;
    grid on;
    axis equal;

    % Number of ramps
    n = length(deflection_angles);
    
    % Initialize arrays for ramp points and shock angles
    ramp_points = zeros(2, n+1);
    shock_angles = zeros(1, n+1);
    
    % First ramp point is origin
    ramp_points(:, 1) = [0; 0];
    
    % Calculate shock angles for each section
    for i = 1:n
        shock_angles(i) = betas(i) + total_deflection_angles(i) - deflection_angles(i);
    end
    % Final normal shock angle
    shock_angles(n+1) = total_deflection_angles(n) + 90;
    
    % Find intersection point of first two shocks
    first_ramp_end = [cosd(total_deflection_angles(1)); sind(total_deflection_angles(1))];
    intersection = find_intersection(ramp_points(:, 1), shock_angles(1), first_ramp_end, shock_angles(2));
    
    % Calculate ramp points
    for i = 2:n+1
        if i == 2
            % First ramp point after origin
            ramp_points(:, i) = [cosd(total_deflection_angles(i-1)); sind(total_deflection_angles(i-1))];
        else
            % Find intersection of previous ramp angle and next shock
            ramp_points(:, i) = find_intersection(ramp_points(:, i-1), total_deflection_angles(i-1), intersection, shock_angles(i));
        end
    end
    
    % Plot ramps
    plot(ramp_points(1, :), ramp_points(2, :), "LineWidth", 2, "Color", "k");
    
    % Plot shocks and add beta angle labels
    for i = 1:n+1
        shock_line = [ramp_points(:, i), intersection];
        plot(shock_line(1, :), shock_line(2, :), "LineWidth", 2, "LineStyle", ":", "Color", "r");
        
        % Add beta angle labels for shocks (excluding the final normal shock)
        if i <= n
            % Find midpoint of shock line for label placement
            mid_point = 0.6 * ramp_points(:, i) + 0.4 * intersection;
            
            % Ensure label is on the right side 
            % Calculate vector perpendicular to shock line
            shock_vector = intersection - ramp_points(:, i);
            perp_vector = [shock_vector(2); -shock_vector(1)];
            perp_vector = perp_vector / norm(perp_vector) * 0.03;
            
            % Check if perpendicular vector points to the right, if not, reverse it
            if perp_vector(1) < 0
                perp_vector = -perp_vector;
            end
            
            % Add the beta angle label
            text(mid_point(1) + perp_vector(1), mid_point(2) + perp_vector(2), ...
                 ['\beta_' num2str(i) ' = ' num2str(betas(i), '%.1f') '°'], ...
                 'FontSize', 9, 'HorizontalAlignment', 'left');
        else
            % Label for normal shock
            mid_point = 0.6 * ramp_points(:, i) + 0.4 * intersection;
            text(mid_point(1) + 0.03, mid_point(2), 'Normal Shock', ...
                 'FontSize', 9, 'HorizontalAlignment', 'left');
        end
    end
    
    % Add deflection angle labels for ramps
    for i = 1:n
        if i == 1
            % First ramp
            mid_point = 0.5 * (ramp_points(:, i) + ramp_points(:, i+1));
        else
            % Subsequent ramps - place label closer to the ramp junction
            mid_point = 0.7 * ramp_points(:, i) + 0.3 * ramp_points(:, i+1);
        end
        
        % Calculate ramp vector
        ramp_vector = ramp_points(:, i+1) - ramp_points(:, i);
        
        % Calculate perpendicular vector (pointing to the right)
        perp_vector = [ramp_vector(2); -ramp_vector(1)];
        perp_vector = perp_vector / norm(perp_vector) * 0.03;
        
        % Check if perpendicular vector points to the right, if not, reverse it
        if perp_vector(1) < 0
            perp_vector = -perp_vector;
        end
        
        % Add the deflection angle label
        text(mid_point(1) + perp_vector(1), mid_point(2) + perp_vector(2), ...
             ['\delta_' num2str(i) ' = ' num2str(deflection_angles(i), '%.1f') '°'], ...
             'FontSize', 9, 'HorizontalAlignment', 'left');
    end
    
    % Add labels
    title('Inlet Geometry with Shock Waves');
    xlabel('x');
    ylabel('y');
    legend('Inlet Ramps', 'Shock Waves', 'Location', 'northwest');
    
    % Get current axis limits
    ax = gca;
    x_limits = get(ax, 'XLim');
    y_limits = get(ax, 'YLim');
    
    % Expand limits by 0.15 in each direction for more room around labels
    new_x_limits = [x_limits(1)-0.15, x_limits(2)+0.15];
    new_y_limits = [y_limits(1)-0.15, y_limits(2)+0.15];
    
    % Apply new limits
    set(ax, 'XLim', new_x_limits);
    set(ax, 'YLim', new_y_limits);
    
    % Set fixed tick spacing of 0.2 for both axes
    xticks(floor(new_x_limits(1)/0.2)*0.2:0.2:ceil(new_x_limits(2)/0.2)*0.2);
    yticks(floor(new_y_limits(1)/0.2)*0.2:0.2:ceil(new_y_limits(2)/0.2)*0.2);
    
    % Ensure gridlines match ticks
    grid on;
end

function intersection = find_intersection(p1, angle1, p2, angle2)
    % Calculate slopes
    m1 = tand(angle1);
    m2 = tand(angle2);

    % Calculate y-intercepts
    b1 = p1(2) - m1 * p1(1);
    b2 = p2(2) - m2 * p2(1);

    % Find intersection
    x_int = (b2 - b1) / (m1 - m2);
    y_int = m1 * x_int + b1;

    intersection = [x_int; y_int];
end