function [error_msg, start_point, z, iteration_counter, gnorm, dx] = gradient_descent( axis, obj_func, start_point, alpha, tol, maxiter, is_point_within_range, report )

%% minimum allowed perturbation
dxmin = 1e-6;

%% initialize error_msg, gradient norm, optimization vector, iteration counter, perturbation
error_msg = ''; gnorm = inf; iteration_counter = 0; dx = inf;

z = obj_func(start_point(1),start_point(2));

% gradient descent algorithm:
while gnorm >= tol && iteration_counter <= maxiter && dx >= dxmin
        
    %% plot current start point
    plot3(axis, start_point(1), start_point(2), z, 'r*', 'LineWidth', 2 );
    
    %% calculate gradient:
    eval_gradient = compute_gradient_and_evaluate(obj_func, start_point, report);
    gnorm = norm(eval_gradient);
    
    %% take step:
    next_point = start_point' - ( alpha * eval_gradient );
    
    %% check steps
    if ~isfinite(next_point)
        display(['Number of iterations: ' num2str(iteration_counter)])
        error('x is inf or NaN')
    end
    
    if ~is_point_within_range(next_point)
        error_msg = sprintf( "Point(%d,%d) is Out Of Range! Solution: Decrease 'Alpha' or Increase 'Range'", next_point );
        return
    end
    
    next_z = obj_func(next_point(1),next_point(2));

    %% plot path from old to new point
    plot3( axis, [start_point(1) next_point(1) ],[start_point(2) next_point(2)], [z next_z], 'ko-', 'LineWidth', 1 );
    pause( 0.1 );
    hold( axis, 'on' );
    
    %% update termination metrics and general values
    iteration_counter = iteration_counter + 1;
    dx = norm( next_point - start_point );
    z = next_z;
    start_point = next_point';
end

%% calculate gradient of the objective function
function eval_gradient = compute_gradient_and_evaluate(obj_func, point, report)
    syms x y
    %% computing symbolic gradient
    gradient_vector = gradient(obj_func, [x,y]);
    
    %% evaluate computed gradient
    eval_gradient = [
        eval( subs(gradient_vector(1), [x y], point) )
        eval( subs(gradient_vector(2), [x y], point) )
    ];
    %% dump
    if report
        fprintf( 1, '  -OldPoint(%5.7f,%5.7f)- ===>', point );
        fprintf( 1, '  -NewPoint(%5.7f,%5.7f)-\n', eval_gradient );
    end
    
    
