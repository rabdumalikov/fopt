function [error_msg, start_point, z, iteration_counter, gnorm, dx] = newton_multidimensional( axis, obj_func, start_point,  tol, maxiter, is_point_within_range )

%% minimum allowed perturbation
dxmin = 1e-6;
 
%% initialize gradient norm, optimization vector, iteration counter, perturbation
error_msg = ''; gnorm = inf; iteration_counter = 0; dx = inf;
 
z = obj_func( start_point(1), start_point(2) );

% gradient descent algorithm:
while gnorm >= tol && iteration_counter <= maxiter && dx >= dxmin
    
    %% plot start point
    plot3( axis, start_point(1), start_point(2), z, 'r*', 'LineWidth', 2 );
    
    %% computing next point  
    A = compute_hessian_and_evaluate( obj_func, start_point );
    
    [L, D] = ldl(A);
    
    D = replace_negative_entries_with_delta( D );

    new_mtrx_A = L * D * L';

    newton_step = -(new_mtrx_A \ compute_gradient_and_evaluate( obj_func, start_point ));
    
    % figure out alpha
    alpha = linesearch( obj_func, start_point, newton_step );
    
    next_point = start_point' + alpha * newton_step;        
        
    %% check steps
    if ~isfinite(next_point)
        display(['Number of iterations: ' num2str(iteration_counter)])
        error('x is inf or NaN')
    end
     
    if ~is_point_within_range(next_point)
        error_msg = sprintf( "Point(%d,%d) is Out Of Range! Solution: Decrease 'Alpha' or Increase 'Range'", next_point );
        return
    end
    
    next_z = obj_func( next_point(1),next_point(2) );
    
    %% plot next point
    plot3( axis, next_point(1), next_point(2), next_z, 'r*', 'LineWidth', 4 );
    %% plot line between old and new point
    plot3( axis, [start_point(1) next_point(1) ],[start_point(2) next_point(2)], [z next_z], 'ko-', 'LineWidth', 1);
    pause( 0.1 );
    hold( axis, 'on' );
    
    %% update termination metrics and general values
    iteration_counter = iteration_counter + 1;
    dx = norm(next_point-start_point);
    z = next_z;
    start_point = next_point';
    
end

%% plot last point 
plot3( axis, next_point(1),next_point(2), next_z,'ko', 'LineWidth', 10 );
pause( 1 );


%% replace negative entries of matrix with some delta
function result = replace_negative_entries_with_delta( D )
    
    delta = 0.1;
    
    [m,n] = size(D);
    for i=1:m
        if D(i,i) <= 0
            D(i,i) = delta;
        end
    end
      
    result = D;
    
%% calculate gradient of the objective function
function eval_gradient = compute_gradient_and_evaluate(obj_func, point)
    syms x y
    %% computing symbolic gradient
    gradient_vector = gradient(obj_func, [x,y]);
    
    %% evaluate computed gradient
    eval_gradient = [
        eval( subs(gradient_vector(1), [x y], point))
        eval( subs(gradient_vector(2), [x y], point))
    ];
    %% dump
    fprintf( 1, '  -OldPoint(%5.7f,%5.7f)- ===>', point );
    fprintf( 1, '  -NewPoint(%5.7f,%5.7f)-\n', eval_gradient );
    
    
%% calculate gradient of the objective function
function eval_hessian = compute_hessian_and_evaluate(obj_func, point)
    syms x y
    %% computing symbolic hessian
    hessian_mtrx = hessian(obj_func, [x,y]);
    
    %% evaluate computed gradient
    eval_hessian = [
        eval(subs(hessian_mtrx(1,1), [x y], point)), eval(subs(hessian_mtrx(1,2), [x y], point));
        eval(subs(hessian_mtrx(2,1), [x y], point)), eval(subs(hessian_mtrx(2,2), [x y], point))
    ];
    %% dump
    fprintf ( 1, '  -OldPoint(%5.7f,%5.7f)- ===>', point );
    fprintf ( 1, '  -NewPoint(%5.7f,%5.7f,%5.7f,%5.7f)-\n', eval_hessian );
    
%% hk - is search direction
function lambda = linesearch( obj_func, start_point, newton_step )
    lambda = 1;
 
    obj_func_for_vector = @(start_point) obj_func(start_point(1),start_point(2));
    
    while obj_func_for_vector( start_point' + lambda * newton_step ) >= obj_func_for_vector( start_point ) && lambda ~= 0
        lambda = lambda / 2;
    end
    
    if lambda == 0
        lambda = 0.001; % fallback
    end
 