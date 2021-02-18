%  Copyright (c) 2017-present Rustam Abdumalikov
%
%  "FunctionOptimizer" application
%
% Distributed under the Boost Software License, Version 1.0. (See
% accompanying file LICENSE_1_0.txt or copy at
% http://www.boost.org/LICENSE_1_0.txt)

function traces = newton_multidimensional( objective_func, start_point, maxiter, tol, step_size, activate_logs )

%% minimum allowed perturbation
dxmin = 1e-6;
 
%% gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; iteration_counter = 0; dx = inf;
 
z = objective_func( start_point(1), start_point(2) );

traces = [ start_point, z ];

% gradient descent algorithm:
while gnorm >= tol && iteration_counter <= maxiter && dx >= dxmin
    
    %% computing next point  
    A = compute_hessian_and_evaluate( objective_func, start_point, activate_logs );
    
    [L, D] = ldl(A);
    
    D = replace_negative_entries_with_delta( D );

    new_mtrx_A = L * D * L';

    newton_step = -1 * (new_mtrx_A \ compute_gradient_and_evaluate( objective_func, start_point, activate_logs ));
    
    % figure out 'alpha'
    alpha = linesearch( objective_func, start_point, newton_step, step_size );
    
    next_point = start_point' + alpha * newton_step;                
    
    %% update termination metrics and general values    
    next_z = objective_func(next_point(1), next_point(2));
     
    iteration_counter = iteration_counter + 1;
    dx = norm( next_point - start_point );
    z = next_z;
    start_point = next_point';
    
    %% save point    
    traces = [ traces; [ start_point, z ] ];
    
end

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
function eval_gradient = compute_gradient_and_evaluate(obj_func, point, activate_logs )
    syms x y
    %% computing symbolic gradient
    gradient_vector = gradient(obj_func, [x,y]);
    
    %% evaluate computed gradient
    eval_gradient = [
        eval( subs(gradient_vector(1), [x y], point))
        eval( subs(gradient_vector(2), [x y], point))
    ];
    %% dump
    if activate_logs
        fprintf( 1, '  -OldPoint(%5.7f,%5.7f)- ===>', point );
        fprintf( 1, '  -NewPoint(%5.7f,%5.7f)-\n', eval_gradient );
    end
    
    
%% calculate gradient of the objective function
function eval_hessian = compute_hessian_and_evaluate(obj_func, point, activate_logs)
    syms x y
    %% computing symbolic hessian
    hessian_mtrx = hessian(obj_func, [x,y]);
    
    %% evaluate computed gradient
    eval_hessian = [
        eval(subs(hessian_mtrx(1,1), [x y], point)), eval(subs(hessian_mtrx(1,2), [x y], point));
        eval(subs(hessian_mtrx(2,1), [x y], point)), eval(subs(hessian_mtrx(2,2), [x y], point))
    ];
    %% dump
    if activate_logs 
        fprintf ( 1, '  -OldPoint(%5.7f,%5.7f)- ===>', point );
        fprintf ( 1, '  -NewPoint(%5.7f,%5.7f,%5.7f,%5.7f)-\n', eval_hessian );
    end
    
%% hk - is search direction
function lambda = linesearch( obj_func, start_point, newton_step, step_size )
    lambda = step_size; % setting it to '1' doesn't work well.
 
    obj_func_for_vector = @(start_point) obj_func(start_point(1),start_point(2));
    
    while obj_func_for_vector( start_point' + lambda * newton_step ) >= obj_func_for_vector( start_point ) && lambda ~= 0
        lambda = lambda / 2;
    end
    
    if lambda == 0
        lambda = 0.001; % fallback
    end
 