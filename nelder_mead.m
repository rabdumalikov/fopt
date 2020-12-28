function opt_value = nelder_mead( axis, start_point, objective_func, maxiter )

compute_b = @(c, n) c /( n * sqrt(2) ) * ( sqrt(n+1) - 1 );
compute_a = @(b, c) b + c / sqrt(2);

iteration_counter = 0;

[ ~, dimension ] = size ( start_point );

% constants
c = 10;

alpha = 1
beta  = 2
gamma = 1/2
delta = 1/2

b = compute_b( c, dimension + 1 );
a = compute_a( b, c );

x0 = start_point;
      
simplex_vertices = [ x0; x0 + [ a b ]; x0 + [ b a ] ];
plotit( axis, simplex_vertices, objective_func, 'red' )

while iteration_counter <= maxiter
    
    simplex_vertices = sort_by_function_values( objective_func, simplex_vertices );
       
    idx_for_best = 1;
    idx_for_bad = dimension;
    idx_for_worst = dimension + 1;
    
    best_point = simplex_vertices(idx_for_best,:);
    bad_point = simplex_vertices(idx_for_bad,:);
    worst_point = simplex_vertices(idx_for_worst,:);
    
    centroid = compute_centroid( dimension, simplex_vertices(1:dimension,:) )
    
    plot3(axis, centroid(1), centroid(2), objective_func( centroid ), 'r*', 'LineWidth', 2 );

    % Transformation
    %% calc reflection
    reflection_point = centroid + alpha * ( centroid - worst_point );

    worst_value = objective_func( worst_point );
    bad_value = objective_func( bad_point );
    reflection_value = objective_func( reflection_point );
    best_value = objective_func( best_point );
    
    point_to_replace = [];
    colour = []
    
    if best_value <= reflection_value && reflection_value < bad_value
        point_to_replace = reflection_point;
        colour = 'yellow'
        
    %% calc expansion
    elseif reflection_value < best_value
        expansion_point = centroid + beta * ( reflection_point - centroid );
        expantion_value = objective_func( expansion_point );
        
        if expantion_value < reflection_value
            point_to_replace = expansion_point;
        else
            point_to_replace = reflection_point;
        end
        
        colour = 'green';
        
    %% outside contraction
    elseif bad_value <= reflection_value && reflection_value < worst_value
        ocontraction_point = centroid + gamma * ( reflection_point - centroid );
        ocontraction_value = objective_func( ocontraction_point );

        if ocontraction_value <= reflection_value
            point_to_replace = ocontraction_point;
            colour = 'black';
        else
            %% shrinking 
            simplex_vertices = shrink( idx_for_best, dimension, simplex_vertices, best_point, delta );
            plotit( axis, simplex_vertices, objective_func, 'red' );
        end
    %% inside contraction
    elseif reflection_value >= worst_value
        icontraction_point = centroid - gamma * ( reflection_point - centroid );
        icontraction_value = objective_func( icontraction_point );

        if icontraction_value < worst_value
           point_to_replace = icontraction_point;
           colour = 'cyan';
        else
            %% shrinking 
            simplex_vertices = shrink( idx_for_best, dimension, simplex_vertices, best_point, delta );
            plotit( axis, simplex_vertices, objective_func, 'red' );
        end
    end
    
    if point_to_replace
        simplex_vertices(idx_for_worst,:) = point_to_replace;
        plotit( axis, simplex_vertices, objective_func, colour )
    end
    
    iteration_counter = iteration_counter + 1;
end

arr = evaluate_points( objective_func, simplex_vertices );
opt_value = min( arr );


function result = shrink( idx_for_best, dimension, simplex_vertices, best_point, delta )
    %% shrinking 
    result = simplex_vertices;
    
    for i=idx_for_best+1:dimension+1        
        result(i,:) = best_point + delta * ( simplex_vertices(i,:) - best_point );
    end
    
function result = evaluate_points( obj_func, points )
    [~, dim] = size(points);
    result = zeros(1, dim+1);
    for i=1:dim+1
        result(i) = feval( obj_func, points(i,:) )
    end


function res = compute_centroid( dimension, arr )
    res = sum( arr )/dimension;

function points = sort_by_function_values( obj_func, points )
    [~, dim] = size(points);
        
    eval_pts = evaluate_points( obj_func, points );   
    [ ~, index ] = sort( eval_pts );
    points = points(index,:);   

function plotit( axis, points, func, colour )
    x1 = points(1,:);
    x2 = points(2,:);
    x3 = points(3,:);
    
    plot3( axis, [x1(1) x2(1) ],[x1(2) x2(2)], [func(x1) func(x2)], 'ko-', 'color', colour, 'LineWidth', 2 );    
    plot3( axis, [x1(1) x3(1) ],[x1(2) x3(2)], [func(x1) func(x3)], 'ko-', 'color', colour, 'LineWidth', 2 );
    plot3( axis, [x2(1) x3(1) ],[x2(2) x3(2)], [func(x2) func(x3)], 'ko-', 'color', colour, 'LineWidth', 2 );
    
    pause( 0.2 );
    hold( axis, 'on' );