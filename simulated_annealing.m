
    
function intermediate = simulated_annealing( axis, obj_func, start_point, tol, maxiter, is_point_within_range, report )

%% constants
Tinit = 100;        % initial temperature
alpha = 0.8;
max_consec_rejections = 100;
max_success = 20;

% counters etc
iteration_counter = 0;
amount_successes = 0;
finished = 0;
consec = 0;
T = Tinit;
start_point_value = obj_func(start_point);
total_amount_of_iterations = 0;

idx = 1
while ~finished;
    iteration_counter = iteration_counter + 1; % just an iteration counter
    
    intermediate( idx, : ) = start_point;
    idx = idx + 1;
    
    next_point = generate_new_point_in_range(start_point, is_point_within_range);
    next_point_value = obj_func(next_point);
    
    delta_f = next_point_value - start_point_value;
    
    %% new solution is better, so accept it.
    if delta_f < 0
        start_point = next_point;
        start_point_value = next_point_value;
        amount_successes = amount_successes + 1;
        consec = 0;
    else
        %% accept new solution, even bad one, based on probability.
        if( rand <= accept_solution( delta_f, T ) ) 
            start_point = next_point;
            start_point_value = next_point_value;
            amount_successes = amount_successes + 1;
        else
            consec = consec+1;
        end
    end
    
    %% Stop / decrement T criteria
    if iteration_counter >= maxiter || amount_successes >= max_success;
        
        total_amount_of_iterations = total_amount_of_iterations + iteration_counter;
        
        if T < tol || consec >= max_consec_rejections;
            finished = 1;
            break;
        else
            % decrease T according to cooling schedule
            T = lower_temperature( T, alpha ); 

            iteration_counter = 1;  
            amount_successes = 1;
        end
    end
end

minimum = start_point;
fval = start_point_value;

if report
    fprintf(1, '\n  Initial temperature:     \t%g\n', Tinit);
    fprintf(1, '  Final temperature:       \t%g\n', T);
    fprintf(1, '  Consecutive rejections:  \t%i\n', consec);
    fprintf(1, '  Number of function calls:\t%i\n', total_amount_of_iterations);
    fprintf(1, '  Total final loss:        \t%g\n', fval);
end

function result = lower_temperature( T, alpha )
    result = alpha * T;

function point = generate_new_point_in_range( current_point, is_point_within_range )
    point = generate_new_point( current_point );
    while ~is_point_within_range( point ) 
        point = generate_new_point(current_point); 
    end
    
function point = generate_new_point( current_point )
    point = current_point + ( randperm(length(current_point)) == length(current_point) ) * randn / ( rand * 5 );

function result = accept_solution( delta_f, current_temperature )
    result = exp( -delta_f/current_temperature);