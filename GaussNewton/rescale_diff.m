function rescale_value = rescale_diff( norm_value, min_range_value, max_range_value )

% Rescale normalized values on an array between -1 and 1 back 
rescale_value = 0.5*norm_value .* (max_range_value - min_range_value);

return