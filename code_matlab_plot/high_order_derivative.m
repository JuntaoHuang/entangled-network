
function df = high_order_derivative(f, x)
    % Input:     
    % f - vector of function values at x points
    % x - vector of x points
    
    % Output:
    % df - derivative values at x points using fourth-order finite differences
    
    h = x(2) - x(1); % Assuming uniform spacing

    % Pre-allocate df array
    df = zeros(size(f));

    % Apply fourth-order centered difference for interior points
    for i = 3:length(x)-2
        df(i) = (-f(i+2) + 8*f(i+1) - 8*f(i-1) + f(i-2)) / (12*h);
    end

    % Apply second-order forward/backward difference at boundaries
    df(1) = (-3*f(1) + 4*f(2) - f(3)) / (2*h);
    df(2) = (-3*f(2) + 4*f(3) - f(4)) / (2*h);
    df(end-1) = (3*f(end-1) - 4*f(end-2) + f(end-3)) / (2*h);
    df(end) = (3*f(end) - 4*f(end-1) + f(end-2)) / (2*h);
end
