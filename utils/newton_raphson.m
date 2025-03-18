function [err, iterations, solution] = ...
                    newton_raphson(f,df,guess,tol,max_iter)
% Newton-Raphson Iterative Method in 1 dimension
%  INPUTS: function, function derivative, initial guess,
%          solution tolerance, maximum number of iterations
% OUTPUTS: absolute error in final guess, the number of total iterations
%          to converge on a solution, the solution

iterations = 1;

% until the maximum iterations is met
while (iterations <= max_iter)
    % find the absolute error in this guess
    err=f(guess)/df(guess);

    % come up with a new guess
    guess = guess - err;

    % if the guess is precise enough, this is the solution
    if abs(err) < tol
        solution = guess;
        return;
    % if not, iterate and try again
    else
        iterations = iterations + 1;
    end
end

% if maximum iterations is reached, return the last guess
solution = guess;