% stopping after a specified number of iterations
classdef iteration_stopping

    properties
        maxiter
        iter_save
    end

    methods

        function obj = iteration_stopping(input)
                 obj.maxiter = input.maxiter;
                 obj.iter_save = input.iter_save;
        end

    end

end