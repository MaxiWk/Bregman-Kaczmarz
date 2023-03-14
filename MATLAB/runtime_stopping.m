% stopping after a specified computation time
classdef runtime_stopping

    properties
        max_time
        iter_save
        maxiter = 1e9;
    end

    methods

        function obj = runtime_stopping(input)
                 obj.max_time = input.max_time;
                 obj.iter_save = input.iter_save;
        end

    end

end