% stopping after a specified residual value 
classdef res_stopping

    properties
        iter_save
        max_res
        maxiter = 1e9;
    end

    methods

        function obj = res_stopping(input)
                 obj.iter_save = input.iter_save;
                 obj.max_res = input.max_res;
        end

    end

end