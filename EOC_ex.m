function rate = EOC(e)
            % Determine rate of convergence for e_h,
            % e_2h , given row vector of errors where
            % each entry corresponds to a smaller h/2 value.
            n = length(e);
            rate = log( e((1:n-1))./(e(2:n)) )./log(2);
            
end