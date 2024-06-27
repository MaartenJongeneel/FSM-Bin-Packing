function w = vee(W)
%vee Provides the vector representation of a input so(3) or se(3) matrix 
%
% Syntax:  w = vee(R)
%
% Inputs:
%    W        - Skew-symmetric matrix in so(3)
%             - Twist matrix in se(3)
%
% Output:
%    w        - vector representation of R in so(3) -> w in R^3
%             - vector representation of R in se(3) -> w in R^6
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 18-July-2023
%--------------------------------------------------------------------------

    w = [];
    if ( size(W,1) == 3  && size(W,2) == 3 )
        w = zeros(3,1);
        R0 = W+W';
        for idx=1:3
           for idy=1:3
              if abs( R0(idx,idy) ) < 5e-4
                  R0(idx,idy) = 0;
              end
           end
        end
        if all( ( all(R0) == 0 ) )
            w(1,1) = real( W(3,2) );
            w(2,1) = real( W(1,3) );
            w(3,1) = real( W(2,1) );
        else
            if isnan(W(1,1))
                w = [NaN;NaN;NaN];
            else
                fprintf('Input Matrix W is not Skew-Symmetric \n');
                W
            end
        end
    elseif (  size(W,1) == 4  && size(W,2) == 4 )
        w(1:3) = W(1:3,4);
        w(4:6) = vee(W(1:3,1:3));
    end

end
