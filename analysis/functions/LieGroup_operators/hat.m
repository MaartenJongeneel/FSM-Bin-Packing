function S = hat(w)
%hat Given the vector representation of an element of so(3) or se(3) matrix 
%
% Syntax:  R = hat(w)
%
% Inputs:
%    w        - vector representation w in R^3 of so(3)
%             - vector representation w in R^6 of se(3)
%
% Output:
%    S        - Skew-symmetric matrix S in so(3)
%             - Twist matrix S in se(3)
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 18-July-2023
%--------------------------------------------------------------------------

    % if the vector is given as a row, transpose it
    if (size(w,2) == 3) || (size(w,2) == 6)
        w = w';
    end
    
    S = [];

    if (size(w,1) == 3)
        
        S = zeros(3);
        S(2,1) =  w(3);      S(1,2) = -w(3);
        S(3,1) = -w(2);      S(1,3) =  w(2);
        S(3,2) =  w(1);      S(2,3) = -w(1);
  
    elseif size(w,1) == 6 
        S = [hat(w(4:6)), w(1:3); zeros(1,4)];
        
    else
        fprintf('hat(w): input vector w is of wrong size: %dx%d  \n', size(w,1),size(w,2));
    end

end
