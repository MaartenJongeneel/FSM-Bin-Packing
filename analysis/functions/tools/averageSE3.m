function H_ = averageSE3(H,method)
%averageSE3 Provides the mean homogeneous matrix
% given a time series of noisy homogeneous matices as input, the function
% remove outliers and computes the average homogeneous matrix in SE(3)
%
% Syntax:  H_ = averageSE3(H)
%
% Inputs:
%    H        - Time series of noisy homogeneous matrices H of size 4x4xN
%    method   - Averaging method: mean (default) or median
% Output:
%    H_       - Mean homogeneous matrix
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 12-July-2023
%--------------------------------------------------------------------------

    o_ = zeros(3,1);

    if nargin < 2
        m = 'mean';
    elseif nargin == 2 
       if strcmp( method,'mean')
           m = 'mean';
       elseif strcmp( method,'median')
           m = 'median';
       else
           fprintf('[warning] averageSE3: wrong method chosen, mean is used as default \n');
           m = 'mean';
       end
    end
    
    % remove outliers to the position and average it
    if strcmp( m,'mean')
        o_ = mean( rmoutliers( squeeze( H(1:3,4,:) )' ) )';
    elseif strcmp( m,'median')
        o_ = median( rmoutliers( squeeze( H(1:3,4,:) )' ) )';
    end
    
    % compute the mean rotation matrix
    R_ = averageSO3( H(1:3,1:3,:) ) ;
    
    H_ = [  R_ ,o_; 
          0 0 0 1];
end

