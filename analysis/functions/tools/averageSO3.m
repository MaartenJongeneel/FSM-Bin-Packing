function R_ = averageSO3(R)
%averageSO3 Provides the mean rotation matrix
% given a time series of noisy rotation matices as input, the function
% remove outliers and computes the average rotation matrix in SO(3)
%
% Syntax:  R_ = averageSO3(R)
%
% Inputs:
%    R        - Time series of noisy rotation matrices R of size 3x3xN
%
% Output:
%    R_       - Average rotation matrix
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 12-July-2023
%--------------------------------------------------------------------------
    
    th = 1e-6;
    RotVec = zeros(size(R,3) ,3);
    R_guess = R(:,:,1);

    while true
        % project all the rotation matrices to the tangent space of R_guess
        for ii=1:size(R,3)   
            RotVec(ii,:) = vee( logm( R(:,:,ii) * R_guess' ) )';
        end 
        % remove outliers
        RotVec_ = rmoutliers( RotVec );
        % average the vectors on the tangent space of R_guess
        phi_ = mean( RotVec_ );
        % bring back the averaged vector on the tangent space to the
        % manifold SO(3). This is the approximation of the averaged
        % rotation matrix.
        R_guess = expm( hat( phi_ ) ) * R_guess;
        % If the norm of the averaged rotation vector is close enough to
        % zero, then R_guess is a good approximation of the average
        % rotation matrix
        if norm(phi_) <= th
            break
        end
        
    end
    
    R_ = R_guess;
end

