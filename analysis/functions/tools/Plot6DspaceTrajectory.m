function fig = Plot6DspaceTrajectory(H,n,s,d,label,figTitle,lineWidth)
%Plot6DspaceTrajectory() 3D plot of the frames along the path defined by H
%   The 3D path, both in position and orientation, given by the sequence of
%   homogeneous matrices H is drawn in a plot3
%
% Syntax:  fig = Plot6DspaceTrajectory(H,n,s,d,label,figTitle,lineWidth)
%
% Inputs:
%    H               - Time series of homogeneous matrices of size 4x4xN
%    n               - Figure number: n = 1 by default 
%    s               - scale: determines the lenght of the frame axes. The 
%                             scale is 0.001 [m] by default
%    d               - Number of steps between frames to draw. Default n=1
%    label           - Label applied to the axes. The default is 'none'
%    figTitle        - Title of the figure
%    lineWidth       - Tickness of the arrows 
%
% Outputs:
%    fig             - Handle of the figure window
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 15-September-2023
%--------------------------------------------------------------------------

    if nargin == 1
        n = 1;
        s = 0.001;
        d = 1;
        label = 'none';
        figTitle = '';
        lineWidth = 3;
    elseif nargin == 2
        s = 0.001;
        d = 1;
        label = 'none';
        figTitle = '';
        lineWidth = 3;
    elseif nargin == 3
        d = 1;
        label = 'none';
        figTitle = '';
        lineWidth = 3;
    elseif nargin == 4
        label = 'none';
        figTitle = '';
        lineWidth = 3;
    elseif nargin == 5
        figTitle = '';
        lineWidth = 3;
    elseif nargin == 6  
        lineWidth = 3;
    end
    fig = figure(n);
    title(figTitle, 'interpreter', 'latex')
    
    hold on 
    for ii=1:d:size(H,3)
       PlotFrame(H(:,:,ii),label, s, lineWidth); 
%        PlotFrame([H(1:3,1:3,ii), [H(1:2,4,ii); H(3,4,ii)*5]; 0 0 0 1],label, s, lineWidth); 
    end
    hold off
    grid on
    axis equal



end

