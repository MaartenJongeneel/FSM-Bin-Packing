function fig = Plot6DtimeTrajectory(H, t, n, figTitle,c)
%Plot6DtimeTrajectory - Plots the six components
%
% Syntax:  fig = Plot6DtimeTrajectory(H, t, n, figTitle, verbose)
%
% Inputs:
%    H               - Time series of homogeneous matrices of size 4x4xN
%    t               - Time vector associated to each homegeneous matrix
%    n               - Figure number: n = 1 by default 
%    figTitle        - Title of the figure
%    c               - Plot color
%
% Outputs:
%    fig             - Handle of the figure window
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 13-July-2023
%--------------------------------------------------------------------------

    if nargin == 1
        n = 1;
        t = 0:1:size(H,3)-1;
        c = 'b.';
        figTitle = '';
    elseif nargin == 2 
        n = 1;
        c = 'b.';
        figTitle = '';
    elseif nargin == 3
        c = 'b.';
        figTitle = '';
    elseif nargin == 4
        c = 'b.';
    end
    
    ax = ['x' ; 'y'; 'z' ];
    
    fig = figure(n);
    sgtitle(figTitle)
    for jj=1:size(H,3)   
        RotVec(jj,:) = vee( logm( H(1:3,1:3,jj) ) )';
    end
    for ii=0:2
        subplot(3,2,ii*2+1)
        plot(t,squeeze(H(ii+1,4,:)),c); 
        hold on
        grid on
        xlabel('time [s]')
        ylabel(['$p_' ax(ii+1) '$ [m]'])
        
        subplot(3,2,(ii+1)*2)
        plot(t,(RotVec(:,ii+1)),c)
        hold on
        grid on
        xlabel('time [s]')
        ylabel(['$\theta_' ax(ii+1) '$ [rad]'])

    end



end
