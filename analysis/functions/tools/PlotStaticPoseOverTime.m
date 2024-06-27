function fig = PlotStaticPoseOverTime(H, t, n, figTitle, verbose)
%PlotStaticPoseOverTime - Plots the six components
%
% Syntax:  fig = PlotStaticPoseOverTime(H, t, n, figTitle, verbose)
%
% Inputs:
%    H               - Time series of homogeneous matrices of size 4x4xN
%    t               - Time vector associated to each homegeneous matrix
%    n               - Figure number: n = 1 by default 
%    figTitle        - Title of the figure
%    verbose         - Plot further information as:
%                         + data without outliers
%                         + mean and median
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
        fn = 1;
        time = 0:1:size(H,3)-1;
        full = 0;
    elseif nargin == 2 
        fn = 1;
        time = t;
        full = 0;
    elseif nargin == 3
        fn = n;
        time = t;
        full = 0;
    elseif nargin == 4
        fn = n;
        time = t;
        title = figTitle;
        full = 0;
    elseif nargin == 5
        fn = n;
        time = t;
        title = figTitle;
        full = verbose;
    end
    
    ax = ['x' ; 'y'; 'z' ];
    
    fig = figure(fn);
    sgtitle(title)
    if full
       H_ = averageSE3(H);
       H__= averageSE3(H,'median'); 
       phi_  = vee( logm( H_(1:3,1:3)  ) )';
       phi__ = vee( logm( H__(1:3,1:3) ) )';
    end
    for jj=1:size(H,3)   
        RotVec(jj,:) = vee( logm( H(1:3,1:3,jj) ) )';
    end
    for ii=0:2
        subplot(3,2,ii*2+1)
        if full
            plot(time,squeeze(H(ii+1,4,:)), 'm.'); 
            hold on
            [H_o, TFrm ] = rmoutliers(squeeze(H(ii+1,4,:)));
            plot(time(~TFrm),H_o, 'b.');
            plot([time(1) time(end)], H_(ii+1,4)*[1 1] , 'g-','LineWidth', 2);
            plot([time(1) time(end)], H__(ii+1,4)*[1 1], 'r-','LineWidth', 2);
            hold off
        else
            plot(time,squeeze(H(ii+1,4,:)), 'b.'); 
        end
        grid on
        xlabel('time [s]')
        ylabel(['p' ax(ii+1) ' [m]'])
        
        subplot(3,2,(ii+1)*2)
        if full 
            plot(time,(RotVec(:,ii+1)),'m.')
            hold on
            [RotVec_, TFrm ] = rmoutliers(squeeze(RotVec(:,ii+1)));
            plot(time(~TFrm),RotVec_,'b.') 
            plot([time(1) time(end)],phi_(1,ii+1)  *[1 1], 'g-','LineWidth', 2);
            plot([time(1) time(end)],phi__(1,ii+1) *[1 1], 'r-','LineWidth', 2);
            hold off
            grid on
        else
            plot(time,(RotVec(:,ii+1)),'b.')
        end
        xlabel('time [s]')
        ylabel(['theta' ax(ii+1) ' [rad]'])

    end
    if full
       legend('outliers', 'data w/o outliers', 'mean', 'median') 
    end


end

