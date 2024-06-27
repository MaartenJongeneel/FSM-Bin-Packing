function q = PlotFrame(H,n,s, lineWidth)
% PlotFrame() Draws a 3D orthogonal reference frame placed and
% oriented as indicated in the homogeneous matrix wTa wrt the world frame
%
% Syntax:  PlotFrame(H,n,s)
%
% Inputs:
%    H               - Time series of homogeneous matrices of size 4x4xN
%    n               - label applied to the axes. The default is 'none'
%    s               - scale: determines the lenght of the frame axes. The 
%                             scale is 0.001 [m] by default
%    lineWidth       - Tickness of the arrows
%
% Outputs:
%    q               - Handle to the quiver3 object
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 15-September-2023
%--------------------------------------------------------------------------

    if nargin == 1
        n = 'none';
        s = 0.001;
        lineWidth = 3;
    elseif nargin == 2
        s = 0.001;
        lineWidth = 3;
    elseif nargin == 3
        lineWidth = 3;
    end

    pos = repmat(H(1:3,4)',3,1);

    q = quiver3(pos(1:3,1),pos(1:3,2),pos(1:3,3),H(1,1:3)',H(2,1:3)',H(3,1:3)');
    q.AutoScaleFactor = s;

    % shift the label of Y axis a bit in order to avoid labels overlap
    H(1:3,2)=H(1:3,2)+0.1*ones(3,1);
    if isnumeric(n)
        n=num2str(n);
        p=strcat('{',n,'}$');
        labels = [strcat('$X_',p);strcat('$Y_',p);strcat('$Z_',p)];

    else
        if strcmp(n,'none')
            labels = ['','',''];
        else
            p=strcat('{',n,'}$');
            labels = [strcat('$X_',p);strcat('$Y_',p);strcat('$Z_',p)];
        end
    end


    text(pos(1:3,1)+s*H(1,1:3)',pos(1:3,2)+s*H(2,1:3)',pos(1:3,3)+s*H(3,1:3)',labels);

    cd_tail = uint8([blkdiag([255 255],[255 255],[255 255]);255*ones(1,6)]);
%     cd_tail = zeros(4,6); cd_tail(3,:) = 255; cd_tail(4,:) = 255;
    cd_tail = uint8(cd_tail);
    set(q.Tail,'ColorBinding','interpolated','ColorData',cd_tail,'LineWidth',lineWidth);
    cd_head = uint8([blkdiag([255 255 255],[255 255 255],[255 255 255]);255*ones(1,9)]);
%     cd_head = zeros(4,9); cd_head(3,:) = 255; cd_head(4,:) = 255;
    cd_head = uint8(cd_head);
    set(q.Head,'ColorBinding','interpolated','ColorData',cd_head,'LineWidth',lineWidth);

end
