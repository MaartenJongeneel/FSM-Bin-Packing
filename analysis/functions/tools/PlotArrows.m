function q = PlotArrows(p,n,s,c)
% PlotArrows() Draws a 3D orthogonal reference frame placed and
% oriented as indicated in the homogeneous matrix wTa wrt the world frame
%
% Syntax:  PlotArrows(p,n,s)
%
% Inputs:
%    p               - Matrix of all positions Nx3
%    n               - Matrix of vector normals Nx3
%    s               - scale: determines the lenght of the frame axes. The 
%                             scale is 0.001 [m] by default
%    c               - color: 'r', 'g', 'b'
%
% Outputs:
%    q               - Handle to the quiver3 object
% 
% Author: Alexander A. Oliva, Ph.D., Postdoctoral researcher
% Eindhoven University of Technology (TU/e), Mechanical Engineering Dept.
% email address: a.a.oliva@tue.nl  
% July 2023; Last revision: 13-July-2023
%--------------------------------------------------------------------------

    if nargin == 3
        asf = s;
        c = 'b';
    else
        asf = 0.001;
    end
    
    l = 1;
    
    if size(p,1) == 3
        pos = p';
        l = size(p,2);
    elseif size(p,2) == 3
        pos = p;
        l = size(p,1);
    else
        display('vector does not have the right dimension');
        return
    end
    
    if size(n,1) == 3
        vel = n';
    elseif size(n,2) == 3
        vel = n;
    else
        display('vector does not have the right dimension');
        return
    end


    q = quiver3(pos(:,1),pos(:,2),pos(:,3),vel(:,1),vel(:,2),vel(:,3));
    q.AutoScaleFactor = asf;
    
    if c == 'r'
        cd_tail = uint8([255*ones(1,2); zeros(2,2); 255*ones(1,2)]);
        cd_head = uint8([255*ones(1,3); zeros(2,3); 255*ones(1,3)]);
    elseif c == 'g'
        cd_tail = uint8([zeros(1,2); 255*ones(1,2);zeros(1,2); 255*ones(1,2)]);
        cd_head = uint8([zeros(1,3); 255*ones(1,3);zeros(1,3); 255*ones(1,3)]);
    elseif c == 'b'
        cd_tail = uint8([zeros(2,2); 255*ones(2,2)]);
        cd_head = uint8([zeros(2,3); 255*ones(2,3)]);
    end

%     cd_tail = uint8([zeros(2,2); 255*ones(2,2)]);
    cd_tail = repmat(cd_tail,1,l);
    set(q.Tail,'ColorBinding','interpolated','ColorData',cd_tail,'LineWidth',3);
%     cd_head = uint8([zeros(2,3); 255*ones(2,3)]);
    cd_head = repmat(cd_head,1,l);
    set(q.Head,'ColorBinding','interpolated','ColorData',cd_head,'LineWidth',3);
    
    axis equal
    hold on
    plot3(pos(:,1),pos(:,2),10*pos(:,3),'bo','MarkerFaceColor','b','MarkerSize', 5)

end
