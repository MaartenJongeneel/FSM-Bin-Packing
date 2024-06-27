function [X_vel] = H2X_vel(H)
% from ^A H_B to ^A X_B

R = H(1:3,1:3);
o = H(1:3,4);

X_vel = [     R    ,     hat(o)*R;
         zeros(3,3),       R    ];
end