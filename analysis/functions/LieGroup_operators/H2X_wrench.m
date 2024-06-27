function [X_wrench] = H2X_wrench(H)
% from ^A H_B to _A X^B
R = H(1:3,1:3);

X_wrench = [R               , zeros(3,3);
            hat(H(1:3,4))*R,        R  ];
end