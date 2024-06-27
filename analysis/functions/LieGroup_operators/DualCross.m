function C = DualCross(v)
C = [hat(v(4:6)), zeros(3);...
     hat(v(1:3)), hat(v(4:6))];    
end