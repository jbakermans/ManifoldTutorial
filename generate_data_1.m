function X_D = generate_data_1()
    % Make syntethic data in 2d
    t = linspace(0, 2*pi, 200);
    x = -(721*sin(t))/4 + 196/3*sin(2*t) - 86/3*sin(3*t) - 131/2*sin(4*t) + 477/14*sin(5*t) + 27*sin(6*t) - 29/2*sin(7*t) + 68/5*sin(8*t) + 1/10*sin(9*t) + 23/4*sin(10*t) - 19/2*sin(12*t) - 85/21*sin(13*t) + 2/3*sin(14*t) + 27/5*sin(15*t) + 7/4*sin(16*t) + 17/9*sin(17*t) - 4*sin(18*t) - 1/2*sin(19*t) + 1/6*sin(20*t) + 6/7*sin(21*t) - 1/8*sin(22*t) + 1/3*sin(23*t) + 3/2*sin(24*t) + 13/5*sin(25*t) + sin(26*t) - 2*sin(27*t) + 3/5*sin(28*t) - 1/5*sin(29*t) + 1/5*sin(30*t) + (2337*cos(t))/8 - 43/5*cos(2*t) + 322/5*cos(3*t) - 117/5*cos(4*t) - 26/5*cos(5*t) - 23/3*cos(6*t) + 143/4*cos(7*t) - 11/4*cos(8*t) - 31/3*cos(9*t) - 13/4*cos(10*t) - 9/2*cos(11*t) + 41/20*cos(12*t) + 8*cos(13*t) + 2/3*cos(14*t) + 6*cos(15*t) + 17/4*cos(16*t) - 3/2*cos(17*t) - 29/10*cos(18*t) + 11/6*cos(19*t) + 12/5*cos(20*t) + 3/2*cos(21*t) + 11/12*cos(22*t) - 4/5*cos(23*t) + cos(24*t) + 17/8*cos(25*t) - 7/2*cos(26*t) - 5/6*cos(27*t) - 11/10*cos(28*t) + 1/2*cos(29*t) - 1/5*cos(30*t);
    y = -(637*sin(t))/2 - 188/5*sin(2*t) - 11/7*sin(3*t) - 12/5*sin(4*t) + 11/3*sin(5*t) - 37/4*sin(6*t) + 8/3*sin(7*t) + 65/6*sin(8*t) - 32/5*sin(9*t) - 41/4*sin(10*t) - 38/3*sin(11*t) - 47/8*sin(12*t) + 5/4*sin(13*t) - 41/7*sin(14*t) - 7/3*sin(15*t) - 13/7*sin(16*t) + 17/4*sin(17*t) - 9/4*sin(18*t) + 8/9*sin(19*t) + 3/5*sin(20*t) - 2/5*sin(21*t) + 4/3*sin(22*t) + 1/3*sin(23*t) + 3/5*sin(24*t) - 3/5*sin(25*t) + 6/5*sin(26*t) - 1/5*sin(27*t) + 10/9*sin(28*t) + 1/3*sin(29*t) - 3/4*sin(30*t) - (125*cos(t))/2 - 521/9*cos(2*t) - 359/3*cos(3*t) + 47/3*cos(4*t) - 33/2*cos(5*t) - 5/4*cos(6*t) + 31/8*cos(7*t) + 9/10*cos(8*t) - 119/4*cos(9*t) - 17/2*cos(10*t) + 22/3*cos(11*t) + 15/4*cos(12*t) - 5/2*cos(13*t) + 19/6*cos(14*t) + 7/4*cos(15*t) + 31/4*cos(16*t) - cos(17*t) + 11/10*cos(18*t) - 2/3*cos(19*t) + 13/3*cos(20*t) - 5/4*cos(21*t) + 2/3*cos(22*t) + 1/4*cos(23*t) + 5/6*cos(24*t) + 3/4*cos(26*t) - 1/2*cos(27*t) - 1/10*cos(28*t) - 1/3*cos(29*t) - 1/19*cos(30*t);
    % Scale between 0 and 1, then flip (for better end result)
    x = (x - min(x)) / (max(x) - min(x));
    y = (y - min(y)) / (max(y) - min(y));
    % Define rotation functions
    rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
    roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
    rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
    % Do some rotations to create a 3d trajectory
    X_D = rotz(deg2rad(30))*roty(deg2rad(-30))*rotx(deg2rad(15))*...
        [x; y; 0.1*rand(size(x))];
end