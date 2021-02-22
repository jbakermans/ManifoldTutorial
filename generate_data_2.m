function [X_n, X_m] = generate_data_2()
    % Make syntethic data in 2d
    t = linspace(0, 2*pi, 200);    
    x = 1/40*sin(23/9 - 100*t) + 1/6*sin(139/46 - 97*t) + 1/6*sin(41/10 - 96*t) + 1/8*sin(13/3 - 95*t) + 1/6*sin(3 - 93*t) + 1/6*sin(45/13 - 92*t) + 1/6*sin(9/2 - 91*t) + 1/5*sin(9/7 - 90*t) + 2/7*sin(11/7 - 89*t) + 1/7*sin(20/9 - 88*t) + 1/4*sin(11/8 - 87*t) + 1/4*sin(7/13 - 86*t) + 1/17*sin(19/5 - 85*t) + 1/5*sin(17/4 - 84*t) + 1/4*sin(29/8 - 82*t) + 1/5*sin(29/10 - 81*t) + 1/4*sin(4/9 - 79*t) + 3/7*sin(5/8 - 78*t) + 1/11*sin(1/5 - 77*t) + 1/3*sin(63/16 - 74*t) + 4/9*sin(21/5 - 73*t) + 4/9*sin(25/6 - 72*t) + 1/4*sin(11/6 - 71*t) + 1/35*sin(46/13 - 70*t) + 1/4*sin(13/6 - 69*t) + 1/5*sin(1/2 - 66*t) + 3/7*sin(9/2 - 65*t) + 1/3*sin(14/3 - 64*t) + 6/11*sin(5/2 - 63*t) + 1/6*sin(1/3 - 62*t) + 2/7*sin(23/7 - 61*t) + 5/11*sin(2/7 - 59*t) + 1/3*sin(33/7 - 58*t) + 1/6*sin(19/6 - 57*t) + 5/9*sin(1/7 - 56*t) + 11/12*sin(25/6 - 53*t) + 5/6*sin(33/10 - 51*t) + 4/7*sin(16/7 - 50*t) + 19/13*sin(13/5 - 48*t) + 3/7*sin(61/30 - 47*t) + 4/9*sin(24/7 - 46*t) + 11/10*sin(1/7 - 45*t) + 1/3*sin(1/11 - 44*t) + 12/7*sin(29/7 - 43*t) + 2/7*sin(35/8 - 42*t) + 3/8*sin(11/10 - 41*t) + 3/4*sin(7/3 - 40*t) + 31/16*sin(8/9 - 39*t) + 12/7*sin(1/10 - 38*t) + 4/5*sin(7/9 - 37*t) + 12/13*sin(10/3 - 36*t) + 16/9*sin(41/9 - 34*t) + 7/5*sin(26/9 - 33*t) + 11/10*sin(13/10 - 29*t) + 12/5*sin(25/6 - 27*t) + 50/17*sin(7/6 - 26*t) + 13/8*sin(55/18 - 25*t) + 23/8*sin(33/8 - 24*t) + 31/7*sin(4/3 - 23*t) + 19/7*sin(37/9 - 22*t) + 15/4*sin(22/5 - 21*t) + 37/7*sin(9/2 - 20*t) + 17/8*sin(14/3 - 19*t) + 3/4*sin(5/2 - 18*t) + 43/8*sin(16/5 - 17*t) + 37/4*sin(3/5 - 16*t) + 52/9*sin(14/9 - 15*t) + 11/4*sin(9/5 - 13*t) + 17/10*sin(13/3 - 11*t) + 201/7*sin(10/3 - 10*t) + 19/5*sin(49/16 - 9*t) + 101/8*sin(24/7 - 8*t) + 110/13*sin(95/24 - 7*t) + 112/3*sin(21/5 - 5*t) + 23/4*sin(19/6 - 4*t) + 435/11*sin(11/6 - 3*t) + 3992/13*sin(67/22 - t) - 1314/11*sin(2*t + 3/4) - 125/9*sin(6*t + 15/14) - 26/5*sin(12*t + 1/2) - 141/20*sin(14*t + 2/3) - 17/4*sin(28*t + 3/4) - 13/6*sin(30*t + 2/7) - 10/7*sin(31*t + 1/2) - 7/6*sin(32*t + 7/6) - 3/7*sin(35*t + 15/14) - 2/7*sin(49*t + 5/8) - 3/5*sin(52*t + 1/24) - 2/5*sin(54*t + 4/3) - 3/8*sin(55*t + 3/2) - 2/5*sin(60*t + 5/8) - 1/4*sin(67*t + 1/4) - 3/5*sin(68*t + 1/7) - 1/11*sin(75*t + 1) - 2/7*sin(76*t + 5/9) - 1/4*sin(80*t + 3/5) - 3/8*sin(83*t + 7/10) - 1/4*sin(94*t + 10/11) - 1/23*sin(98*t + 1/3) - 1/12*sin(99*t + 3/2) - 685/8;
    y = 1/5*sin(3 - 100*t) + 1/6*sin(7/3 - 99*t) + 2/7*sin(2/3 - 98*t) + 1/6*sin(19/6 - 97*t) + 1/7*sin(13/4 - 96*t) + 1/18*sin(31/7 - 95*t) + 1/14*sin(7/4 - 93*t) + 1/6*sin(12/5 - 92*t) + 2/9*sin(13/4 - 91*t) + 2/7*sin(29/14 - 90*t) + 1/6*sin(8/3 - 89*t) + 1/5*sin(23/24 - 88*t) + 1/7*sin(7/6 - 87*t) + 2/9*sin(23/7 - 86*t) + 1/6*sin(7/3 - 85*t) + 1/7*sin(1/5 - 83*t) + 1/4*sin(1/2 - 82*t) + 1/3*sin(10/3 - 81*t) + 1/5*sin(11/9 - 80*t) + 3/8*sin(23/9 - 79*t) + 1/4*sin(9/7 - 78*t) + 1/3*sin(23/22 - 77*t) + 1/7*sin(38/11 - 76*t) + 1/5*sin(11/5 - 75*t) + 1/6*sin(1/12 - 74*t) + 1/4*sin(7/5 - 72*t) + 3/8*sin(35/9 - 71*t) + 2/7*sin(15/7 - 70*t) + 6/11*sin(21/10 - 69*t) + 2/7*sin(5/9 - 68*t) + 1/2*sin(3/5 - 67*t) + 1/7*sin(37/10 - 66*t) + 2/9*sin(6/5 - 65*t) + 3/7*sin(17/7 - 64*t) + 6/7*sin(14/3 - 63*t) + 2/5*sin(13/4 - 62*t) + 2/7*sin(91/23 - 61*t) + 2/7*sin(8/5 - 60*t) + 9/10*sin(13/7 - 59*t) + 6/13*sin(1/14 - 58*t) + 1/10*sin(29/8 - 56*t) + 1/24*sin(23/7 - 55*t) + 23/22*sin(15/7 - 54*t) + 3/4*sin(17/5 - 52*t) + 6/11*sin(19/6 - 51*t) + 6/11*sin(7/10 - 50*t) + 2/3*sin(19/8 - 49*t) + 2/3*sin(11/7 - 48*t) + 1/7*sin(6/13 - 47*t) + 2/7*sin(50/11 - 46*t) + 4/7*sin(8/5 - 45*t) + 15/14*sin(4/3 - 44*t) + 19/20*sin(44/15 - 42*t) + 1/2*sin(12/5 - 41*t) + 4/7*sin(1 - 40*t) + 1/7*sin(27/7 - 39*t) + 13/5*sin(11/7 - 38*t) + 13/10*sin(16/5 - 36*t) + 7/5*sin(8/9 - 35*t) + 8/9*sin(19/20 - 33*t) + 19/13*sin(25/6 - 32*t) + 9/7*sin(10/7 - 31*t) + sin(13/5 - 30*t) + 21/10*sin(11/10 - 28*t) + 3/2*sin(13/8 - 27*t) + 3/4*sin(23/8 - 26*t) + 11/4*sin(17/6 - 25*t) + 19/7*sin(11/4 - 23*t) + 11/4*sin(17/6 - 21*t) + 17/4*sin(16/11 - 20*t) + 29/11*sin(3/8 - 18*t) + 1/8*sin(47/12 - 17*t) + 58/7*sin(11/8 - 16*t) + 27/2*sin(13/4 - 15*t) + 53/8*sin(4/5 - 14*t) + 186/11*sin(22/9 - 12*t) + 97/13*sin(12/5 - 10*t) + 191/19*sin(1/34 - 9*t) + 92/7*sin(21/5 - 8*t) + 17*sin(24/7 - 7*t) + 547/7*sin(6/5 - 6*t) + 805/9*sin(13/6 - 5*t) + 266/3*sin(18/7 - 4*t) + 81/5*sin(2/7 - 3*t) + 1138/7*sin(2/3 - 2*t) + 956/3*sin(22/5 -t) - 169/12*sin(11*t) - 2*sin(29*t) - 39/4*sin(13*t + 11/8) - 30/7*sin(19*t + 6/5) - 13/8*sin(22*t + 7/6) - 49/10*sin(24*t + 4/9) - 4/3*sin(34*t + 7/10) - 5/9*sin(37*t + 4/7) - 1/3*sin(43*t + 1/2) - 4/7*sin(53*t + 5/4) - 1/11*sin(57*t + 1/9) - 3/4*sin(73*t + 13/9) - 1/5*sin(84*t + 2/3) - 1/10*sin(94*t + 5/7) - 28/3;
    % Scale between 0 and 1
    x = (x - min(x)) / (max(x) - min(x));
    y = (y - min(y)) / (max(y) - min(y));
    % Define rotation functions
    rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
    roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
    rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;
    % Also do a bit of squeezing and stretching
    T = [0.9153, -0.4796; -0.2268, -0.5838]; %-1 + 2*rand(2)
    % Do some rotations to create a 3d trajectory
    X_n = rotz(deg2rad(15))*roty(deg2rad(-15))*rotx(deg2rad(45))*...
        [x; y; 0.1*rand(size(x))];
    X_m = rotz(deg2rad(0))*roty(deg2rad(0))*rotx(deg2rad(-45))*...
        [T*[x; y;]; 0.1*rand(size(x))];
end