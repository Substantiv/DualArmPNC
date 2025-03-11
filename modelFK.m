function [p, pc, R] = modelFK(x, gamma, q)
% modelFK  Computes the forward kinematics of the robot
% Inputs:
%   x      - Scalar, the position of the robot in the global coordinate frame along the x-axis.
%   gamma  - Scalar, the rotational angle of the robot (rotation around the z-axis).
%   q      - 6x1 vector, joint angles q = [q11, q12, q13, q21, q22, q23], 
%            representing the rotation angles for each joint.
% Outputs:
%   p      - 3x6 matrix, representing the positions of each joint. 
%   pc     - 3x6 matrix, representing the positions of the center of mass (COM) for each joint. 
%   R      - 3x7 matrix, representing the rotation matrices for each joint.

% Define joint angles from input q
q11 = q(1); q12 = q(2); q13 = q(3);
q21 = q(4); q22 = q(5); q23 = q(6);

% Define masses
mt = 0.2; m0 = 10.18; m1 = 0.4; m2 = 1; m3 = 0.6; m4 = m1; m5 = m2; m6 = m3;
m = [m1 m2 m3 m4 m5 m6];

% Link lengths
a1 = 1.5; a2 = 1.5; a3 = 1.5; a4 = a1; a5 = a2; a6 = a3;
b = 0.258;

% Define inertia matrices for each link
I0xx = 10.4; I0yy = 10.4; I0zz = 10.4;
I1xx = 2e-4; I1yy = 2e-4; I1zz = 2e-4;
I2xx = 3.5e-4; I2yy = 3.5e-4; I2zz = 2e-4;
I3xx = 3.5e-4; I3yy = 3.7e-4; I3zz = 2e-4;

I1 = [I1xx 0 0; 0 I1yy 0; 0 0 I1zz];
I2 = [I2xx 0 0; 0 I2yy 0; 0 0 I2zz];
I3 = [I3xx 0 0; 0 I3yy 0; 0 0 I3zz];
I4 = I1; I5 = I2; I6 = I3;
It = [2.08e-4, 0, 0; 0, 2.08e-4, 0; 0, 0, 3.33e-4];


% Rotation matrix (Inertial Coordinate -> Base Coordinate)
Rb_I = [1, 0, 0;
        0, cos(gamma), -sin(gamma);
        0, sin(gamma), cos(gamma)];

% Transformation matrix (Inertial Coordinate -> Base Coordinate)
Tb_I = [1, 0, 0, x;
        0, cos(gamma), -sin(gamma), 0.5;
        0, sin(gamma), cos(gamma), 0.5;
        0, 0, 0, 1];

% Transformation matrices for each link
Ta1 = [-1, 0, 0; 0, 1, 0; 0, 0, -1];

% Rotation matrix (Base Coordinate -> Left1 Coordinate)
R1_b1 = Ta1 * [cos(q11), -sin(q11), 0;
               sin(q11), cos(q11),  0;
               0, 0, 1];

% Position of the first joint in base frame
P1_b1 = Ta1 * [a1 * cos(q11); a1 * sin(q11); 0];

% Transformation matrix (Base Coordinate -> Left1 Coordinate)
A1_b1 = [R1_b1, P1_b1; zeros(1, 3), 1];

% Transformation matrices for subsequent links
A2_1 = [cos(q12), -sin(q12), 0, a2 * cos(q12);
        sin(q12), cos(q12), 0, a2 * sin(q12);
        0, 0, 1, 0;
        0, 0, 0, 1];

A3_2 = [cos(q13), -sin(q13), 0, a3 * cos(q13);
        sin(q13), cos(q13), 0, a3 * sin(q13);
        0, 0, 1, 0;
        0, 0, 0, 1];

Ae1_3 = [1, 0, 0, 0;
          0, 1, 0, 0;
          0, 0, 1, b / 2;
          0, 0, 0, 1];

% Combining transformations for all joints
A1_0 = Tb_I * A1_b1;
A2_0 = A1_0 * A2_1;
A3_0 = A2_0 * A3_2;
Ae1_0 = A3_0 * Ae1_3;

% Positions of joints
p1 = A1_0(1:3, 4);
p2 = A2_0(1:3, 4);
p3 = A3_0(1:3, 4);
pe1 = Ae1_0(1:3, 4);

% Transformation matrix for the Center of Mass (COM)
A1c_1 = [1, 0, 0, a1 / 2;
          0, 1, 0, 0;
          0, 0, 1, 0;
          0, 0, 0, 1];
A2c_2 = [1, 0, 0, a2 / 2;
          0, 1, 0, 0;
          0, 0, 1, 0;
          0, 0, 0, 1];
A3c_3 = [1, 0, 0, a3 / 2;
          0, 1, 0, 0;
          0, 0, 1, 0;
          0, 0, 0, 1];

% Transformation matrices for COM
A1c_0 = A1_0 * A1c_1;
A2c_0 = A2_0 * A2c_2;
A3c_0 = A3_0 * A3c_3;

% Rotation matrices for COM
R1 = A1c_0(1:3, 1:3);
R2 = A2c_0(1:3, 1:3);
R3 = A3c_0(1:3, 1:3);

% Rotation matrices of joints
R1_0 = A1_0(1:3, 1:3);
R2_0 = A2_0(1:3, 1:3);
R3_0 = A3_0(1:3, 1:3);

% Transformation matrices for subsequent joints (4-6)
R4_b2 = [cos(q21), -sin(q21), 0;
          sin(q21), cos(q21), 0;
          0, 0, 1];

P4_b2 = [a4 * cos(q21); a4 * sin(q21); 0];

A4_b2 = [R4_b2, P4_b2;
          zeros(1, 3), 1];

A5_4 = [cos(q22), -sin(q22), 0, a5 * cos(q22);
        sin(q22), cos(q22), 0, a5 * sin(q22);
        0, 0, 1, 0;
        0, 0, 0, 1];

A6_5 = [cos(q23), -sin(q23), 0, a6 * cos(q23);
        sin(q23), cos(q23), 0, a6 * sin(q23);
        0, 0, 1, 0;
        0, 0, 0, 1];

Ae2_6 = [1, 0, 0, 0;
          0, 1, 0, 0;
          0, 0, 1, b / 2;
          0, 0, 0, 1];

% Final transformation matrices
A4_0 = Tb_I * A4_b2;
A5_0 = A4_0 * A5_4;
A6_0 = A5_0 * A6_5;
Ae2_0 = A6_0 * Ae2_6;

% Positions of joints 4-6
p4 = A4_0(1:3, 4);
p5 = A5_0(1:3, 4);
p6 = A6_0(1:3, 4);
pe2 = Ae2_0(1:3, 4);

% Storing positions in arrays
p = [p1 p2 p3 p4 p5 p6];

% Center of mass positions for joints
A4c_4 = [1, 0, 0, a1 / 2; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A5c_5 = [1, 0, 0, a2 / 2; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
A6c_6 = [1, 0, 0, a3 / 2; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];

A4c_0 = A4_0 * A4c_4;
A5c_0 = A5_0 * A5c_5;
A6c_0 = A6_0 * A6c_6;

% Rotation matrices for joints 4-6
R4 = A4c_0(1:3, 1:3);
R5 = A5c_0(1:3, 1:3);
R6 = A6c_0(1:3, 1:3);

% Rotation matrices for joints 4-6 in base frame
R4_0 = A4_0(1:3, 1:3);
R5_0 = A5_0(1:3, 1:3);
R6_0 = A6_0(1:3, 1:3);

% Center of mass positions for joints
p1c = A1c_0(1:3, 4);
p2c = A2c_0(1:3, 4);
p3c = A3c_0(1:3, 4);
p4c = A4c_0(1:3, 4);
p5c = A5c_0(1:3, 4);
p6c = A6c_0(1:3, 4);

pc = [p1c p2c p3c p4c p5c p6c];

% Final rotation matrices
Ae1_b = [1, 0, 0, 0;
          0, 1, 0, 0;
          0, 0, 1, b / 2;
          0, 0, 0, 1];
Ab_0 = A3_0 * Ae1_3 * Ae1_b;
pb = Ab_0(1:3, 4);
R7 = Ab_0(1:3, 1:3);

% Final rotation matrices
R = [R1 R2 R3 R4 R5 R6 R7];

end
