%% Main Code to solve system of equations using gaussian elimination, calls on other code: {gauss_elimination}

% Define Matrix A
A = [pi, -sqrt(2), -1, 1, 0;
     exp(1), -1, 1, 2, 1;
     1, 1, -sqrt(3), 1, 2;
     -1, -1, 1, -sqrt(5), 3];

% Call the function gauss_elimination
solution = gauss_elimination(A);

% Display the solution solves for x1...x4
disp('Solutions for X by Gaussian Elimination:');
disp(solution);


% Solve system of equations using matrices to check the gaussian method
A=[pi, -sqrt(2), -1, 1; 
    exp(1), -1, 1, 2;
    1, 1, -sqrt(3), 1;
    -1, -1, 1, -sqrt(5);  ];

B=[0;1;2;3];


x=inv(A)*B;
disp('Solutions for X by matrix multiplication:');
disp(x);
% Code returns same values----