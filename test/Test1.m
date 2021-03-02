clear;
close all;
clc;
addpath('../src');

% number of used threads
maxNumCompThreads('automatic'); 

eps_deg  = (1.e-3)^2; % eigenvalue degeneration tolerance in [MeV^2]
N_lambda = 5;         % number of artificial lambda iteration

% Hamiltonian with positive spectrum given below will be constructed
n = 3000;
PositiveSpectrum = 1 + linspace( 0 , 10 , n );       
[ h , Delta ] = randHamiltonianMatix( PositiveSpectrum );

% lambdas are generated to simulate the root-finding in lambda iterations
lambdas = 10 * randn( 1 , N_lambda ); % [MeV]










% proposed method
tic;
TargetValues1 = zeros( 1 , N_lambda );

chfbsolver = CHFBsolver( h , Delta , eps_deg );
for i = 1 : N_lambda
    TargetValues1(i) = chfbsolver.lambdaIteration( lambdas(i) );
end
t1 = toc;






% straightforward method
tic;
TargetValues2 = zeros( 1 , N_lambda );

for i = 1 : N_lambda
    
    [Q,E] = eig( [h,Delta;Delta,-h] - lambdas(i)*[eye(n),zeros(n);zeros(n),-eye(n)] );
    [Q,E] = CHFBsolver.sortem( Q , E );
     
    TargetValues2(i) = norm( Q( n+1:n+n , 1:n ) , 'fro' )^2;
end
t2 = toc;








% short report
thr = maxNumCompThreads;
fprintf('Used number of threads:      %d.\n' , thr );
fprintf('Number of lambda iterations: %d.\n' , N_lambda );
fprintf('Problem size: n = %d, i.e. H is %dx%d symmetric Hamiltonian matrix.\n' , n , 2*n , 2*n );
fprintf('\n');
fprintf('Time elapsed (proposed method):        %f s.\n' , t1 );
fprintf('Time elapsed (straightforward method): %f s.\n' , t2 );
fprintf('Time elapsed ratio:                    %f.\n'   , t2/t1 );
fprintf('\n');
fprintf('Comparison of target values: \n');
for i = 1 : N_lambda
    lambda = lambdas(i);
    f1     = TargetValues1(i);
    f2     = TargetValues2(i);
    diff   = abs(f1-f2);
    fprintf('lambda = %10.5f MeV, f1 = %10.5f, f2 = %10.5f., |f1-f2| = %10.5e\n' , lambda , f1 , f2 , diff );
end










function [ h , Delta ] = randHamiltonianMatix( PositiveSpectrum )
    
    n = length(PositiveSpectrum);
    
    [U,~] = qr( randn(n,n) );
    [V,~] = qr( randn(n,n) );
    
    [thetas,~] = cart2pol( randn(n,1) , randn(n,1) );
    C = diag( cos(thetas) );
    S = diag( sin(thetas) );

    Q1 = U*C*transp(V);
    Q2 = U*S*transp(V);
    
    H = transp([Q1,-Q2;Q2,Q1]) * [ +diag(PositiveSpectrum) ,  zeros(n,n)                ; ...
                                    zeros(n,n)             , -diag(PositiveSpectrum)  ] * ...
               [Q1,-Q2;Q2,Q1];
    
    
    h     = H( 1:n ,   1:n   );
    Delta = H( 1:n , n+1:n+n );
    
    h     = 0.5 * ( h     + transp(h)     );
    Delta = 0.5 * ( Delta + transp(Delta) );
    
    return;
end