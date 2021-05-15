clear;
close all;
clc;
addpath('../src');

% number of used threads
maxNumCompThreads('automatic'); 

eps_deg  = 1.e-3; % eigenvalue degeneration tolerance in [MeV]
N_lambda = 5;     % number of artificial lambda iteration

% reading h and Delta matrices corresponding to O16 from HFBTHO
h     = struct2array( load( '../data/O16/HFBTHO/h.mat'     ) );
Delta = struct2array( load( '../data/O16/HFBTHO/Delta.mat' ) );
n     = length(h);

% lambdas are generated to simulate the root-finding in lambda iterations
lambdas = 10 * randn( 1 , N_lambda ); % [MeV]

% radnom permutation of rows and colums of h and Delta
P     = randperm(n);
h     = h(P,P);
Delta = Delta(P,P);





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
    
    [Q,E] = eig( [h,Delta;Delta,-h] - [lambdas(i)*eye(n),zeros(n);zeros(n),-lambdas(i)*eye(n)] );
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










