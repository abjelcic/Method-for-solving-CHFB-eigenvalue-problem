clear;
close all;


%ova linija br. threadova na optimalno
maxNumCompThreads('automatic'); 



degentol = 1.e-6; 
N_lambda = 5; %number of artificial lambda iteration

%Hamiltonian with positive spectrum given below will be constructed
nd = 1000;
Spectrum = [ -1*ones(1,nd) , ...
             -2*ones(1,nd) , ...
              3*ones(1,nd) ];

                







Lambdas = randn(1,N_lambda);

n = length(Spectrum);
h = RandSymmetricMatrix( Spectrum );
Delta = 1.e-12*ones(n,n);
H = [ h , Delta ; Delta , -h ];





tic
%% Moja metoda

hd   = h + 1j*Delta;
hdhd = hd * hd';
hdhd = 0.5*( hdhd + hdhd' );

TargetValues1 = [];
for lambda = Lambdas
    
    h_lambda = h - lambda*eye(n,n);
    hdhd_lambda = hdhd - 2*lambda*h + (lambda^2)*eye(n,n);
    
    [Q,D] = eig(hdhd_lambda);
    [Q,D] = sortem(Q,D);
    D = diag(D);
    Q1 = real(Q);
    Q2 = imag(Q);

    % Odkomentiraj ovo i uvjeri se da ed matrica ima blokdiag formu za
    % degenerirani slucaj
    % ed = [Q1,-Q2;Q2,Q1]' * [h,Delta;Delta,-h] * [Q1,-Q2;Q2,Q1];
    % surf(log10(abs(ed)))
    % return;

    hQ1     = h_lambda*Q1;
    hQ2     = h_lambda*Q2;
    DeltaQ1 = Delta*Q1;
    DeltaQ2 = Delta*Q2;


    i = 1;
    while i <= n
        ev = D(i);

        j = i;
        while j<=n && abs(ev-D(j))<degentol
            j = j + 1;
        end
        j = j - 1;

        % sad su od i-te do j-te (uljkucivo) sv. vrijednosti degenerirane!



        % tu bi se mozda dalo nesto ustedi operacija, ali za sad ne brinem
        e = + Q1(:,i:j)'*hQ1(:,i:j)     - Q2(:,i:j)'*hQ2(:,i:j) ...
            + Q1(:,i:j)'*DeltaQ2(:,i:j) + Q2(:,i:j)'*DeltaQ1(:,i:j);

        d = + Q1(:,i:j)'*DeltaQ1(:,i:j) - Q2(:,i:j)'*DeltaQ2(:,i:j) ...
            - Q1(:,i:j)'*hQ2(:,i:j)     - Q2(:,i:j)'*hQ1(:,i:j);

        e = 0.5*(e+e');
        d = 0.5*(d+d');



        [qq,ee] = eig( [ e , d ; d , -e ] );
        [qq,ee] = sortem(qq,ee);

        % qq ima oblik qq=[C,-S;S,C], gdje su u slucaju nn=1 stvarno
        % C i S obicni brojevi cos i sin
        nn = j-i+1;
        C = qq(1:nn,1:nn);
        S = qq(nn+1:nn+nn,1:nn);

        Q2(:,i:j) = Q2(:,i:j)*C + Q1(:,i:j)*S;

        i = j + 1;    
    end

    f1 = norm( Q2 , 'fro' )^2;
    TargetValues1 = [ TargetValues1 ; f1 ];
end
%%
t1 = toc;



tic;
%% Naivna metoda
TargetValues2 = [];
for lambda = Lambdas
    [Q,E] = eig( H - lambda*[eye(n),zeros(n);zeros(n),-eye(n)] );
    [Q,E] = sortem(Q,E);
    f2 = norm( Q( n+1:end , 1:n ) , 'fro' )^2; 
    TargetValues2 = [ TargetValues2 ; f2 ];
end
%%
t2 = toc;





thr = maxNumCompThreads;
fprintf('Koristim %d threadova procesora.\n', thr);
fprintf('Velicina problema n=%d, dakle %dx%d Hamiltonijan je u pitanju!\n\n', n, 2*n, 2*n );
fprintf('Napravljene su %d lambda iteracije.\n', N_lambda);
fprintf('Moja   metoda: %f sec\n', t1   );
fprintf('Naivna metoda: %f sec\n', t2   );
fprintf('Moj je %f puta brzi! \n', t2/t1);
fprintf('\n');
fprintf('\n');
fprintf('Usporedba target values sa obje metode: \n');
[ TargetValues1 , TargetValues2 , TargetValues1-TargetValues2 ]






function A = RandSymmetricMatrix( Spectrum )
    Q = RandOrhogonalMatrix( length(Spectrum) );
    A = Q'*diag(Spectrum)*Q;
    A = 0.5*(A+A');
    return;
end

function Q = RandOrhogonalMatrix( n )
    [Q,~,~] = svd( randn(n,n) );
    return;
end

function [ P2 , D2 ]=sortem( P , D )
    D2 = diag(sort(diag(D),'descend'));
    [~,ind] = sort(diag(D),'descend');
    P2 = P(:,ind);
end

function H = RandHamiltonianMatix( PositiveSpectrum )
    
    n = length(PositiveSpectrum);
    [U,~,~] = svd(randn(n,n));
    [V,~,~] = svd(randn(n,n));
    
    C = sqrt( rand(n,1) );
    S = sqrt( ones(n,1) - C.^2 );
    C = diag(C);
    S = diag(S);

    Q1 = U*C*V';
    Q2 = U*S*V';
    
    H = [Q1,-Q2;Q2,Q1]' * [ diag(PositiveSpectrum) , zeros(n,n) ; ...
                            zeros(n,n)             , -diag(PositiveSpectrum)  ] * ...
        [Q1,-Q2;Q2,Q1];
    
    H = 0.5*(H+H');
    
    return;
end

function A = randsymm( n )
    A = zeros(n,n);
    for i = 1 : n
        for j = i+1 : n
            A(i,j) = randn();
            A(j,i) = A(i,j);
        end
        A(i,i) = randn();
    end
    return;
end

