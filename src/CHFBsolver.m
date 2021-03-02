classdef CHFBsolver
    % CHFBsolver
    %
    % How to use:
    % 1.) Instantiate a class with given matrices h,Delta and degeneration
    % tolerance eps_deg using obj = CHFBsolver( h , Delta , eps_deg ).
    %
    % 2.) Once instantiated, the method obj.lambdaIteration( lambda )
    % returns the function value f(lambda) for given lambda.
    
    properties( Access = private )
        eps_deg;
        h;
        Delta;
        X; % X = ( h + 1j*Delta )*( h - 1j*Delta ) = ( h^2 + Delta^2 ) + 1j*( Delta*h - h*Delta )
    end
    
    methods( Access = public )
    
        function obj = CHFBsolver( h , Delta , eps_deg )
            % CHFBsolver( h , Delta , eps_deg )
            % Construct an instance of this class.
            % h and Delta must be real symmetric nxn matrices.
            % eps_deg must be real positive number.
            
            assert( size(h,1) == size(h,2)      , 'h must be square matrix'         );
            assert( issymmetric(h) && isreal(h) , 'h must be real symmetric matrix' );
            
            assert( size(Delta,1) == size(Delta,2)      , 'Delta must be square matrix'         );
            assert( issymmetric(Delta) && isreal(Delta) , 'Delta must be real symmetric matrix' );
            
            assert( size(h,1) == size(Delta,1) , 'h and Delta must be matrices of the same order' );

            assert( isreal(eps_deg) && size(eps_deg,1)==1 && size(eps_deg,2)==1 && eps_deg>0 , 'eps_deg must be positive real number' );

            
            obj.eps_deg = eps_deg;
            
            obj.h     = h;
            obj.Delta = Delta;
            
            % step 0.)
            ReX   = h^2 + Delta^2;
            ImX   = Delta*h;
            ImX   = ImX - transp(ImX);
            X     = complex( ReX , ImX );
            obj.X = 0.5 * ( X + ctranspose(X) );
            
            
            return;
        end
        
        function f = lambdaIteration( obj , lambda )
            % f = obj.lambdaIteration( lambda )
            % Returns f(lambda) for given lambda during the lambda
            % iterations for instantianted object
            % obj with parameters h, Delta and eps_deg.
            
            
            n = size(obj.h,1);
            
            h_lambda = obj.h - lambda*eye(n,n);
            
            % step 1.)
            hdhd_lambda = obj.X - 2*lambda*obj.h + (lambda^2)*eye(n,n);
            
            % step 2.)
            [Q,E] = eig( hdhd_lambda );
            [Q,E] = CHFBsolver.sortem( Q , E );
            E  = diag(E);
            Q1 = real(Q);
            Q2 = imag(Q);
            
            hQ1     = h_lambda  * Q1;
            hQ2     = h_lambda  * Q2;
            DeltaQ1 = obj.Delta * Q1;
            DeltaQ2 = obj.Delta * Q2;

            
            i = 1;
            while( i <= n )
                ev = E(i);
                
                % step 3.)
                j = i;
                while( j<=n && abs(ev-E(j))<obj.eps_deg )
                    j = j + 1;
                end
                j = j - 1;
                % E(i) = E(i+1) = ... = E(j-1) = E(j)


                % step 4.)
                e = + transp(Q1(:,i:j))*hQ1(:,i:j)     - transp(Q2(:,i:j))*hQ2(:,i:j)      ...
                    + transp(Q1(:,i:j))*DeltaQ2(:,i:j) + transp(Q2(:,i:j))*DeltaQ1(:,i:j);

                d = + transp(Q1(:,i:j))*DeltaQ1(:,i:j) - transp(Q2(:,i:j))*DeltaQ2(:,i:j)  ...
                    - transp(Q1(:,i:j))*hQ2(:,i:j)     - transp(Q2(:,i:j))*hQ1(:,i:j);

                e = 0.5 * ( e + transp(e) );
                d = 0.5 * ( d + transp(d) );


                % step 5.)
                [qq,ee] = eig( [ e , d ; d , -e ] );
                [qq,~ ] = CHFBsolver.sortem( qq , ee );

                nb = j-i+1;
                C = qq(    1 : nb    , 1 : nb );
                S = qq( nb+1 : nb+nb , 1 : nb );
                
                
                % step 6.)
                Q2(:,i:j) = Q2(:,i:j)*C + Q1(:,i:j)*S;

                
                i = j + 1;    
            end

            
            % step 7.)
            f = norm( Q2 , 'fro' )^2;
            return;
        
        end

    end

    methods( Static )
    
        function [ P , D ] = sortem( P , D )
            [D,ind] = sort( diag(D) , 'descend' );
            D = diag(D);
            P = P(:,ind);
            return;
        end
        
    end
    
end

