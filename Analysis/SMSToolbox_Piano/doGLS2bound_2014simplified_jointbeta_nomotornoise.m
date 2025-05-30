% actual bGLS!!!
% for JOINT model, including boundaries for beta to improve estimates

function [x,sM,sT,LL]=doGLS2bound_2014simplified_jointbeta_nomotornoise(A,b,L,H)
assert(size(A,2)==3);
ITER=1;
% important: A and b must have 0 mean!!!
z=mldivide (A,b); % start with normal regression 
z2=z;

% iterate to find the maximal likelihood estimator.

for I=1:ITER, % this is enough iteration, if not accurate, increase
%z2
    % calculate the residual noise covariance    
    d=A*z2-b;
    K=cov(d(1:(end-1)),d(2:end));
    

    %K
    K11=(K(1,1)+K(2,2))/2;
    K12=K(1,2);

    if K12>0
        sM=0;
    else
        sM=sqrt(-K12);
    end
    
    if K11<3*(sM^2)
        sM=sqrt(K11/3);
    end
    sM=0; %NOOOTE I KILLED HERE MOTOR VARIANCE.
    
    sT=sqrt(K11 - sM^2);
    
    
    
    K11=(sT^2)+ 2*(sM^2);
    K12=-sM^2;


    % calculate GLS  with known covariance
    NN=length(b);
    CC=diag(K11*ones(1,NN),0)+ diag(K12*ones(1,NN-1),1) + diag(K12*ones(1,NN-1),-1);
    iC=inv(CC);
    z2=inv((A')*iC*A)*(A')*iC*b;
    %z2';
    if (z2(2))<0.1
        z2(2)=0.1;
    end
    if (z2(2))>0.9
        z2(2)=0.9;
    end

    if (z2(2))>0
        if z2(1)<L*(z2(2))
            z2(1)=L*(z2(2));
        end
        if z2(1)>H*(z2(2))
            z2(1)=H*(z2(2));
        end
    else
        if z2(1)<H*(z2(2))
            z2(1)=H*(z2(2));
        end
        if z2(1)>L*(z2(2))
            z2(1)=L*(z2(2));
        end
    end
end
%z2
% output the residual noise covriance
x=z2;

% Residual vector
residual = A*x - b;

% Mean vector (should match size)
mu = zeros(size(b));

% Validate CC: symmetry & positive definiteness
if ~isequal(CC, CC') || any(eig(CC) <= 0)
    warning('⚠️ Covariance matrix was not positive definite — regularizing.');
    CC = CC + eye(size(CC)) * 1e-6;
end

% Validate residual size matches mu and CC
n = length(mu);
if ~isequal(size(residual), [n,1])
    error('❌ Residual and mean vector sizes do not match.');
end
if ~isequal(size(CC), [n,n])
    error('❌ Covariance matrix size does not match residual vector.');
end

% Calculate pdf safely
try
    pdf_val = mvnpdf(residual, mu, CC);
catch err
    warning('❌ mvnpdf failed: %s', err.message);
    pdf_val = 1e-300 * (1 + rand()*9);  % tiny positive value to avoid -Inf in LL
end

% Final floor to avoid log2(0) or log2(negative)
if ~isfinite(pdf_val) || pdf_val <= 0
    pdf_val = 1e-300 * (1 + rand()*9);
end

LL = log2(pdf_val);