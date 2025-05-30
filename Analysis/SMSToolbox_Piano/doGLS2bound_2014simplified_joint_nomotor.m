% actual bGLS!!!
% for JOINT model, including boundaries for alpha to improve estimates

function [x,sM,sT,LL]=doGLS2bound_2014simplified_joint_nomotor(A,b,L,H)
assert(size(A,2)==3);
ITER=20;

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
    
    if (1-z2(2))>0
        if z2(3)<L*(1-z2(2))
            z2(3)=L*(1-z2(2));
        end
        if z2(3)>H*(1-z2(2))
            z2(3)=H*(1-z2(2));
        end
    else
        if z2(3)<H*(1-z2(2))
            z2(3)=H*(1-z2(2));
        end
        if z2(3)>L*(1-z2(2))
            z2(3)=L*(1-z2(2));
        end
    end
end
%z2
% output the residual noise covriance
x=z2;

LL = log2(mvnpdf(A*x-b,zeros(size(b)),CC));

end