eq=input('No of equations:\n');
A=input('Input Matrix : \n');
it=input('Max. Iteration :\n');
er=input('Max Error :\n');
x1=input('Input 1 for Direct Power Method:\nInput 2 for Inverse Power Method:\nInput 3 for Shifted Power Method:\nInput 4 for QR Method:\n');

if (x1==1)
    disp('Direct power method')
    x=ones(eq,1);
    iter=0;
    err=100000;
    e1=1;
    while(err>er)
        b=A*x;
        eigen_value=max(abs(b));
        eigen_vector=b/eigen_value;
        x=eigen_vector;
        e2=eigen_value;
        err=abs((100*(e2-e1))/e2);
        e1=e2;
        
        iter=iter+1;
        if (iter>it)
            break
        end
    end
    
    disp("Eigenvalue");
    disp(eigen_value);
    disp("Eigenvector");
    disp(x);
    disp('Iterations ');
    disp(iter);
end
if (x1==2)
    disp('Inverse power method')
    x=ones(eq,1);
    iter=0;
    err=100000;
    e1=1;
    while(err>er)
        b=inv(A)*x;
        eigen_value=max(abs(b));
        eigen_vector=b/eigen_value;
        x=eigen_vector;
        e2=eigen_value;
        err=abs((100*(e2-e1))/e2);
        e1=e2;
        
        iter=iter+1;
        if (iter>it)
            break
        end
    end
    
    disp("Eigenvalue");
    disp(1/eigen_value);
    disp("Eigenvector");
    disp(x);
    disp('Iterations');
    disp(iter);
end
if (x1==3)
    s=input('Value of shift: \n');
    disp('Shifted power method')
    x=ones(eq,1);
    iter=0;
    err=100000;
    e1=1;
    while(err>er)
        
        b=inv(A-s*eye(eq))*x;
        eigen_value=max(abs(b));
        eigen_vector=b/eigen_value;
        x=eigen_vector;
        e2=eigen_value;
        err=abs((100*(e2-e1))/e2);
        e1=e2;
        
        iter=iter+1;
        if (iter>it)
            break
        end
    end
    
    disp("Eigenvalue");
    disp((1/eigen_value)+s);
    disp("Eigenvector");
    disp(x);
    disp('Iterations');
    disp(iter);
end

if(x1==4)
    n=eq;
    err=er;
    err=err/100;
    if (rank(A)<n)
        fprintf("\n\nMatrix is singular , try again with another matrix!");
        return;
    end
    m=n;
    ops=0;
    diagA=zeros(n,1);
    errA=ones(n,1);
    for i=1:n
        diagA(i)=A(i,i);
    end
    while(ops<9999)
        z = zeros(n,1);
        Q = zeros(m,n);
        R = zeros(n,n);
        for i = 1:m
            Q(i,1) = A(i,1) / norm(A(:,1));
        end
        for i = 2:n
            sum = zeros(m,1);
            for k = 1:i-1
                sum = sum + (A(:,i)'*Q(:,k))*Q(:,k);
            end
            z = A(:,i) - sum;
            for j = 1:m
                Q(j,i) = z(j)/norm(z);
            end
        end
        
        Qinv= inv(Q);
        R=Qinv*A;
        if(ops>0)
            for i=1:n
                errA(i)= abs((A(i,i)-diagA(i))/A(i,i));
            end
        end
        for i=1:n
            diagA(i)=A(i,i);
            if(norm(errA)<=err) 
                break; 
            end
        end
        ops=ops+1;
        A= R*Q ;
    end
    fprintf("\nEigenvalues\n");
    disp(diagA);
end