clear;
fid = fopen("INPUT.txt");
line =fgetl(fid);
size = sscanf(line, '%f ');
A = zeros(size,size+1);
for i=1:1:size
    line = fgetl(fid);
    A(i,1:1:size+1) = sscanf(line, '%f ');
end
fprintf("Enter the number in front of the method you want to use\n");
fprintf("1. Gauss Elimination(without pivoting)\n2. Gauss Elimination(with pivoting)\n3. GE (with scaling and pivoting)\n4. LU decomposition by using GE (without pivoting)\n5. LU decomposition by using GE (with pivoting)\n6. LU decomposition by Crout's method (without pivoting)\n7. Cholesky decomposition (for symmetric positive definite matrix)\n");
method = input("");

if method == 1
    for i=1:1:size-1
        if A(i,i)~=0
            for j=i+1:1:size
                A(j,:) = A(j,:)-(A(j,i)/A(i,i))*A(i,:);
            end
        else
            for k=i+1:1:size
                if A(k,i)~=0
                    B = A(k,:);
                    A(k,:) = A(i,:);
                    A(i,:) = B;
                    break;
                end
            end
            for j=i:1:size
                A(j,:) = A(j,:)-(A(j,i)/A(i,i))*A(i,:);
            end
        end
    end
    X = zeros(size,1);
    X(size,1) = A(size,size+1)/A(size,size);
    for i=size-1:-1:1
        for j=size:-1:i+1
            X(i,1) = A(i,size+1)-X(j,1)*A(i,j);
        end
        X(i,1) = X(i,1)/A(i,i);
    end
    fprintf("X:\n");
    X;
    filename = "output1.txt";
    outf = fopen (filename, "w");
    
    disp( X);
    fclose(outf);
end

if method == 2
    for i=1:1:size-1
        if A(i,i)~=0
            A(i,:) = A(i,:)/A(i,i);
            for j=i+1:1:size
                A(j,:) = A(j,:)-A(j,i)*A(i,:);
            end
        else
            for k=i+1:1:size
                if A(k,i)~=0
                    B = A(k,:);
                    A(k,:) = A(i,:);
                    A(i,:) = B;
                    break;
                end
            end
            A(i,:) = A(i,:)/A(i,i);
            for j=i:1:size
                A(j,:) = A(j,:)-A(j,i)*A(i,:);
            end
        end
    end
    X = zeros(size,1);
    X(size,1) = A(size,size+1)/A(size,size);
    for i=size-1:-1:1
        for j=size:-1:i+1
            X(i,1) = A(i,size+1)-X(j,1)*A(i,j);
        end
        
        X(i,1) = X(i,1)/A(i,i);
    end
    fprintf("X:\n");
    X;
    filename = "output1.txt";
    outf = fopen (filename, "w");
    disp(X);
    fprintf("Permutation Matrix:\n");
    disp(A);
    fclose(outf);
end


if method == 6
    U = zeros(size,size);
    L = zeros(size,size);
    for j=1:1:size
        L(j,1)=A(j,1);
    end
    U(1,:)=A(1,1:size)/L(1,1);
    U(1,1)=1;
    for k=2:size
        for j=2:size
            for i=j:size
                s1=0;
                for l = 1:j-1
                    s1 = s1 + L(i,l)*U(l,j);
                end
                L(i,j)=A(i,j)- s1;
            end
            s2=0;
            for l = 1:k-1
                s2 = s2 + L(k,l)*U(l,j);
            end
            U(k,j)=(A(k,j)- s2)/L(k,k);
        end
    end
    Y = zeros(size,1);
    Y(1,1) = A(1,size+1)/L(1,1);
    for i=2:1:size
        for j=1:1:i-1
            Y(i,1) = A(i,size+1)-Y(j,1)*L(i,j);
        end
        Y(i,1) = Y(i,1)/L(i,i);
    end
    X = zeros(size,1);
    X(size,1) = Y(size,1)/U(size,size);
    for i=size-1:-1:1
        for j=size:-1:i+1
            X(i,1) = Y(i,1)-X(j,1)*U(i,j);
        end
        X(i,1) = X(i,1)/U(i,i);
    end
    L;
    U;
    X;
    filename = "output1.txt";
    outf = fopen(filename, "w");
    fprintf("L:\n");
    disp(L);
    fprintf("U:\n");
    disp( U);
    fprintf("X:\n");
    disp( X);
    fclose(outf);
end

if method == 7
    L = zeros(size,size);
    L(1,1) =sqrt(A(1,1));
    L(2:size,1) = A(2:size,1)/L(1,1);
    
    for j=2:size
        s1=0;
        for l = 1:j-1
            s1 = s1 + L(j,l)*L(j,l);
        end
        L(j,j)=sqrt((A(j,j)- s1));
        for i=j+1:size
            s1=0;
            for l = 1:j-1
                s1 = s1 + L(i,l)*L(j,l);
            end
            L(i,j)=(A(i,j)- s1)/L(j,j);
        end
    end
    Y = zeros(size,1);
    Y(1,1) = A(1,size+1)/L(1,1);
    for i=2:1:size
        for j=1:1:i-1
            Y(i,1) = A(i,size+1)-Y(j,1)*L(i,j);
        end
        Y(i,1) = Y(i,1)/L(i,i);
    end
    X = zeros(size,1);
    X(size,1) = Y(size,1)/L(size,size)';
    for i=size-1:-1:1
        for j=size:-1:i+1
            X(i,1) = Y(i,1)-X(j,1)*L(j,i)';
        end
        X(i,1) = X(i,1)/L(i,i)';
    end
    L;
    L';
    X;
    filename = "output1.txt";
    outf = fopen(filename, "w");
    
    fprintf("Cholesky factor LC :\n");
    disp(L);
    fprintf("X:\n");
    disp(X);
    fclose(outf);
end