choice = input('Choose the method of solution by selecting number shown with method: Muller-1, Bairstow-2\n','s');
if(choice=='1')
    degree = input('Enter the degree of Polynomial f(x):\n');
    coef_of_polynomial = input('Enter (degree+1) number of coefficients in increasing order of degree separated by a space:\n','s');
    coef_1 = str2num(coef_of_polynomial);
    
    i=0;
    str=num2str(coef_1(1));
    while(i<degree)
        i=i+1;
        str = strcat(str,'+(',num2str(coef_1(i+1)),')*x^',num2str(i));
    end
    
    str = strcat('@(x)',str);
    fun = str2func(str);
    i=0;
    am = input('Enter the first starting point\n');
    bm = input('Enter the second starting point\n');
    cm = input('Enter the third starting point\n');
    
    disp('Enter the stopping criteria:');
    max_iter = input('Enter the maximum number of iterations allowed\n');
    rel_err = input('Enter the maximum relative approximate error allowed(in %)\n');
    
    x1=am;
    x2=bm;
    x3=cm;
    err=100;
    iter=0;
    stop_err = 0;
    stop_iter=0;
    stop_value=0;
    i=1;
    
    while((err>rel_err) && (iter<max_iter) )
        c = fun(x3);
        a = ((fun(x3)-fun(x2))/(x3-x2)-(fun(x3)-fun(x1))/(x3-x1))/(x2-x1);
        b = (fun(x3)-fun(x1))/(x3-x1)-a*(x1-x3);
        x1=x2;
        x2=x3;
        if(b>0)
            x3=x3-(2*c/(b+sqrt(b*b-4*a*c)));
        else
            x3=x3-(2*c/(b-sqrt(b*b-4*a*c)));
            p=x3;
        end
        err=abs((x3-x2)/x3);
        errg(i)=err;
        iter=iter+1;
        i=i+1;
        
        if(err<rel_err)
            stop_err=1;
        end
        if(iter==max_iter)
            stop_iter=1;
        end
    end
    answer = x3;
    fprintf('Root is %f+%fi\n',real(p),imag(p));
    
    if(stop_err)
        disp('Iterations stopped: Maximum relative error stopping criteria met');
    end
    if(stop_iter)
        disp('Iterations stopped: Maximum number of iterations reached');
    end
    
    figure(2)
    ezplot(fun);
    grid on;
end

if(choice=='2')
    func = input('Enter the polynomial f(x): ','s');
    f = str2sym(func);
    r = input('Enter starting value of r: ');
    s = input('Enter starting value of s: ');
    relApproxErr = input('Enter the maximum relative approximate error allowed(in %): ');
    max_iter = input('Enter the maximum number of iterations allowed: ');
    y = matlabFunction(f);
    n = polynomialDegree(f) + 1;
    m = n - 1;
    a = sym2poly(f);
    a = fliplr(a);
    M = zeros;
    N = zeros;
    while n > 0
        for j = 1: max_iter
            b = zeros;
            c = zeros;
            b(n) = a(n);
            b(n-1) = a(n-1) + r*b(n);
            c(n) = b(n);
            c(n-1) = b(n-1) + r*c(n);
            
            for i = n-2: -1: 1
                b(i) = a(i) + r*b(i+1) + s*b(i+2);
                c(i) = b(i) + r*c(i+1) + s*c(i+2);
            end
            
            syms x y
            eqn1 = c(2)*x + c(3)*y == -b(1);
            eqn2 = c(3)*x + c(4)*y == -b(2);
            sol = solve([eqn1, eqn2], [x, y]);
            delr = sol.x;
            dels = sol.y;
            errPercentr = abs((delr/(r+delr))*100);
            errPercents = abs((dels/(r+dels))*100);
            N(j) = errPercents;
            M(j) = errPercentr;
            if errPercents > relApproxErr || errPercentr > relApproxErr
                r = r + delr;
                s = s + dels;
            end
            if errPercents < relApproxErr && errPercentr < relApproxErr
                r = r + delr;
                s = s + dels;
                break;
            end
        end
        m = m - 2;
        if (r^2 + 4*s) < 0
            x1 = r/2;
            x2 = (-(r^2 + 4*s))/2;
            fprintf('Roots of the function are: %f + %fi, %f - %fi\n', x1, x2, x1, x2);
        else
            x1 = (r + (r^2 + 4*s)^(1/2))/2;
            x2 = (r - (r^2 + 4*s)^(1/2))/2;
            fprintf('Roots of the function are: %f, %f\n', x1, x2);
        end
        if m > 2
            a = ones;
            for k = n : -1: 3
                a(k-2) = b(k);
            end
        end
        n = n - 2;
        if m == 2
            if (b(4)^2 - 4*b(5)*b(3)) < 0
                x3 = -b(4)/2*b(5);
                x4 = (-(b(4)^2 - 4*b(5)*b(3)))^(1/2)/2*b(5);
                fprintf('Roots of the function are: %f + %fi, %f - %fi\n', x3, x4, x3, x4);
            else
                x3 = (-b(4) + (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                x4 = (-b(4) - (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                fprintf('Roots of the function are: %f, %f\n', x3, x4);
            end
            break;
        end
        if m == 1
            x1 = -b(3)/b(4);
            fprintf('Roots of the function are: %f\n', x1);
            break;
        end
    end
    figure(1);
    ezplot(f)
    grid on   
end