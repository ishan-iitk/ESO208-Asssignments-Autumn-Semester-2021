choice = input('Choose the method of solution by selecting number shown with method:\n Bisection-1,\n False Position-2,\n Modified False Position-3,\n Newton-Raphson-4,\n Secant-5 \n','s');
if(choice == '1')
    clear
    format long
    b = input('Enter the function f(x)\n','s');
    str_1 = strcat('@(x)', b);
    f_bis = str2func(str_1);
    xl = input('Enter the first starting point\n');
    xu = input('Enter the second starting point\n');
    disp('Enter stopping criteria:');
    relative_error = input('Enter the maximum relative error allowed (in %)\n');
    maximum_iterations = input('Enter the maximum number of iterations  allowed\n');
    
    c1 = 0;
    iter = 0;
    xr = (xl+xu)/2;
    error = abs(100*(xr-c1)/xr);
    i=1;
    
    while ((error > relative_error) && (iter<maximum_iterations))
        if(f_bis(xl)*f_bis(xu)<0)
            xr = (xl+xu)/2;
            if(f_bis(xl)*f_bis(xr)<0)
                xu = xr;
            else
                xl = xr;
            end
        end
        error = abs(100*(xr-c1)/xr);
        c1 = xr;
        p=xr;
        errg(i)=error;
        iter = iter + 1;
        i=i+1;    
    end
    fprintf('The root is: %f\n',p);
    figure(1)
    plot(errg);
    title('Plot of Relative Error vs Iterations');
    xlabel('Number of Iterations');
    ylabel('Relative Error');
    figure(2)
    ezplot(f_bis);
    grid on;
end

if(choice == '2')
    clear
    fun_2 = input('Enter the function f(x):\n','s');
    string_fun= strcat('@(x)',fun_2);
    f_fal = str2func(string_fun);
    xl = input('Enter the first starting point\n');
    xu = input('Enter the second starting point\n');
    disp('Enter stopping criteria:');
    relative_error = input('Enter the maximum relative error allowed(in %)\n');
    maximum_iterations = input('Enter the maximum number of iterations  allowed\n');

    x0=xl;
    x1=xu;
    error = 100;
    x2 = x1 - (f_fal(x1)*(x1 - x0))/(f_fal(x1) - f_fal(x0));
    i=1;
    if(f_fal(xl)*f_fal(xu)>0)
        disp('False position method cannot be used to solve this equation');
    else
    
    while ((error > relative_error) && (i<=maximum_iterations))
        
        if (f_fal(x2)*f_fal(x0)<0)
            x1=x2;
        else
            x0=x2;  
        end
        
        p=x2;
        x2 = x1 - f_fal(x1)*(x1 - x0)/(f_fal(x1) - f_fal(x0));
        
        error = abs(100*(x2-x1)/x2);
        errg(i) = error;
        
        i=i+1;
       
    end
    
    fprintf('The root is: %f\n',p);
    figure(1)
    ezplot(f_fal);
    figure(2)
    plot(errg);
    title('Plot of Relative Error vs Iterations');
    xlabel('Number of Iterations');
    ylabel('Relative Error');
    end
end

if(choice == '3')
    clear
    func_3 = input('Enter the function f(x)\n','s');
    string_1 = strcat('@(x)',func_3);
    f_modfal = str2func(string_1);
    xl = input('Enter the first starting point\n');
    xu = input('Enter the second starting point\n');
    disp('Enter stopping criteria:');
    relative_error = input('Enter the maximum relative error allowed (in %)\n');
    maximum_iterations = input('Enter the maximum number of iterations allowed\n');
    
    iter = 0;
    error = 100;
    xr = xu - f_modfal(xu)*(xu - xl)/(f_modfal(xu) - f_modfal(xl));
    i=1;
    fl=f_modfal(xl);
    fu=f_modfal(xu);
    il=0;
    iu=0;
    
    while ((error > relative_error) && (iter<maximum_iterations))
        xrold=xr;
        xr = xu - f_modfal(xu)*(xu - xl)/(f_modfal(xu) - f_modfal(xl));
        fr=f_modfal(xr);
        iter=iter+1;
        test=fl*fr;
        if (test<0)
            xu=xr;
            fu=f_modfal(xu);
            iu=0;
            il=il+1;
            if(il>=2)
                fl=fl/2;
            end
        else
            if(test>0)
                xl=xr;
                fl=f_modfal(xl);
                il=0;
                iu=iu+1;
                if(iu>=2)
                    fu=fu/2;
                end
            end
        end
        
        p=xr;
        error = abs(100*(xr-xu)/xr);
        errg(i) = error;
        i=i+1;    
    end
    
    fprintf('The root is: %f\n',p);
    figure(1)
    ezplot(f_modfal);
    figure(2)
    plot(errg);
    title('Plot of Relative Error vs Iterations');
    xlabel('Number of Iterations');
    ylabel('Relative Error');
end

if(choice == '4')
    clear
    format longg
    func_ask4 = input('Enter the function f(x)\n','s');
    str1 = strcat('@(x)',func_ask4);
    func_ask41 = input('Enter the first derivative of function f(x)\n','s');
    str2 = strcat('@(x)',func_ask41);
    f_newton1 = str2func(str1);
    f_newton2 = str2func(str2);
    a_newton = input('Enter the starting point\n');
    disp('Enter stopping criteria:');
    relative_error = input('Enter the maximum relative error allowed (in %)\n');
    maximum_iterations = input('Enter the maximum number of iterations  allowed\n');
    
    x1=a_newton;
    x2 = x1 - f_newton1(x1)/f_newton2(x1);
    err = 100;
    iter = 0;
    i=1;
    
    while ((err > relative_error) && (iter<maximum_iterations))
        x2 = x1 - f_newton1(x1)/f_newton2(x1);
        err = abs(100*(x2-x1)/x2);
        x1 = x2;
        p=x2;
        errg(i) = err;
        iter = iter+1;
        i=i+1;
        error(iter)=err;
    end
    
    fprintf('Root is: %f\n',p); 
    figure(1)
    ezplot(f_newton1)
    figure(2)
    plot(errg);
    title('Plot of Relative Error vs Iterations');
    xlabel('Number of Iterations');
    ylabel('Relative Error');
end

if(choice == '5')
    format longg
    func_ask5 = input('Enter the function f(x)\n','s');
    str1 = strcat('@(x)',func_ask5);
    f_secant = str2func(str1);
    a_secant = input('Enter the first starting point\n');
    b_secant = input('Enter the second starting point\n');
    disp('Enter stopping criteria:\n');
    relative_error = input('Enter the maximum relative error allowed (in %)\n');
    
    maximum_iterations = input('Enter the maximum number of iterations  allowed\n');
    x0 = a_secant;
    x1 = b_secant;
    x2 = x1 - f_secant(x1)*(x1 - x0)/(f_secant(x1) - f_secant(x0));
    iter=0;
    err = 100;
    stop_val=0;
    i=1;
    
    while ((err > relative_error) && (iter<=maximum_iterations))
        x2 = x1 - f_secant(x1)*(x1 - x0)/(f_secant(x1) - f_secant(x0));
        err = abs(100*(x2-x1)/x2);
        x0 = x1;
        x1 = x2;
        errg(i)=err;
        p=x2;
        i=i+1;
        iter = iter + 1;
    end
    
    fprintf('Root is: %f\n',p);
    figure(1)
    ezplot(f_secant)
    figure(2)
    plot(errg);
    title('Plot of Relative Error vs Iterations');
    xlabel('Number of Iterations');
    ylabel('Relative Error');
end