#include<iostream>
#include<stdexcept>
#include<cmath>

template <typename FunType>
double FindMinimum(FunType f, double x0, double eps = 1e-8,
double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4){
    if(eps<=0 || hinit<=0 || hmax<=0 || lambda<=0) throw std::domain_error("Invalid parameters");
    double a=x0-hinit;
    double b=x0-hinit;
    double c=x0;
    bool found=false;
    while(std::fabs(hinit)<hmax){
        if(f(c+hinit)<f(c)){
            b=c+hinit;
            a=c-hinit;
        }
        else if(f(c-hinit)<f(c)){
            b = c - hinit;
            a = b - hinit;
        }
        else{
            a=c-hinit;
            b=c+hinit;
            found=true;
            break;
        }
        c=b;
        hinit=hinit*lambda;
    }
    if(!found){
        throw std::logic_error("Minimum has not found");
    }
    double phi=(1+std::sqrt(5))/2;
    c=b-(b-a)/phi;
    double d=a+(b-a)/phi;
    double u=f(c);
    double v=f(d);
    while(std::fabs(b-a)>eps){
        if(u<v){
            b=d;
            d=c;
            c=a+(c-a)/phi;
            v=u;
            u=f(c);
        }
        else{
            a=c;
            c=d;
            d=b-(b-d)/phi;
            u=v;
            v=f(d);
        }
    }
    return (a+b)/2;
}

double quadraticFunction(double x) {
    return x * x;
}

double logMinimumFunction(double x) {
    return 2 * std::log(x + 1) + 3 * x + 5; 
}

double sinFunction(double x) {
    return std::sin(x);
}

double noMinimumFunction(double x) {
    return x;  
}

int main() {
    try {
      
        try {
            double minQuadratic = FindMinimum(quadraticFunction, 1.0);
            std::cout << "Minimum of quadratic function: " << minQuadratic << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in quadratic function: " << e.what() << std::endl;
        }

        
        try {
            double minLog = FindMinimum(logMinimumFunction, -1.0);
            std::cout << "Minimum of log-containing function: " << minLog << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in log-containing function: " << e.what() << std::endl;
        }

        
        try {
            double minSin = FindMinimum(sinFunction, 1.0);
            std::cout << "Minimum of sin function: " << minSin << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in sin function: " << e.what() << std::endl;
        }

        
        try {
            FindMinimum(noMinimumFunction, 1.0);
        } catch (const std::exception &e) {
            std::cerr << "Error in function without minimum: " << e.what() << std::endl;
        }

        
        try {
            FindMinimum(quadraticFunction, 1.0, -1.0);
        } catch (const std::exception &e) {
            std::cerr << "Error with invalid parameters: " << e.what() << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}