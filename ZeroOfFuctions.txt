#include<iostream>
#include<stdexcept>
#include<cmath>
#include <limits>
#include<functional>
#include<vector>
#include<complex>
#include<algorithm>

template <typename FunType>
bool BracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5,
double hmax = 1e10, double lambda = 1.4){
    if(hinit <=0 || hmax<=0 || lambda<=0 ) throw std::domain_error("Invalid parameters");
    double temp_a=x0;
    double f1=f(temp_a);
    double h=hinit;
    while(std::fabs(h)<hmax){
        double temp_b=temp_a+h;
        double f2=f(temp_b);
        while(!std::isfinite(f2)){
            h=h/(2*(1+lambda));
            if(std::fabs(h)<=std::fabs(temp_a)*std::numeric_limits<double>::epsilon()){
                return false;
            }
            temp_b=temp_a+h;
            f2=f(temp_b);
        }
        if(f1*f2<=0){
            if(temp_a>temp_b) std::swap(temp_a, temp_b);
            a=temp_a;
            b=temp_b;
            return true;
        }
        h=h*lambda;
        a=b;
        f1=f2;
    }

    temp_a=x0;
    f1=f(temp_a);
    h=hinit;
    while(std::fabs(h)<hmax){
        double temp_b=temp_a-h;
        double f2=f(temp_b);
        while(!std::isfinite(f2)){
            h=h/(2*(1+lambda));
            if(std::fabs(h)<=std::fabs(temp_a)*std::numeric_limits<double>::epsilon()){
                return false;
            }
            temp_b=temp_a+h;
            f2=f(temp_b);
        }
        if(f1*f2<=0){
            if(temp_a>temp_b) std::swap(temp_a, temp_b);
            a=temp_a;
            b=temp_b;
            return true;
        }
        h=h*lambda;
        a=b;
        f1=f2;
    }
    return false;
}

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b,
RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100){
    if(eps<=0 || maxiter<=0 ) throw std::domain_error("Invalid parameters");
    if(f(a)*f(b)>0) throw std::range_error("Root must be bracketed");

    std::function<double(double)> phi= [](double x){return x/(1+std::fabs(x));};

    double f1,f2,c,c_old;

    if(mode==Unmodified){
        f1=f(a);
        f2=f(b);
        c=a;
        c_old=b;
        int brojac=0;
      while(std::fabs(c-c_old)>eps){
          if(brojac-maxiter==0) throw std::logic_error("Given accuracy has not achieved");
            c_old=c;
            c=(a*f2-b*f1)/(f2-f1);
            double f3=f(c);
            if(std::fabs(f3)<eps){
                return c;
            }
            if(f1*f3<0){
                b=a;
                f2=f1;
            }
            a=c;
            f1=f3;
            brojac++;
        }
        return c; 
    }
    else if(mode==Illinois){
        f1=f(a);
        f2=f(b);
        c=a;
        c_old=b;
        int brojac=0;
         while(std::fabs(c-c_old)>eps){
          if(brojac-maxiter==0) throw std::logic_error("Given accuracy has not achieved");
            c_old=c;
            c=(a*f2-b*f1)/(f2-f1);
            double f3=f(c);
            if(std::fabs(f3)<eps){
                return c;
            }
            if(f1*f3<0){
                b=a;
                f2=f1;
            }
            else{
                f2=f2/2;
            }
            a=c;
            f1=f3;
            brojac++;
        }
        return c; 
    }
    else if(mode==Slavic){
        f1=phi(f(a));
        f2=phi(f(b));
        c=a;
        c_old=b;
        int brojac=0;
      while(std::fabs(c-c_old)>eps){
          if(brojac-maxiter==0) throw std::logic_error("Given accuracy has not achieved");
            c_old=c;
            c=(a*f2-b*f1)/(f2-f1);
            double f3=phi(f(c));
            if(std::fabs(f3)<eps){
                return c;
            }
            if(f1*f3<0){
                b=c;
                f2=f3;
            }
            a=c;
            f1=f3;
            brojac++;
        }
        return c; 
    }
    else if(mode==IllinoisSlavic){
        f1=phi(f(a));
        f2=phi(f(b));
        c=a;
        c_old=b;
        int brojac=0;
      while(std::fabs(c-c_old)>eps){
          if(brojac-maxiter==0) throw std::logic_error("Given accuracy has not achieved");
            c_old=c;
            c=(a*f2-b*f1)/(f2-f1);
            double f3=phi(f(c));
            if(std::fabs(f3)<eps){
                return c;
            }
            if(f1*f3<0){
                b=a;
                f2=f1;
            }
            else{
                f2=f2/2;
            }
            a=c;
            f1=f3;
            brojac++;
        }
        return c; 
    }else{
        return 0;
    }

} 

int sgn(double val) {
    return (0.0 < val) - (val < 0.0);
}

template <typename FunType>
double RiddersSolve(FunType f, double a, double b, double eps = 1e-10,
int maxiter = 100){
    if(eps<=0 || maxiter<=0 ) throw std::domain_error("Invalid parameters");
    if(f(a)*f(b)>0) throw std::range_error("Root must be bracketed");
    double f1=f(a);
    double f2=f(b);
    int brojac=0;
    while(std::fabs(b-a)>=eps){
        if(brojac-maxiter==0) throw std::logic_error("Given accuracy has not achieved");
        double c=(a+b)/2;
        double f3=f(c);
        if(std::fabs(f3)<eps){
            return c;
        }
        double d=c+f3*(c-a)*sgn(f1-f2)/std::sqrt(f3*f3-f1*f2);
        double f4=f(d);
        if(std::fabs(f4)<eps){
            return d;
        }
        if(f3*f4<=0){
            a=c;
            b=d;
            f1=f3;
            f2=f4;
        }
        else if(f1*f4<=0){
            b=d;
            f2=f4;
        }
        else{
            a=d;
            f1=f4;
        }
        brojac++;
    }
    return (a+b)/2;
}

template <typename FunTip1, typename FunTip2>
double NewtonRaphsonSolve(FunTip1 f, FunTip2 fprim, double x0, double eps = 1e-10,double damping = 0, int maxiter = 100)
{
    
    if(eps < 0 || maxiter < 0) throw std::domain_error("Invalid parameters");
    else if( damping < 0 || damping >= 1) throw std::domain_error("Invalid parameters");
    double delta_x=std::numeric_limits<double>::infinity();
    double v=f(x0);
    double d=fprim(x0);
    int brojac=0;

    while(std::fabs(delta_x) > eps) {

        if(std::abs(fprim(x0)) <eps || brojac == maxiter || !std::isfinite(x0)) throw std::logic_error("Convergence has not achieved");
        if(std::fabs(v) <= eps) return x0;
        delta_x = v / d;
        double w= v;
        v = f(x0 - delta_x);
        d = fprim(x0 -delta_x);
        while(std::abs(v) > std::abs(w) || !std::isfinite(v) || d == 0)
        {
            delta_x *= damping;
            v = f(x0 -delta_x);
            d = fprim(x0 - delta_x);
            if (std::abs(delta_x) < eps) {
                throw std::logic_error("Convergence has not achieved");
            }
        }
        x0 -= delta_x;
        brojac++;
    }
    return x0;
}


bool operator==(std::complex<double>c1, std::complex<double> c2) {
    return(std::abs(c1.real() - c2.real()) < std::numeric_limits<double>::epsilon() 
    && std::abs(c1.imag() - c2.imag()) < std::numeric_limits<double>::epsilon());
}

std::complex<double>operator*(std::complex<double>c1, std::complex<double>c2) {
    return {c1.real()*c2.real()-c1.imag()*c2.imag(), c1.real()*c2.imag() + c1.imag()*c2.real()};
}

std::complex<double>operator *(double x, std::complex<double>c) {
    std::complex<double> c2(x,0);
    return c2 * c;
}

std::complex<double>operator *(std::complex<double>c,double x) {
    return x * c;
}

std::pair<std::complex<double>,bool> Laguerre(std::vector<double>p, int n, std::complex<double>x, double eps, int maxiter){
    std::complex<double> delta_x=std::numeric_limits<double>::infinity();
    int k=1;
    while(std::fabs(delta_x)>eps && (k<maxiter)){
        std::complex<double> f=p[n];
        std::complex<double> d=0;
        std::complex<double> s=0;
        for(int i=n-1; i>=0; i--){
            s=s*x+2*d;
            d=d*x+f;
            f=f*x+p[i];
        }
        if(std::fabs(f)<=eps){
            return {x, true};
        }
        std::complex<double> r = std::sqrt((n-1)*((n-1)*d*d-n*f*s));
        if(std::abs(d+r)>std::abs(d-r)){
            delta_x=n*f/(d+r);
        }
        else{
            delta_x=n*f/(d-r);
        } 
        x -= delta_x;
        k++;
    }
    if(std::fabs(delta_x)<=eps){
        return {x, true};
    }
    return {x, false};
}

std::pair<std::complex<double>,bool> Laguerre(std::vector<std::complex<double>>p, int n, std::complex<double>x, double eps, int maxiter){
    std::complex<double> delta_x=std::numeric_limits<double>::infinity();
    int k=1;
    while(std::fabs(delta_x)>eps && (k<maxiter)){
        std::complex<double> f=p[n];
        std::complex<double> d=0;
        std::complex<double> s=0;
        for(int i=n-1; i>=0; i--){
            s=s*x+2*d;
            d=d*x+f;
            f=f*x+p[i];
        }
        if(std::fabs(f)<=eps){
            return {x, true};
        }
        std::complex<double> r = std::sqrt((n-1)*((n-1)*d*d-n*f*s));
        if(std::abs(d+r)>std::abs(d-r)){
            delta_x=n*f/(d+r);
        }
        else{
            delta_x=n*f/(d-r);
        } 
        x -= delta_x;
        k++;
    }
    if(std::fabs(delta_x)<=eps){
        return {x, true};
    }
    return {x, false};
}
std::complex<double> RandomComplex(double rmin, double rmax, double imin, double imax) {
    double realPart = rmin + ((double)rand() / RAND_MAX) * (rmax - rmin);
    double imagPart = imin + ((double)rand() / RAND_MAX) * (imax - imin);
    return std::complex<double>(realPart, imagPart);
}
std::vector<std::complex<double>>PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10,
int maxiters = 100, int maxtrials = 10){
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    int n=coefficients.size()-1;
    int i=n;
    int brojac=0;
    std::vector<std::complex<double>> v_0;
    while(i>0){
        if(brojac==maxiters) throw std::logic_error("Convergence has not achieved");
        int t=1;
        bool c=false;
        std::complex<double> x;
        while(!c && (t<maxtrials)){
            x=RandomComplex(-10,10,-10,10);
            std::pair<std::complex<double>, bool> temp=Laguerre(coefficients, i, x, eps, maxiters);
            x=temp.first;
            c=temp.second;
            t++;
        }
        if(!c){
            throw std::logic_error("Convergence has not achieved");
        }
        if(std::fabs(x.imag())<=eps){
            x=x.real();
        }
        v_0.push_back(x);

        std::complex<double> v=coefficients[i];
        for(int j=i-1; j>=0; j--){
            std::complex<double> w=coefficients[j];
            coefficients[j]=v;
            v=w+x*v;
        }
        brojac++;
        i--;
    }
    std::sort(v_0.begin(),v_0.end(), [] (std::complex<double> x, std::complex<double> y) {
        if (x.real() < y.real()) return true;
        else if (x.real() > y.real()) return false;
        return x.imag() < y.imag();
    });

    return v_0;
}

std::vector<std::complex<double>> PolyRoots(std::vector<double> coefficients, double eps = 1e-10,
int maxiters = 100, int maxtrials = 10){
    if(eps < 0 || maxiters < 0 || maxtrials < 0) throw std::domain_error("Invalid parameters");
    
    int n=coefficients.size()-1;
    int i=n;
    int brojac=0;
    std::vector<std::complex<double>> v_0(n+1);
    while(i>=1){
        if(brojac==maxiters) throw std::logic_error("Convergence has not achieved");
        int t=1;
        bool c=false;
        std::complex<double>x;
        while(!c && (t<maxtrials)){
            x=RandomComplex(-10,10,-10,10);
            std::pair<std::complex<double>,bool> temp=Laguerre(coefficients, i, x, eps, maxiters);
            x=temp.first;
            c=temp.second;
            t++;
        }
        if(!c){
            throw std::logic_error("Convergence has not achieved");
        }
        if(std::fabs(x.imag())<=eps){
            x=x.real();
            v_0[i]=x;
            double v=coefficients[i];
            for(int j = i - 1; j >= 0; j--) {
                double w = coefficients[j];
                coefficients[j] = v;
                v = w + x.real() * v;
            }
            i--;
        }
        else{
            v_0[i]=x;
            v_0[i-1]=std::conj(x);
            double alpha=2*x.real();
            double beta=std::fabs(x)*std::fabs(x);
            double u=coefficients[i];
            double v=coefficients[i-1]+alpha*u;
            for(int j=i-2;j>=0; j--) {
                double w=coefficients[j];
                coefficients[j]=u;
                u=v;
                v=w+alpha*v-beta*coefficients[j];
            }
            i-=2;
        }
        brojac++;
    }
    v_0.erase(v_0.begin());
    std::sort(v_0.begin(),v_0.end(),[](std::complex<double>c1, std::complex<double>c2) {
        if(std::abs(c1.real() - c2.real()) > std::numeric_limits<double>::epsilon()) return c1.real() < c2.real();
        return c1.imag() < c2.imag();
    });
    return v_0;
}


double testFunction1(double x) {
    return x * x - 4;
}

double testFunction2(double x) {
    return x * x * x - 6 * x * x + 11 * x - 6;
}

double testFunction3(double x) {
    return sin(x);
}

std::complex<double> testComplexFunction(std::complex<double> z) {
    return z * z - 4.0;
}

int main() {
    try {
        std::cout << "BracketRoot Test:" << std::endl;
        double a, b;

        try {
            bool bracketSuccess1 = BracketRoot(testFunction1, 0, a, b);
            if (bracketSuccess1)
                std::cout << "Bracket for testFunction1: [" << a << ", " << b << "]" << std::endl;
            else
                std::cout << "Bracket for testFunction1 not found." << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in BracketRoot for testFunction1: " << e.what() << std::endl;
        }

        try {
            bool bracketSuccess2 = BracketRoot(testFunction2, 0, a, b);
            if (bracketSuccess2)
                std::cout << "Bracket for testFunction2: [" << a << ", " << b << "]" << std::endl;
            else
                std::cout << "Bracket for testFunction2 not found." << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in BracketRoot for testFunction2: " << e.what() << std::endl;
        }

        try {
            bool bracketSuccess3 = BracketRoot(testFunction3, 0, a, b);
            if (bracketSuccess3)
                std::cout << "Bracket for testFunction3: [" << a << ", " << b << "]" << std::endl;
            else
                std::cout << "Bracket for testFunction3 not found." << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in BracketRoot for testFunction3: " << e.what() << std::endl;
        }

        std::cout << "\nRegulaFalsiSolve Test:" << std::endl;

        try {
            auto quadraticFunction = [](double x) { return x * x - 4; };
            double a_boundary = 1.0;
            double b_boundary = 4.5;
            double root = RegulaFalsiSolve(quadraticFunction, a_boundary, b_boundary, Unmodified);
            std::cout << "Root of the quadratic function: " << root << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in RegulaFalsiSolve: " << e.what() << std::endl;
        }
        std::cout << "\nRegulaFalsiSolve Non-Converging Example:" << std::endl;

        auto nonConvergingFunction = [](double x) { return 1.0 / x; };
        double a_boundary_non_converging = -1.5;
        double b_boundary_non_converging = 1.5;

        try {
            double root_non_converging = RegulaFalsiSolve(nonConvergingFunction, a_boundary_non_converging, b_boundary_non_converging, Unmodified);
            std::cout << "Root of the non-converging function: " << root_non_converging << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in RegulaFalsiSolve Non-Converging Example: " << e.what() << std::endl;
        }
        
        std::cout << "\nNewtonRaphsonSolve Test:" << std::endl;

        try {
            double newtonRoot1 = NewtonRaphsonSolve(testFunction1, [](double x) { return 2 * x; }, 1);
            double newtonRoot2 = NewtonRaphsonSolve(testFunction2, [](double x) { return 3 * x * x - 12 * x + 11; }, 1);
            double newtonRoot3 = NewtonRaphsonSolve(testFunction3, [](double x) { return cos(x); }, 1);

            std::cout << "Root for testFunction1 using Newton's method: " << newtonRoot1 << std::endl;
            std::cout << "Root for testFunction2 using Newton's method: " << newtonRoot2 << std::endl;
            std::cout << "Root for testFunction3 using Newton's method: " << newtonRoot3 << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error in NewtonRaphsonSolve: " << e.what() << std::endl;
        }

        
        std::cout << "\nPolyRoots Test:" << std::endl;

        try {
            std::vector<double> polyCoeff1 = {1, 0, -4};
            std::vector<double> polyCoeff2 = {1, -6, 11, -6};
            

            auto roots1 = PolyRoots(polyCoeff1);
            auto roots2 = PolyRoots(polyCoeff2);
            

            std::cout << "Roots for polyCoeff1: ";
            for (auto root : roots1)
                std::cout << root << " ";
            std::cout << std::endl;

            std::cout << "Roots for polyCoeff2: ";
            for (auto root : roots2)
                std::cout << root << " ";
            std::cout << std::endl;

            
        } catch (const std::exception &e) {
            std::cerr << "Error in PolyRoots with real coefficients: " << e.what() << std::endl;
        }

        
        std::cout << "\nPolyRoots Test with Complex Coefficients:" << std::endl;

        try {
            std::vector<std::complex<double>> complexCoeff1 = {1, 0, -4};
            std::vector<std::complex<double>> complexCoeff2 = {1, -6, 11, -6};
            

            auto complexRoots1 = PolyRoots(complexCoeff1);
            auto complexRoots2 = PolyRoots(complexCoeff2);
            

            std::cout << "Roots for complexCoeff1: ";
            for (auto root : complexRoots1)
                std::cout << root << " ";
            std::cout << std::endl;

            std::cout << "Roots for complexCoeff2: ";
            for (auto root : complexRoots2)
                std::cout << root << " ";
            std::cout << std::endl;

            
        } catch (const std::exception &e) {
            std::cerr << "Error in PolyRoots with complex coefficients: " << e.what() << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
