#include<iostream>
#include<stdexcept>
#include<cmath>
#include<vector>

template <typename FunType>
std::pair<double, bool> RombergIntegration(FunType f, double a, double b,double eps = 1e-8, 
    int nmax = 1000000, int nmin = 50){
        if(eps<0 || nmin<0 || nmax<0 ||nmax<nmin) throw std::domain_error("Bad parameter");
        int znak=1;
        if(a>b){
            std::swap(a,b);
            znak=-1;
        }
        int N=2;
        std::vector<double> I;
        double h=(b-a)/N;
        double s=(f(a)+f(b))/2;
        double I_old=s;
        for(int i=1; N<=nmax; i++){
            for(int j=1; j<=N/2; j++){
                s+=f(a+(2*j-1)*h);
            }
            I.push_back(h*s);
            double p=4;
            for(int k=I.size()-2; k>=0; k--){
                I[k]=(p*I[k+1]-I[k])/(p-1);
                p*=4;
            }
            if(std::fabs(I[0]-I_old)<=eps && N>=nmin) return std::make_pair(znak*I[0], true);
            I_old=I[0];
            h/=2;
            N*=2;
        }
        return std::make_pair(I_old*znak, false);
}

template <typename FunType>
std::pair<double, bool> TanhSinhIntegration(FunType f, double a, double b,double eps = 1e-8, 
    int nmax = 1000000, int nmin = 20, double range = 3.5){
        if(eps<0 || nmin<0 || nmax<0 ||nmax<nmin || range<0) throw std::domain_error("Bad parameter");
        int znak=1;
        if(a>b){
            std::swap(a,b);
            znak=-1;
        }
        int N=2;
        double h=2*range/N;
        double p=(b+a)/2;
        double q=(b-a)/2;
        double s=0;
        double I_old=s;
        while(N<nmax){
            for(int i=1; i<=N/2; i++){
                double t=-range+(2*i-1)*h;
                double u=M_PI*std::sinh(t)/2;
                double v=f(p+q*std::tanh(u));
                if(std::isfinite(v)){
                    s=s+q*M_PI*std::cosh(t)*v/(2*(std::cosh(u)*std::cosh(u)));
                }
            }
            double I=h*s;
            if(N>=nmin && std::abs(I-I_old)<=eps){
                return std::make_pair(I*znak, true);
            }
            I_old=I;
            N=2*N;
            h=h/2;
        }
        return std::make_pair(I_old*znak, false);
}

std::pair<double,bool>operator+ (std::pair<double,bool>a, std::pair<double,bool>b) {
    return {a.first+b.first, a.second && b.second};
}

template<typename FunType>
std::pair<double,bool> AdaptiveAux(FunType f, double a, double b, double eps, double f1, double f2, double f3, int R){
    if(!std::isfinite(f1)) f1 = 0;
    if(!std::isfinite(f2)) f2 = 0;
    if(!std::isfinite(f3)) f3 = 0;
    double c=(a+b)/2;
    double I_1=(b-a)*(f1+4*f3+f2)/6;
    double f4=f((a+c)/2);
    double f5=f((c+b)/2);
    if(!std::isfinite(f4)) f4 = 0;
    if(!std::isfinite(f5)) f5 = 0;
    double I_2=(b-a)*(f1+4*f4+2*f3+4*f5+f2)/12;
    if(std::fabs(I_1-I_2)<=eps) return {I_2, true};
    if(R<=0) return {I_2, false};
    return AdaptiveAux(f,a, c, eps, f1, f3, f4, R-1)+AdaptiveAux(f, c, b, eps, f3, f2, f5, R-1);
}

template <typename FunType>
std::pair<double, bool> AdaptiveIntegration(FunType f, double a, double b,
    double eps = 1e-10, int maxdepth = 30, int nmin = 1){
        if(eps<0 || maxdepth<0 || nmin<0) throw std::domain_error("Bad parameter");
        int znak=1;
        if(a>b){
            std::swap(a,b);
            znak=-1;
        }
        double h=(b-a)/nmin;
        std::pair<double, bool>s = {0, true};
        for(int i=1; i<=nmin; i++){
            s=s+AdaptiveAux(f, a, a+h, eps, f(a), f(a+h), f(a+h/2), maxdepth);
            a+=h;
        }
        s.first=s.first*znak;
        return s;
}
int main ()
{
    
    auto sin_function = [](double x) { return std::sin(x); };
    auto exp_function = [](double x) { return std::exp(x); };
    auto quadratic_function = [](double x) { return 3 * x * x; };
    auto inverse_sqrt_function = [](double x) { if(x==0) return 0.0; else return 1/std::sqrt(x); };

    try{
        
    std::cout << "Romberg Integration:" << std::endl;
    auto result_sin_romberg = RombergIntegration(sin_function, 0, M_PI);
    std::cout << "sin(x) Result: " << result_sin_romberg.first << ", Success: " << result_sin_romberg.second << std::endl;

    auto result_exp_romberg = RombergIntegration(exp_function, 0, 1);
    std::cout << "exp(x) Result: " << result_exp_romberg.first << ", Success: " << result_exp_romberg.second << std::endl;

    auto result_quadratic_romberg = RombergIntegration(quadratic_function, 1, 3);
    std::cout << "3 * x^2 Result: " << result_quadratic_romberg.first << ", Success: " << result_quadratic_romberg.second << std::endl;

    
    std::cout << "\nAdaptive Integration:" << std::endl;
    auto result_sin_adaptive = AdaptiveIntegration(sin_function, 0, M_PI);
    std::cout << "sin(x) Result: " << result_sin_adaptive.first << ", Success: " << result_sin_adaptive.second << std::endl;

    auto result_exp_adaptive = AdaptiveIntegration(exp_function, 0, 1);
    std::cout << "exp(x) Result: " << result_exp_adaptive.first << ", Success: " << result_exp_adaptive.second << std::endl;

    auto result_quadratic_adaptive = AdaptiveIntegration(quadratic_function, 1, 3);
    std::cout << "3 * x^2 Result: " << result_quadratic_adaptive.first << ", Success: " << result_quadratic_adaptive.second << std::endl;

    
    std::cout << "\nTanhSinh Integration:" << std::endl;
    auto result_sin_tanhsinh = TanhSinhIntegration(sin_function, 0, M_PI / 2);
    std::cout << "sin(x) Result: " << result_sin_tanhsinh.first << ", Success: " << result_sin_tanhsinh.second << std::endl;

    auto result_exp_tanhsinh = TanhSinhIntegration(exp_function, 0, 1);
    std::cout << "exp(x) Result: " << result_exp_tanhsinh.first << ", Success: " << result_exp_tanhsinh.second << std::endl;

    auto result_quadratic_tanhsinh = TanhSinhIntegration(quadratic_function, 1, 3);
    std::cout << "3 * x^2 Result: " << result_quadratic_tanhsinh.first << ", Success: " << result_quadratic_tanhsinh.second << std::endl;

    
    std::cout << "\nAdditional Test Case:" << std::endl;
    auto result_inverse_sqrt_adaptive = AdaptiveIntegration(inverse_sqrt_function, 0, 1);
    std::cout << "1 / sqrt(x) Result: " << result_inverse_sqrt_adaptive.first << ", Success: " << result_inverse_sqrt_adaptive.second << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what();
    }


    return 0;
}