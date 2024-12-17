#include<iostream>
#include<stdexcept>
#include<vector>
#include<cmath>

using namespace std;


class ChebyshevApproximation{
    
    double xmin;
    double xmax;
    std::vector<double> coef;
    double m;
    //privatni konstrukotr sa vec datim koeficijentima
    ChebyshevApproximation(double xmin, double xmax, std::vector<double> coef){
        if(xmin>=xmax || coef.size()-1<1) throw std::domain_error("Bad parameters");
        
        this->xmin=xmin;
        this->xmax=xmax;
        this->coef=coef;
        this->m=coef.size()-1;
    }
    public:
    template<typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n){
        if(xmin>=xmax || n<1) throw std::domain_error("Bad parameters");
       
        this->xmin=xmin;
        this->xmax=xmax;
        this->m=n;
        coef.resize(n+1);
        std::vector<double> w(n+2);
        std::vector<double> v(n+1);
        for(int i=0; i<=n+1; i++){
            w[i]=std::cos(M_PI*i/(2*n+2));
        }
        for(int i=0; i<=std::floor(n/2); i++){
            v[i]=f((xmin+xmax+(xmax-xmin)*w[2*i+1])/2);
        }
        for(int i=std::floor(n/2)+1; i<=n; i++){
            v[i]=f((xmin+xmax+(xmax-xmin)*w[2*n+1-2*i])/2);
        }
        for(int k=0; k<=n; k++){
            double s=0;
            for(int i=0; i<=n; i++){
                int p=(k*(2*i+1))%(4*n+4);
                if(p>2*n+2){
                    p=4*n+4-p;
                }
                if(p>n+1){
                    s=s-v[i]*w[2*n+2-p];
                }
                else{
                    s=s+v[i]*w[p];
                }
            }
            coef[k]=2*s/(n+1);
        }

    }
    void set_m(int m){
        if(m<=1 || m>coef.size()-1) throw std::domain_error("Bad order");
        this->m=m;
    }
    void trunc(double eps){
        if(eps<0) throw std::domain_error("Bad tolerance");
        for(int i = m; i >= 0; i--)
        if(std::fabs(coef[i]) > eps) {
            this->m = i;
            return;
        }
        
        throw std::domain_error("Bad tolerance");
    }
    double operator()(double x) const{
        if(x<xmin || x>xmax) throw std::domain_error("Bad argument");
        double t=(2*x-xmin-xmax)/(xmax-xmin);
        double p=0;
        double q=coef[this->m];
        for(int k=m-1; k>=1; k--){
            double r=2*t*q-p+coef[k];
            p=q;
            q=r;
        }
        return t*q-p+coef[0]/2;

    }
    double derivative(double x) const{
        if(x<xmin || x>xmax) throw std::domain_error("Bad argument");
        double t=(2*x-xmin-xmax)/(xmax-xmin);
        double p=1;
        double q=4*t;
        double s=coef[1]+q*coef[2];
        for(int k=3; k<=m; k++){
            double r=k*(2*t*q/(k-1)-p/(k-2));
            s=s+coef[k]*r;
            p=q;
            q=r;
        }
        return 2*s/(xmax-xmin);
    }
    ChebyshevApproximation derivative() const{
        std::vector<double> coef_der(coef.size());
        double temp=4/(xmax-xmin);
        coef_der[m-1]=temp*m*coef[m];
        coef_der[m-2]=temp*(m-1)*coef[m-1];
        for(int k=m-3; k>=0; k--){
            coef_der[k]=coef_der[k+2]+temp*(k+1)*coef[k+1];
        }
        ChebyshevApproximation f_der(xmin, xmax, coef_der);
        return f_der;
    }
    ChebyshevApproximation antiderivative() const{
        std::vector<double> coef_antider(m+2);
        coef_antider[0]=0;
        double temp=(xmax-xmin)/4;
        for(int k=1; k<=m-1; k++){
            coef_antider[k]=temp/k*(coef[k-1]-coef[k+1]);
        }
        coef_antider[m]=temp/m*coef[m-1];
        coef_antider[m+1]=temp/(m+1)*coef[m];
        ChebyshevApproximation f_antider(xmin, xmax, coef_antider);
        return f_antider;
    }
    double integrate(double a, double b) const{
        if(a<xmin || b>xmax || a>xmax || b<xmin) throw std::domain_error("Bad interval");
        int znak=1;
        if(a>b){
            std::swap(a,b);
            znak=-1;
        }
        auto F(antiderivative());
        return znak*(F(b)-F(a));
    }
    double integrate() const{
        double temp=0;
        for(int i=1; i<=(m+1)/2; i++){
            temp+=(2*coef[2*i]/(1-4*i*i));
        }
        temp*=(xmax-xmin)/2;
        temp+=coef[0]*(xmax-xmin)/2;
        return temp;
    }
};

int main() {
    try {
        
        auto func = [](double x) { return std::sin(x); };
        double xmin = 0.0;
        double xmax = M_PI;
        int n = 10;
        std::cout<< "Function approximated: sin(x)" << std::endl;
        ChebyshevApproximation cheb1(func, xmin, xmax, n);

        cheb1.set_m(7);

        double eps = 1e-6;
        cheb1.trunc(eps);
        
        double x_val = M_PI / 2;
        double result = cheb1(x_val);
        std::cout << "f(" << x_val << ") = " << result << std::endl;
        
        double derivative_result = cheb1.derivative(x_val);
        std::cout << "f'(" << x_val << ") = " << derivative_result << std::endl;

        ChebyshevApproximation cheb_derivative = cheb1.derivative();
        
        ChebyshevApproximation cheb_antiderivative = cheb1.antiderivative();
        
        double a = 0.0;
        double b = M_PI;
        double integral_result = cheb1.integrate(a, b);
        std::cout << "Integral from " << a << " to " << b << " = " << integral_result << std::endl;

        
        double indefinite_integral_result = cheb1.integrate();
        std::cout << "Indefinite integral = " << indefinite_integral_result << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }

    return 0;
}