#include<iostream>
#include<cmath>
#include <limits>
#include<stdexcept>
#include<vector>
#include<algorithm>

template <typename FunType>
std::pair<double, bool> Limit(FunType f, double x0, double h = 0,double eps = 1e-8, double nmax = 20){
    if(h==0) h=0.001*std::max(1.0, std::abs(x0));
    if(eps<=0 || nmax<3 || nmax>30) throw std::domain_error("Invalid parameters");
    if(x0==std::numeric_limits<double>::infinity() || x0*-1==std::numeric_limits<double>::infinity()){
        double x_0=0;
        double y_old=std::numeric_limits<double>::infinity();
        std::vector<double> y(nmax);
        for(int i=0; i<nmax; i++){
            y[i]=f(1/x_0+h);
            double p=2;
            for(int k=i-1; k>=0; k--){
                y[k]=(p*y[k+1]-y[k])/(p-1);
                p*=2;
            }
            if(std::fabs(y[0]-y_old)<eps){
                if(x0<0) y[0]*=-1;
                return std::pair<double,bool>(y[0], true);
            }
            y_old=y[0];
            h=h/2;
        }
        if(x0<0) y[0]*=-1;
        return std::pair<double, bool>(y[0], false);
    }
    else{
        double y_old=std::numeric_limits<double>::infinity();
        std::vector<double> y(nmax);
        for(int i=0; i<nmax; i++){
            y[i]=f(x0+h);
            double p=2;
            for(int k=i-1; k>=0; k--){
                y[k]=(p*y[k+1]-y[k])/(p-1);
                p*=2;
            }
            if(std::fabs(y[0]-y_old)<eps){
                return std::pair<double,bool>(y[0], true);
            }
            y_old=y[0];
            h=h/2;
        }
        return std::pair<double, bool>(y[0], false);
    }
}

int main(){
     try {
        
        auto ExampleFunction=[](double x){return std::sin(x)/x;};
        auto result1 = Limit(ExampleFunction, 0);
        std::cout << "Limit as x approaches 0 for sin(x)/x: ";
        if (result1.second) {
            std::cout << result1.first << " (converged)\n";
        } else {
            std::cout << "Did not converge within the specified tolerance.\n";
        }

        
        auto result2 = Limit([](double x) { return 1 / x; }, std::numeric_limits<double>::infinity());
        std::cout << "Limit as x approaches infinity for 1/x: ";
        if (result2.second) {
            std::cout << result2.first << " (converged)\n";
        } else {
            std::cout << "Did not converge within the specified tolerance.\n";
        }

        
        auto result3 = Limit([](double x) { return 1 / (x - 2); }, 2);
        std::cout << "Limit as x approaches 2 for 1/(x-2): ";
        if (result3.second) {
            std::cout << result3.first << " (converged)\n";
        } else {
            std::cout << "Did not converge within the specified tolerance.\n";
        }
        
        auto result4 = Limit([](double x) { return (std::exp(x) - 1) / x; }, 0.0015);
        std::cout << "Limit as x approaches 0 for (e^x - 1) / x: ";
        if (result4.second) {
            std::cout << result4.first << " (converged)\n";
        } else {
            std::cout << "Did not converge within the specified tolerance.\n";
        }

        
        auto result5 = Limit([](double x) { return std::log(x); }, 1, 1, 1e-8, 15);
        std::cout << "Limit as x approaches 1 for ln(x): ";
        if (result5.second) {
            std::cout << result5.first << " (converged)\n";
        } else {
            std::cout << "Did not converge within the specified tolerance.\n";
        }

        
        auto result6 = Limit([](double x) { return 1 / (x + 1); }, -1);
        std::cout << "Limit as x approaches -1 for 1 / (x + 1): ";
        if (result6.second) {
            std::cout << result6.first << " (converged)\n";
        } else {
            std::cout << "Did not converge within the specified tolerance.\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}