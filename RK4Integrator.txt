#include<iostream>
#include<stdexcept>
#include<cmath>
#include<vector>

template<typename FunTip>
double RK4Step(FunTip f, double x, double y, double h){
    double k1=f(x,y);
    double k2=f(x+h/2, y+h*k1/2);
    double k3=f(x+h/2, y+h*k2/2);
    double k4=f(x+h, y+h*k3);
    return y+h*(k1+2*k2+2*k3+k4)/6;
}


template <typename FunType>
std::vector<std::pair<double, double>> RK4Integrator(FunType f, double x0,
double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false){
    if (eps <= 0 || h == 0)throw std::domain_error("Invalid parameters");
    if((h>0 && xmax<x0) || (h<0 && xmax>x0)) return {{x0,y0}};
    std::vector<std::pair<double, double>> c;
    if(adaptive==false){
        double x=x0;
        double y=y0;
        if(h>0){
            while(x<=xmax+eps){
                c.push_back({x,y});
                double k1=f(x,y);
                double k2=f(x+h/2, y+h*k1/2);
                double k3=f(x+h/2, y+h*k2/2);
                double k4=f(x+h, y+h*k3);
                y=y+h*(k1+2*k2+2*k3+k4)/6;
                x=x+h;
            }
        }
        else{
            while(x>=xmax-eps){
                c.push_back({x,y});
                double k1=f(x,y);
                double k2=f(x+h/2, y+h*k1/2);
                double k3=f(x+h/2, y+h*k2/2);
                double k4=f(x+h, y+h*k3);
                y=y+h*(k1+2*k2+2*k3+k4)/6;
                x=x+h;
            }
        }
    }
    else{
        if(h>0){
            double x=x0;
            double y=y0;
            c.push_back({x,y});
            while(x<=xmax+eps){
                double u=RK4Step(f, x, y, h/2);
                double v=RK4Step(f, x+h/2, u, h/2);
                double w=RK4Step(f, x, y, h);
                double temp=std::fabs(w-v)/h;
                if(temp<=eps){
                    x=x+h;
                    y=v;
                    c.push_back({x,y});
                }
                h=h*std::min(5.0, 0.9*std::pow(eps/temp, 1/4.));
            }
            if(std::fabs(c[c.size()-1].first-xmax)>eps){
                c.erase(c.begin());
                if(!c.empty()){
                 h=xmax-c[c.size()-1].first;
                double u=RK4Step(f, x, y, h/2);
                double v=RK4Step(f, x+h/2, u, h/2);
                double w=RK4Step(f, x, y, h);
                c[c.size()-1]={xmax,v};   
                }
            }
        }
        else{
            double x=x0;
            double y=y0;
            c.push_back({x,y});
            while(x>=xmax-eps){
                double u=RK4Step(f, x, y, h/2);
                double v=RK4Step(f, x+h/2, u, h/2);
                double w=RK4Step(f, x, y, h);
                double temp=std::fabs(w-v)/((-1)*h);
                if(temp<=eps){
                    x=x+h;
                    y=v;
                    c.push_back({x,y});
                }
                h=h*std::min(5.0, 0.9*std::pow(eps/temp, 1/4.));
            }
            if(std::fabs(c[c.size()-1].first-xmax)>eps){
                c.erase(c.begin());
                if(!c.empty()){
                 h=xmax-c[c.size()-1].first;
                double u=RK4Step(f, x, y, h/2);
                double v=RK4Step(f, x+h/2, u, h/2);
                double w=RK4Step(f, x, y, h);
                c[c.size()-1]={xmax,v};   
                }
                
            }
        }
    }
    return c;
}


template <typename FunType>
std::vector<std::pair<double, std::vector<double>>> RK4SystemIntegrator(FunType f, double x0, const std::vector<double>& y0, double xmax, double h) {
    if (y0.size() != f(x0, y0).size()) throw std::range_error("Incompatible formats");
    if ((h > 0 && xmax < x0) || (h < 0 && xmax > x0)) return {{x0, y0}};

    int n = y0.size();
    std::vector<double> y(y0);
    std::vector<std::pair<double, std::vector<double>>> c;

    std::vector<double> K1(n), K2(n), K3(n), K4(n), u(n), v(n);

    while ((h > 0 && x0 <= xmax) || (h < 0 && x0 >= xmax)) {
        c.push_back({x0, y});

        for (int k = 0; k < n; k++) {
            K1[k] = h * f(x0, y)[k];
            u[k] = y[k] + K1[k] / 2;
        }

        for (int k = 0; k < n; k++) {
            K2[k] = h * f(x0 + h / 2, u)[k];
            v[k] = y[k] + K2[k] / 2;
        }

        for (int k = 0; k < n; k++) {
            K3[k] = h * f(x0 + h / 2, v)[k];
            u[k] = y[k] + K3[k];
        }

        for (int k = 0; k < n; k++) {
            K4[k] = h * f(x0 + h, u)[k];
        }

        for (int k = 0; k < n; k++) {
            y[k] = y[k] + (K1[k] + 2 * K2[k] + 2 * K3[k] + K4[k]) / 6;
        }

        x0 = x0 + h;
    }

    return c;
}




double f(double x, double y) {
    return 3 * x + 2;
}


std::vector<double> f_system(double x, const std::vector<double>& y) {
    std::vector<double> res(2);
    res[0] = y[0] + y[1];
    res[1] = -2 * y[0] + y[1];
    return res;
}
int main() {
    
    double x0 = 0;
    double y0 = 1;
    double xmax = 1;
    double h = 0.1;
    std::vector<std::pair<double, double>> c = RK4Integrator(f, x0, y0, xmax, h);
    std::cout << "RK4Integrator (y' = 3x + 2):" << std::endl;
    std::cout << "x\ty" << std::endl;
    for (const auto &p : c) {
        std::cout << p.first << "\t" << p.second << std::endl;
    }

    
    x0 = 0;
    std::vector<double> y0_system = {1, 0};
    xmax = 1;
    h = 0.1;
    std::vector<std::pair<double, std::vector<double>>> c_system = RK4SystemIntegrator(f_system, x0, y0_system, xmax, h);
    std::cout << std::endl;
    std::cout << "RK4SystemIntegrator (y1' = y1 + y2, y2' = -2*y1 + y2):" << std::endl;
    std::cout << "x\ty1\ty2" << std::endl;
    for (const auto &p : c_system) {
        std::cout << p.first << "\t" << p.second[0] << "\t" << p.second[1] << std::endl;
    }

    return 0;
}