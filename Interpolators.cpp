#include<iostream>
#include <limits>
#include<stdexcept>
#include<cmath>
#include <utility>
#include<vector>
#include<algorithm>


bool Fpequal(double x, double y){
    double eps=10*std::numeric_limits<double>::epsilon()*(std::fabs(x)+std::fabs(y));
    return std::fabs(x-y)<=eps;
}
class AbstractInterpolator{
    public:
    AbstractInterpolator(const std::vector<std::pair<double, double>> &data):dataset(data){
        std::sort(dataset.begin(), dataset.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b){
            return a.first<b.first; 
        });
        for(int i=1; i<dataset.size(); i++){
            if(Fpequal(dataset[i-1].first, dataset[i].first)) throw std::domain_error("Invalid data set");
        }
    }
    virtual double operator()(double x) const = 0;
    protected:
    int Locate(double x) const{
        if(x<=dataset[0].first){
            index=0;
          return index;  
        } 
        if(x>dataset[dataset.size()-1].first){
            index=dataset.size()-1;
          return index+1;  
        } 
        //prethodni interval
        if(index>0 && x>dataset[index-1].first && x<=dataset[index].first){
            index--;
            return index+1;
        } 
        //zapamceni interval
        if(index<dataset.size()-1 && x>dataset[index].first && x<=dataset[index+1].first) return index+1;
        //sljedeci interval
        if(index<dataset.size()-2 && x>dataset[index+1].first && x<=dataset[index+2].first){
            index++;
            return index+1;
        }
        auto it=std::lower_bound(dataset.begin(), dataset.end(), x, [](const std::pair<double, double>& a, double b){
            return a.first<b;
        });
        index=std::distance(dataset.begin(), it)-1;
        return index+1;
    }
    private:
    mutable int index=0;
    protected:
    std::vector<std::pair<double, double>> dataset;
};




class LinearInterpolator: public AbstractInterpolator{
    public:
    LinearInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){}
    double operator()(double x) const override{
        int i=Locate(x);
        double y;
        if(i==dataset.size()){
            i-=2;
        }
        else if(i!=0){
            i--;
        }
        y=((dataset[i+1].first-x)/(dataset[i+1].first-dataset[i].first))*dataset[i].second
                +((x-dataset[i].first)/(dataset[i+1].first-dataset[i].first))*dataset[i+1].second;
        return y;
    }
};

class PolynomialInterpolator: public AbstractInterpolator{
    std::vector<double> q;
    
    public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){
        q.resize(dataset.size());
        for(int i=0; i<dataset.size(); i++){
            q[i]=(dataset[i].second);
        }

        for(int j=1; j<dataset.size(); j++){
            for(int i=dataset.size()-1; i>=j; i--){
                q[i]=(q[i]-q[i-1])/(dataset[i].first-dataset[i-j].first);
                
            }
        }
        
    };
    double operator()(double x) const override{
        double f=q[q.size()-1];
        
        for(int i=q.size()-1; i>=1; i--){
            f=f*(x-dataset[i-1].first)+q[i-1];
        }
        return f;
    }
    void AddPoint(const std::pair<double, double> &p){
        for(int i=0; i<dataset.size(); i++){
            if(Fpequal(dataset[i].first,p.first)) throw std::domain_error("Invalid point");
        }
        dataset.push_back(p);
        std::vector<double>y_n(dataset.size());
        for(int i=0; i<dataset.size(); i++){
            y_n[i]=dataset[i].second;
        }
        for(int j=0; j<dataset.size()-1; j++){
            for(int i=dataset.size()-1; i>=j+1; i--){
                y_n[i]=(y_n[i]-y_n[i-1])/(dataset[i].first-dataset[i-j-1].first);
            }
        }
        q.push_back(y_n[dataset.size()-1]);
    }
    std::vector<double> GetCoefficients() const{
        std::vector<double> p(dataset.size(), 0);
        std::vector<double> w(dataset.size()+1, 0);
        std::vector<double> v(dataset.size()+1);
        w[0]=1;
        for(int i=1; i<w.size();i++){
            w[i]=w[i-1];
            for(int j=i-1; j>=1; j--){
                w[j]=w[j-1]-dataset[i-1].first*w[j];
            }
            w[0]=-dataset[i-1].first*w[0];
        }
        for(int i=0; i<dataset.size(); i++){
            double alpha=1;
            for(int j=0; j<dataset.size(); j++){
                if(j!=i){
                    alpha=alpha*(dataset[i].first-dataset[j].first);
                }
            }
            alpha=dataset[i].second/alpha;
            for(int j=0; j<w.size(); j++){
                v[j]=w[j];
            }
            for(int j=w.size()-2; j>=0; j--){
                v[j]+=dataset[i].first*v[j+1];
                p[j]+=alpha*v[j+1];
            }
        }
        return p;
    }
};

class PiecewisePolynomialInterpolator:public AbstractInterpolator{
    int k;
    public:
    PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data,int order):AbstractInterpolator(data){
        if(order<1 || order>=data.size()) throw std::domain_error("Invalid order");
        k=order;
    }
    double operator()(double x) const override{
        int i=Locate(x);
        
        //uslov za tacke na pocetku niza tacaka
        if((k%2==1 && i-(k-1)/2<1) || (k%2==0 && i-k/2<1)){
            double s=0;
            for(int m=0; m<k+1; m++){
                double p=dataset[m].second;
                for(int j=0; j<k+1; j++){
                    if(j!=m){
                        p=p*(x-dataset[j].first)/(dataset[m].first-dataset[j].first);
                    }
                }
                s=s+p;
            }
            return s;
        }
        //uslov za tacke pri kraju niza tacaka
        else if((k%2==1 && i+(k+1)/2>dataset.size()) || (k%2==0 && i+k/2>dataset.size())){
            double s=0;
            for(int m=dataset.size()-k-1; m<dataset.size(); m++){
                double p=dataset[m].second;
                for(int j=dataset.size()-k-1; j<dataset.size(); j++){
                    if(j!=m){
                        p=p*(x-dataset[j].first)/(dataset[m].first-dataset[j].first);
                    }
                }
                s=s+p;
            }
            return s;
        }
        //ostali
        else{
            
            if(k%2==1){
                    double s=0;
                for(int m=i-(k-1)/2-1; m<i+(k+1)/2; m++){
                    double p=dataset[m].second;
                    for(int j=i-(k-1)/2-1; j<i+(k+1)/2; j++){
                        if(j!=m){
                            p=p*(x-dataset[j].first)/(dataset[m].first-dataset[j].first);
                        }
                    }
                    s=s+p;
                }
                return s;
            }
            else{
                double s=0;
                for(int m=i-k/2-1; m<i+k/2; m++){
                    double p=dataset[m].second;
                    for(int j=i-k/2-1; j<i+k/2; j++){
                        if(j!=m){
                            p=p*(x-dataset[j].first)/(dataset[m].first-dataset[j].first);
                        }
                    }
                    s=s+p;
                }
                return s;
            }
        }
    }
};

class SplineInterpolator: public AbstractInterpolator{
    std::vector<double> r;
    public:
    SplineInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){
        r.resize(dataset.size());
        std::vector<double> alpha(dataset.size()-2);
        r[0]=0;
        r[r.size()-1]=0;
        for(int i=0; i<alpha.size(); i++){
            alpha[i]=2*(dataset[i+2].first-dataset[i].first);
            r[i+1]=3*((dataset[i+2].second-dataset[i+1].second)/(dataset[i+2].first-dataset[i+1].first)
                        -(dataset[i+1].second-dataset[i].second)/(dataset[i+1].first-dataset[i].first));
        }
        for(int i=0; i<alpha.size()-1; i++){
            double temp=(dataset[i+2].first-dataset[i+1].first)/alpha[i];
            alpha[i+1]-=temp*(dataset[i+2].first-dataset[i+1].first);
            r[i+2]-=temp*r[i+1];
        }
        r[r.size()-2]/=alpha[alpha.size()-1];
        
        for(int i=alpha.size()-2; i>=0; i--){
            r[i+1]=(r[i+1]-(dataset[i+2].first-dataset[i+1].first)*r[i+2])/alpha[i];
        }
    }
    double operator()(double x) const override{
        int i=Locate(x);
        if(i==dataset.size()){
            i-=2;
        }
        else if(i!=0){
            i--;
        }
        
        double t=x-dataset[i].first;
        double delta_x=dataset[i+1].first-dataset[i].first;
        double s=(r[i+1]-r[i])/(3*delta_x);
        double q=(dataset[i+1].second-dataset[i].second)/delta_x-delta_x*(r[i+1]+2*r[i])/3;
        return dataset[i].second+t*(q+t*(r[i]+t*s));
    }
};

class BarycentricInterpolator: public AbstractInterpolator{
    std::vector<double> w;
    public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order):AbstractInterpolator(data){
        if(order<0 || order>dataset.size()) throw std::domain_error("Invalid order");
        w.resize(dataset.size(), 0);
        double p=0;
        for(int i=0; i<w.size(); i++){
            w[i]=0;
            for(int k=std::max(0,i-order); k<=std::min(i, (int)w.size()-order); k++){
                p=1.0;
                for(int j=k; j<=k+order; j++){
                    if(j!=i){
                        p=p/(dataset[i].first-dataset[j].first);
                    }
                }
                if(k%2==1){
                    p=-p;
                }
            }
            w[i]+=p;
        }
    }
    double operator()(double x) const override{
        double p=0;
        double q=0;
        for(int i=0; i<w.size(); i++){
            if(Fpequal(x,dataset[i].first)) return dataset[i].second;
            double u=w[i]/(x-dataset[i].first);
            p=p+u*dataset[i].second;
            q=q+u;
        }
        return p/q;
    }
    std::vector<double> GetWeights() const{
        return w;
    }
};

class TrigonometricInterpolator:public AbstractInterpolator{
    public:
    TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data):AbstractInterpolator(data){
        if(std::abs(dataset[0].second-dataset[dataset.size()-1].second)>1e-10)throw std::domain_error("Function is not periodic");
    }
    double operator()(double x) const override{
        double omega=2*M_PI/(dataset[dataset.size()-1].first-dataset[0].first);
        double y=0;
        if(dataset.size()%2==0){
            int m=(dataset.size()-2)/2;
            for(int k=0; k<2*m+1; k++){
                double temp=1;
                for(int j=0; j<2*m+1; j++){
                    if(j!=k){
                        if(Fpequal(std::sin(omega/2)*(dataset[k].first-dataset[j].first),0))throw std::domain_error("divided by zero");
                        temp*=(std::sin(omega/2)*(x-dataset[j].first))/(std::sin(omega/2)*(dataset[k].first-dataset[j].first));
                    }
                }
                y=y+dataset[k].second*temp;
            }
        }
        else{
            int m=(dataset.size()-1)/2;
            double phi=0;
            for(int k=0; k<2*m; k++){
                double temp=1;
                double alpha_k=0;
                for(int j=0; j<2*m; j++){
                    if(j!=k){
                        alpha_k=alpha_k+dataset[j].first;
                    }
                }
                alpha_k=(2*phi/omega)-alpha_k;
                for(int j=0; j<2*m; j++){
                    if(j!=k)
                    temp*=(std::sin(omega/2)*(x-dataset[j].first))/(std::sin(omega/2)*(dataset[k].first-dataset[j].first));
                }
                y=y+dataset[k].second*((std::sin(omega/2)*(x-alpha_k))/(std::sin(omega/2)*(dataset[k].first-alpha_k)))*temp;
            }
        }
        return y;
    }
};

int main(){
   
   std::vector<std::pair<double, double>> data;
   

   //lin
   auto quadratic_function = [](double x) { return 7*x * x + 2 * x + 13; };
   for (double x = -5.0; x <= 5.0; x += 0.5) {
        double y = quadratic_function(x);
        data.emplace_back(x, y);
    }
    
    
    try {
        LinearInterpolator lin_interp(data);

        std::cout << "Testing linear interpolation:\n";
        for (double x = 0.25; x <= 4.25; x += 1.0) {
            double y_interp =lin_interp(x);
            double y_actual = quadratic_function(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }
       //konstrukor baca izuzeatk, isti x 
        data.push_back(std::make_pair(3.0, quadratic_function(3.0)));

    } catch (const std::domain_error& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
    }

    data.clear();
    //poly
    auto cubic_function = [](double x) { return x * x * x - 2 * x * x + 3 * x + 1; };
   for (double x = -5.0; x <= 5.0; x += 0.5) {
        double y = cubic_function(x);
        data.emplace_back(x, y);
    }
    
    
    try {
        PolynomialInterpolator poly_interp(data);

        std::cout << "\n\nTesting polynomial interpolation:\n";
        for (double x = 0.25; x <= 4.25; x += 1.0) {
            double y_interp =poly_interp(x);
            double y_actual = cubic_function(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }

        poly_interp.AddPoint({1.23, cubic_function(1.23)});
        poly_interp.AddPoint({-2.62, cubic_function(-2.62)});
        poly_interp.AddPoint({0.9, cubic_function(0.9)});

        std::cout << "\n\nTesting polynomial interpolation (addPoint):\n";
        for (double x = 0.25; x <= 4.25; x += 1.0) {
            double y_interp =poly_interp(x);
            double y_actual = cubic_function(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }

        std::cout << "\n\nTesting polynomial interpolation (getCoefficient):\n";
        auto v=poly_interp.GetCoefficients();
        for(auto x:v){
            std::cout << x << std::endl;
        }


        poly_interp.AddPoint({-4.0,cubic_function(-4.0)});

    } catch (const std::domain_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    data.clear();

    //piecewise
    auto func = [](double x) -> double {
        return 1.0 / ((x - 0.3) * (x - 0.3) + 0.01) + 1.0 / ((x - 0.9) * (x - 0.9) + 0.04);
    };
    data.push_back(std::make_pair(1.0, func(1.0)));
    data.push_back(std::make_pair(0.25, func(0.25)));
    data.push_back(std::make_pair(3, func(3)));
    data.push_back(std::make_pair(2.0, func(2.0)));
    data.push_back(std::make_pair(2.5, func(2.5)));
    data.push_back(std::make_pair(4.0, func(4.0)));
    data.push_back(std::make_pair(3.5, func(3.5)));
    
    
    try {
        PiecewisePolynomialInterpolator ppoly_interp(data, 4);

        std::cout << "\n\nTesting piecewise polynomial interpolation:\n";
        for (double x = 0.5; x <= 3.25; x += 0.5) {
            double y_interp =ppoly_interp(x);
            double y_actual = func(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }

        PiecewisePolynomialInterpolator ppoly_interp2(data, -7);

    } catch (const std::domain_error& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
    }
    data.clear();

    //spline

    for (double x = 0.0; x <= 2 * M_PI; x += 0.1) {
        double y = std::sin(x);
        data.emplace_back(x, y);
    }
    
    try {
        SplineInterpolator spl_interp(data);

        std::cout << "\n\nTesting spline interpolation:\n";
        for (double x = 0.25; x <= 4.25; x += 1.0) {
            double y_interp =spl_interp(x);
            double y_actual = std::sin(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }

    } catch (const std::domain_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    data.clear();

    //bary
    auto func2=[](double x){return (x*x-4)/(x+2);};
    data.push_back(std::make_pair(1.0, func2(1.0)));
    data.push_back(std::make_pair(0.25, func2(0.25)));
    data.push_back(std::make_pair(3.7, func2(3.7)));
    data.push_back(std::make_pair(2.2,func2(2.2)));
    data.push_back(std::make_pair(2.5, func2(2.5)));
    data.push_back(std::make_pair(4.0,func2(4.0)));
    data.push_back(std::make_pair(3.53, func2(3.53)));
    
    try {
        BarycentricInterpolator bar_interp(data,5);

        std::cout << "\n\nTesting barycentric interpolation:\n";
        for (double x = 0.25; x <= 4.25; x += 1.12) {
            double y_interp =bar_interp(x);
            double y_actual = func2(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }
         std::cout << "\n\nTesting barycentric interpolation(getWeights):\n";
        auto v2=bar_interp.GetWeights();
        for(auto x: v2){
            std::cout << x << std::endl;
        }

        BarycentricInterpolator bar_interp2(data,data.size()+5);

    } catch (const std::domain_error& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
    }
    data.clear();
   //trig

   auto func3=[](double x){ return 6*std::cos(x/2)+7*std::sin(x)-3;};
    for (double x = 0.0; x <4 * M_PI; x += 0.8) {
        
        data.push_back(std::make_pair(x,func3(x)));
        
    }
    data.push_back(std::make_pair(4*M_PI, func3(4*M_PI)));
    
    
    try {
        TrigonometricInterpolator trig_interp(data);

        std::cout << "\n\nTesting trigonometric interpolation:\n";
        for (double x = 0.25; x <= 4.25; x += 1.0) {
            double y_interp = trig_interp(x);
            double y_actual = func3(x);
            
            std::cout << "(" << x << ")\n";
            std::cout << "  Interpolated value: " << y_interp << "\n";
            std::cout << "  Actual value:  " << y_actual << "\n";
        }

        data.resize(data.size()-1);

        TrigonometricInterpolator t2(data);

    } catch (const std::domain_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}