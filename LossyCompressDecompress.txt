#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>

void FFT(std::vector<double>& x, std::vector<std::complex<double>>& y, int N, int s=0,int d=0, int t=1){
    if(N==1){
        y[d]=x[s];
    }
    else{
        FFT(x,y,N/2,s,d,2*t);
        FFT(x,y,N/2,s+t,d+N/2,2*t);
        std::complex<double> mi=1;
        std::complex<double> w=std::pow(std::complex<double>(std::cos(2*M_PI/N), std::sin(2*M_PI/N)), -1);
        for(int k=d; k<d+(N/2); k++){
            std::complex<double> u=y[k];
            std::complex<double> v=mi*y[k+(N/2)];
            y[k]=u+v;
            y[k+(N/2)]=u-v;
            mi=mi*w;
        }
    }
}

std::vector<double> LossyCompress(std::vector<double> data, int new_size){
    if (new_size <2 || new_size > data.size()) throw std::range_error("Bad new size");
    if (!((data.size() & (data.size() - 1)) == 0))
    throw std::range_error("Data size must be a power of two");
    std::vector<double> v1(data.size());
    std::vector<std::complex<double>> v2(data.size());
    int N=data.size();
    for(int i=0; i<=N/2-1; i++){
        v1[i]=data[2*i];
    }
    for(int i=N/2; i<=N-1; i++){
        v1[i]=data[2*(N-i)-1];
    }
    FFT(v1,v2,N);
    std::vector<double> y(new_size);
    for(int i=0; i<new_size-1; i++){
        std::complex<double> w =std::complex<double>(std::cos(M_PI / N), std::sin(M_PI / N));
        w = std::pow(w, (-1 * i) / 2.);
        y[i]=(w*v2[i]).real();
    }
    y[new_size-1]=N;
    return y;
}
/*
void InvFFT(std::vector<std::complex<double>>& y, std::vector<std::complex<double>>& x, int N, int s=0,int d=0, int t=1){
    if(N==1){
        x[d]=y[s];
    }
    else{
        InvFFT(y,x,N/2,s,d,2*t);
        InvFFT(y,x,N/2,s+t,d+N/2,2*t);
        std::complex<double> mi=1;
        std::complex<double> w=std::complex<double>(std::cos(2*M_PI/N), std::sin(2*M_PI/N));
        for(int k=d; k<d+(N/2); k++){
            std::complex<double> u=x[k];
            std::complex<double> v=mi*x[k+(N/2)];
            x[k]=(u+v)/2.0;
            x[k+(N/2)]=(u-v)/2.0;
            mi=mi*w;
        }
    }
}
std::vector<double> LossyDecompress(std::vector<double> compressed){
    std::vector<double> rez;
    return rez;
}*/



void invFFT(std::vector<std::complex<double>> &y,std::vector<std::complex<double>> &x, int n, int s = 0, int d = 0,int t = 1) {
  if (n==1)
    x[d]=y[s];
  else {
    invFFT(y,x,n/2,s,d,t*2);
    invFFT(y,x,n/2,s+t,d+n/2,t*2);
    std::complex<double> mi=1;
    std::complex<double> w=std::complex<double>(std::cos(2*M_PI/n), std::sin(2*M_PI/n));
    for (int k=d; k<d+n/2; k++) {
      std::complex<double> c=x[k];
      std::complex<double> c2=mi*x[k+n/2];
      x[k]=(c+c2)/2.0;
      x[k+n/2]=(c-c2)/2.0;
      mi*=w;
    }
  }
}

std::vector<double> LossyDecompress(std::vector<double> compressed) {
  int n=compressed[compressed.size() - 1];
  if (n < 1 || n < compressed.size()) throw std::logic_error("Bad compressed sequence");
  if (!((n & n - 1) == 0)) throw std::range_error("Data size must be a power of two");
  std::vector<std::complex<double>> y(n);
  std::vector<std::complex<double>> x(n);
  y[0] = compressed[0];

  for (int k=1; k<compressed.size()-1; k++) {
    std::complex<double> w =std::complex<double>(std::cos(M_PI/n), std::sin(M_PI/n));
    w=std::pow(w,k/2.);
    y[k]=2.0*w*compressed[k];
  }
  invFFT(y, x, n);
  std::vector<double> result(n);
  for (int i=0; i<=n-1; i++) {
    if (i%2==0)
      result[i]=x[i/2].real();
    else
      result[i]=x[n-(i+1)/2].real();
  }

  return result;
}





int main() {
    
    int sample_rate = 128;  
    int period = 2 * M_PI;  
    std::vector<double> original_data;


    for (int i = 0; i < sample_rate; ++i) {
        double t = static_cast<double>(i) / sample_rate;
        double value = sin(t * period);
        original_data.push_back(value);
    }


    std::cout << "Original Data:" << std::endl;
    for (double value : original_data) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    try {
        
        int new_size = 64;  
        std::vector<double> compressed_data = LossyCompress(original_data, new_size);

        std::cout << "\nCompressed Data:" << std::endl;
        for (double value : compressed_data) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

        std::vector<double> decompressed_data = LossyDecompress(compressed_data);

        std::cout << "\nDecompressed Data:" << std::endl;
        for (double value : decompressed_data) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}