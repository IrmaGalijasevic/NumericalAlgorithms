#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include<algorithm>

using std::vector;
bool Fpequal(double x, double y){
    double eps=10*std::numeric_limits<double>::epsilon()*(std::fabs(x)+std::fabs(y));
    return std::fabs(x-y)<=eps;
}
class Vector {
private:
  vector<double> v;

public:
  explicit Vector(int n) {
    if (n <= 0)
      throw std::range_error("Bad dimension");
    v.resize(n);
  }
  Vector(std::initializer_list<double> l) {
    if (l.size() == 0)
      throw std::range_error("Bad dimension");
    for (auto x : l) {
      v.push_back(x);
    }
  }
  int NElems() const { return v.size(); }
  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }
  double &operator()(int i) {
    if (i < 1 || i > v.size())
      throw std::range_error("Invalid index");
    return v[i - 1];
  }
  double operator()(int i) const {
    if (i < 1 || i > v.size())
      throw std::range_error("Invalid index");
    return v[i - 1];
  }
  double Norm() const {
    double suma = 0;
    for (int i = 0; i < v.size(); i++) {
      suma += v[i] * v[i];
    }
    return sqrt(suma);
  }
  friend double VectorNorm(const Vector &v);
  double GetEpsilon() const {
    return 10 * std::numeric_limits<double>::epsilon() * Norm();
  }
  void Print(char separator = '\n', double eps = -1) const {
    if (eps < 0)
      eps = GetEpsilon();
    for (int i = 0; i < v.size(); i++) {
      if (std::fabs(v[i]) < eps)
        std::cout << 0;
      else
        std::cout << v[i];
      if (separator != '\n' && i != v.size() - 1)
        std::cout << separator;
      else if(i!=v.size()-1)
        std::cout << separator;
    }
  }
  friend void PrintVector(const Vector &v, char separator, double eps);
  friend Vector operator+(const Vector &v1, const Vector &v2);
  Vector &operator+=(const Vector &v) {
    if (v.NElems() != this->NElems())
      throw std::domain_error("Incompatible formats");
    for (int i = 0; i < this->NElems(); i++) {
      this->v[i] += v[i];
    }
    return *this;
  }
  friend Vector operator-(const Vector &v1, const Vector &v2);
  Vector &operator-=(const Vector &v) {
    if (v.NElems() != this->NElems())
      throw std::domain_error("Incompatible formats");
    for (int i = 0; i < this->NElems(); i++) {
      this->v[i] -= v[i];
    }
    return *this;
  }
  friend Vector operator*(double s, const Vector &v);
  friend Vector operator*(const Vector &v, double s);
  Vector &operator*=(double s) {
    for (auto &x : v) {
      x *= s;
    }
    return *this;
  }
  friend double operator*(const Vector &v1, const Vector &v2);
  friend Vector operator/(const Vector &v, double s);
  Vector &operator/=(double s) {
      if(s==0.0) throw std::domain_error("Division by zero");
    for (auto &x : v) {
      x /= s;
    }
    return *this;
  }
  void Chop(double eps=-1){
      if(eps<0) eps=GetEpsilon();
      for(auto &x:v){
          if(std::fabs(x)<eps && !Fpequal(x, eps)) x=0;
      }
  }
  bool EqualTo(const Vector &v, double eps=-1) const{
      if(eps<0) eps=GetEpsilon();
      if(this->NElems()!=v.NElems()) return false;
      for(int i=0; i<this->NElems(); i++){
          if(eps<=std::fabs(v[i]-v.v[i])) return false;
      }
      return true;
  }
};
double VectorNorm(const Vector &v) { return v.Norm(); }
void PrintVector(const Vector &v, char separator = '\n', double eps = -1) {
  v.Print(separator, eps);
}
Vector operator+(const Vector &v1, const Vector &v2) {
  if (v1.NElems() != v2.NElems())
    throw std::domain_error("Incompatible formats");
  Vector v3(v1.NElems());
  for (int i = 0; i < v3.NElems(); i++) {
    v3[i] = v1[i] + v2[i];
  }
  return v3;
}

Vector operator-(const Vector &v1, const Vector &v2) {
  if (v1.NElems() != v2.NElems())
    throw std::domain_error("Incompatible formats");
  Vector v3(v1.NElems());
  for (int i = 0; i < v3.NElems(); i++) {
    v3[i] = v1[i] - v2[i];
  }
  return v3;
}
Vector operator*(double s, const Vector &v) {
  Vector v2(v.NElems());
  for (int i = 0; i < v2.NElems(); i++) {
    v2[i] = v[i] * s;
  }
  return v2;
}
Vector operator*(const Vector &v, double s) { return s * v; }
double operator*(const Vector &v1, const Vector &v2) {
  if (v1.NElems() != v2.NElems())
    throw std::domain_error("Incompatible formats");
  double rez = 0.0;
  for (int i = 0; i < v1.NElems(); i++)
    rez += v1[i]*v2[i];
  return rez;
}
Vector operator/(const Vector &v, double s) {
if(s==0) throw std::domain_error("Division by zero");
  Vector v2(v.NElems());
  for (int i = 0; i < v2.NElems(); i++) {
    v2[i] = v[i] / s;
  }
  return v2;
}

class Matrix {
private:
  vector<vector<double>> mat;

public:
  Matrix(int m, int n) {
    if (m <= 0 || n <= 0)
      throw std::range_error("Bad dimension");
    mat = vector<vector<double>>(m, std::vector<double>(n, 0.0));
  }
  Matrix(const Vector &v) {
    mat = vector<vector<double>>(1, vector<double>(v.NElems(), 0.0));
    for (int i = 0; i < v.NElems(); i++) {
      mat[0][i] = v[i];
    }
  }
  Matrix(std::initializer_list<std::vector<double>> l) {
    if (l.begin()->size() == 0) throw std::range_error("Bad dimension");
    
    int l_size = l.begin()->size();
    for (auto vec : l) {
      if (l_size != vec.size() || vec.size() == 0)
        throw std::logic_error("Bad matrix");
    }
    for (auto vec : l) {
      mat.push_back(vec);
    }
  }
  int NRows() const { return mat.size(); }
  int NCols() const { return mat[0].size(); }
  double *operator[](int i) { return &mat[i][0]; }
  const double *operator[](int i) const { return &mat[i][0]; }
  double &operator()(int i, int j) {
    if (i < 1 || i > mat.size() || j < 1 || j > mat[0].size())
      throw std::range_error("Invalid index");
    return mat[i - 1][j - 1];
  }
  double operator()(int i, int j) const {
    if (i < 1 || i > mat.size() || j < 1 || j > mat[0].size())
      throw std::range_error("Invalid index");
    return mat[i - 1][j - 1];
  }
  double Norm() const {
    double sum = 0;
    for (auto red : mat) {
      for (auto x : red) {
        sum += x*x;
      }
    }
    return std::sqrt(sum);
  }
  friend double MatrixNorm(const Matrix &m);
  double GetEpsilon() const {
    return 10 * std::numeric_limits<double>::epsilon() * Norm();
  }
  void Print(int width = 10, double eps = -1) const {
    if (eps < 0) {
      eps = GetEpsilon();
    }
    int width_neg=width;
    for (int i = 0; i < mat.size(); i++) {
      for (int j = 0; j < mat[0].size(); j++) {
        if (std::fabs(mat[i][j]) < eps) {
          std::cout << std::setw(width) << 0.0;
        } else {
            if(mat[i][j]<0){
                width_neg=width+1;
            }
            else {
                width_neg=width;
            }
          std::cout << std::setw(width_neg) << mat[i][j];
        }
      }
      std::cout << std::endl;
    }
  }
  friend void PrintMatrix(const Matrix &m, int width, double eps);
  friend Matrix operator+(const Matrix &m1, const Matrix &m2);
  Matrix &operator+=(const Matrix &m) {
    if (this->NCols() != m.NCols() || this->NRows() != m.NRows())
      throw std::domain_error("Incompatible formats");

    for (int i = 0; i < this->NRows(); i++) {
      for (int j = 0; j < this->NCols(); j++) {
        mat[i][j] += m[i][j];
      }
    }
    return *this;
  }
  friend Matrix operator-(const Matrix &m1, const Matrix &m2);
  Matrix &operator-=(const Matrix &m) {
    if (this->NCols() != m.NCols() || this->NRows() != m.NRows())
      throw std::domain_error("Incompatible formats");

    for (int i = 0; i < this->NRows(); i++) {
      for (int j = 0; j < this->NCols(); j++) {
        mat[i][j] -= m[i][j];
      }
    }
    return *this;
  }
  friend Matrix operator*(double s, const Matrix &m);
  friend Matrix operator*(const Matrix &m, double s);
  Matrix &operator*=(double s) {
    for (int i = 0; i < this->NRows(); i++) {
      for (int j = 0; j < this->NCols(); j++) {
        mat[i][j] *= s;
      }
    }
    return *this;
  }
  friend Matrix operator*(const Matrix &m1, const Matrix &m2);
  Matrix &operator*=(const Matrix &m) {
    if (this->NCols() != m.NRows())
      throw std::domain_error("Incompatible formats");
    Matrix temp(this->NRows(), m.NCols());

    for (int i = 0; i < this->NRows(); i++) {
      for (int j = 0; j < m.NCols(); j++) {
        for (int k = 0; k < this->NCols(); k++) {
          temp[i][j] += mat[i][k] * m.mat[k][j];
        }
      }
    }
    mat=temp.mat;

    return *this;
  }
  friend Vector operator*(const Matrix &m, const Vector &v);
  friend Matrix Transpose(const Matrix &m);
  void Transpose() {
    if (this->NCols() == this->NRows()) {
      for (int i = 0; i < this->NRows(); i++) {
        for (int j = i + 1; j < this->NCols(); j++) {
          std::swap(mat[i][j], mat[j][i]);
        }
      }
    } else {
      Matrix temp(this->NCols(), this->NRows());
      for (int i = 0; i < this->NRows(); i++) {
        for (int j = 0; j < this->NCols(); j++) {
          temp.mat[j][i] = mat[i][j];
        }
      }
      std::swap(this->mat, temp.mat);
    }
  }
  void Chop(double eps=-1){
      if(eps<0) eps=GetEpsilon();
      for(auto &red:mat){
          for(auto &x:red){
              if(std::fabs(x)<eps && !Fpequal(x, eps)) x=0;
          }
      }
  }
  bool EqualTo(const Matrix& m, double eps=-1)const{
      if(NRows()!= m.NRows() || NCols()!=m.NCols()) return false;
      if(eps<0) eps=GetEpsilon();
      for(int i=0; i<NRows(); i++){
          for(int j=0; j<NCols(); j++){
              if(eps<=std::fabs(mat[i][j]-m[i][j])){ return false;}
          }
      }
      return true;
  }
  friend Matrix LeftDiv(Matrix m1, Matrix m2);  
  friend Vector LeftDiv(Matrix m, Vector v);    
  friend Matrix operator /(const Matrix &m, double s);  
  Matrix &operator /=(double s){
      if(s==0.0) throw std::domain_error("Division by zero");
      for(int i=0; i<NRows(); i++){
          for(int j=0; j<NCols(); j++){
              mat[i][j]/=s;
          }
      }
      return *this;
  }  
  friend Matrix operator /(Matrix m1, Matrix m2);   
  Matrix &operator /=(Matrix m){
      if(m.NRows()!=m.NCols()) throw std::domain_error("Divisor matrix is not square");
      if(NCols()!=m.NCols()) throw std::domain_error("Incompatible formats");
        double eps=m.GetEpsilon();
      for(int k=0; k<m.NCols(); k++){
          int p=k;
          for(int i=k+1; i<m.NRows(); i++){
              if(std::fabs(m[k][p])<std::fabs(m[k][i])) p=i;
          }
          if(std::fabs(m[k][p])<eps && !Fpequal(std::fabs(m[k][p]),eps)) throw std::domain_error("Divisor matrix is singular");
          if(p!=k){
              for(int i=0; i<m.NRows(); i++){
                  std::swap(m[i][k],m[i][p]);
              }
              
          }
          for(int i=k+1; i<m.NCols(); i++){
              double temp=m[k][i]/m[k][k];
              for(int j=k+1; j<m.NRows(); j++){
                  m[j][i]-=temp*m[j][k];
              }
              for(int j=0; j<NRows(); j++){
                  mat[j][i]-=temp*mat[j][k];
              }
          }
      }
      Matrix x(NRows(), m.NRows());
      for(int k=0; k<NRows(); k++){
          for(int i=m.NCols()-1;i>=0; i--){
              long double temp=mat[k][i];
              for(int j=i+1; j<m.NRows(); j++){
                  temp-=m[j][i]*x[k][j];
              }
              x[k][i]=temp/m[i][i];
          }
      }
      *this=x;
      return *this;
  }    
  double Det() const{
      if(NRows()!=NCols()) throw std::domain_error("Matrix is not square");
      double d=1;
      Matrix x(*this);
      double eps=GetEpsilon();
      for(int k=0; k<NRows(); k++){
          int p=k;
          for(int i=k+1; i<NRows(); i++){
              if(std::fabs(x.mat[p][k])<std::fabs(x.mat[i][k]) && !Fpequal(std::fabs(x.mat[i][k]), std::fabs(x.mat[p][k]))) p=i;
          }
          if(std::fabs(x.mat[p][k])<eps && !Fpequal(eps, std::fabs(x.mat[p][k]))) return 0;
          if(p!=k){
              std::swap(x.mat[k],x.mat[p]);
              d=-d;
          }
          d*=x.mat[k][k];
          for(int i=k+1; i<NRows(); i++){
              double temp=x.mat[i][k]/x.mat[k][k];
              for(int j=k+1; j<NRows(); j++){
                  x.mat[i][j]-=temp*x.mat[k][j];
              }
          }
      }
      return d;
  }   
  friend double Det(Matrix m);  
  void Invert(){
       if (NCols() != NRows()) throw std::domain_error("Matrix is not square");
        if (std::fabs(Det()) < std::numeric_limits<double>::epsilon()) throw std::domain_error("Matrix is singular");

    for (int k = 0; k < NRows(); k++) {
      double temp = mat[k][k];
      mat[k][k] = 1;
      for (int j = 0; j < NRows(); j++) {
        mat[k][j] /= temp;
      }
      for (int i = 0; i < NRows(); i++) {
        if (i != k) {
          temp = mat[i][k];
          mat[i][k] = 0;
          for (int j = 0; j < NRows(); j++) {
            mat[i][j]-=temp * mat[k][j];
          }
        }
      }
    }
  }    
  friend Matrix Inverse(Matrix m);  
  void ReduceToRREF(){
      int k=-1;
      int l=-1;
      double eps=GetEpsilon();
      std::vector<bool> w(NCols(),false);
     
      while(k<NRows() && l<NCols()){
          l++;
          k++;
          double v=0;
          int p=0;
          double temp=0;
          while(v<eps && !Fpequal(v, eps) && l<NCols()){
              p=k;
              for(int i=k; i<NRows(); i++){
                  if(v<std::fabs((*this)[i][l]) && !Fpequal(v, std::fabs((*this)[i][l]))){
                      v=std::fabs((*this)[i][l]);
                      p=i;
                  }
              }
              if(v<eps && !Fpequal(v, eps)){
                  l++;
              }
          }
          
          if(l<NCols()){
              w[l]=true;
              if(p!=k){
                  std::swap((*this).mat[p], (*this).mat[k]);
              }
              temp=(*this)[k][l];
              for(int j=l; j<NCols(); j++){
                 (*this)[k][j]/=temp;
              }
              for(int i=0; i<NRows(); i++){
                  if(i!=k){
                      temp=(*this)[i][l];
                      for(int j=l; j<NCols(); j++){
                          (*this)[i][j]-=temp*(*this)[k][j];
                      }
                  }
              }
          }
          
      }
  }  
  friend Matrix RREF(Matrix m); 
  int Rank() const{
      Matrix m(*this);
      int k=-1;
      int l=-1;
      double eps=GetEpsilon();
      std::vector<bool> w(NCols(),false);
     
      while(k<NRows() && l<NCols()){
          l++;
          k++;
          double v=0;
          int p=0;
          double temp=0;
          while(v<eps && !Fpequal(v, eps) && l<NCols()){
              p=k;
              for(int i=k; i<NRows(); i++){
                  if(v<std::fabs(m[i][l]) && !Fpequal(v, std::fabs(m[i][l]))){
                      v=std::fabs(m[i][l]);
                      p=i;
                  }
              }
              if(v<eps && !Fpequal(v, eps)){
                  l++;
              }
          }
          
          if(l<NCols()){
              w[l]=true;
              if(p!=k){
                  std::swap(m.mat[p], m.mat[k]);
              }
              temp=m[k][l];
              for(int j=l; j<NCols(); j++){
                 m[k][j]/=temp;
              }
              for(int i=0; i<NRows(); i++){
                  if(i!=k){
                      temp=m[i][l];
                      for(int j=l; j<NCols(); j++){
                          m[i][j]-=temp*m[k][j];
                      }
                  }
              }
          }
          
      }
      return k;
  } 
  friend int Rank(Matrix m);
};

double MatrixNorm(const Matrix &m) { return m.Norm(); }
void PrintMatrix(const Matrix &m, int width = 10, double eps = -1) {
  if (eps < 0) {
    eps = m.GetEpsilon();
  }

  for (int i = 0; i < m.NRows(); i++) {
    for (int j = 0; j < m.NCols(); j++) {
      if (std::fabs(m[i][j]) < eps) {
        std::cout << std::setw(width) << 0.0;
      } else {
        std::cout << std::setw(width) << m[i][j];
      }
    }
    std::cout << std::endl;
  }
}
Matrix operator+(const Matrix &m1, const Matrix &m2) {
  if (m1.NCols() != m2.NCols() || m1.NRows() != m2.NRows())
    throw std::domain_error("Incompatible formats");
  Matrix m3(m1.NRows(), m1.NCols());
  for (int i = 0; i < m3.NRows(); i++) {
    for (int j = 0; j < m3.NCols(); j++) {
      m3[i][j] = m1[i][j] + m2[i][j];
    }
  }
  return m3;
}
Matrix operator-(const Matrix &m1, const Matrix &m2) {
  if (m1.NCols() != m2.NCols() || m1.NRows() != m2.NRows())
    throw std::domain_error("Incompatible formats");
  Matrix m3(m1.NRows(), m1.NCols());
  for (int i = 0; i < m3.NRows(); i++) {
    for (int j = 0; j < m3.NCols(); j++) {
      m3[i][j] = m1[i][j] - m2[i][j];
    }
  }
  return m3;
}
Matrix operator*(double s, const Matrix &m) {
  Matrix m2(m.NRows(), m.NCols());
  for (int i = 0; i < m.NRows(); i++) {
    for (int j = 0; j < m.NCols(); j++) {
      m2[i][j] = s * m[i][j];
    }
  }
  return m2;
}
Matrix operator*(const Matrix &m, double s) { return s * m; }
Matrix operator*(const Matrix &m1, const Matrix &m2) {
  if (m1.NCols() != m2.NRows())
    throw std::domain_error("Incompatible formats");
  Matrix m3(m1.NRows(), m2.NCols());
  for (int i = 0; i < m1.NRows(); i++) {
    for (int j = 0; j < m2.NCols(); j++) {
      for (int k = 0; k < m1.NCols(); k++) {
        m3[i][j] += m1[i][k] * m2[k][j];
      }
    }
  }
  return m3;
}
Vector operator*(const Matrix &m, const Vector &v) {
  if (m.NCols() != v.NElems())
    throw std::domain_error("Incompatible formats");
  Vector v2(m.NRows());
  for (int i = 0; i < m.NRows(); i++) {
    for(int j=0; j<m.NCols(); j++){
        v2[i]+=m.mat[i][j]*v[j];
    }
  }
  return v2;
}
Matrix Transpose(const Matrix &m) {
  Matrix m2(m.NCols(), m.NRows());
  for (int i = 0; i < m.NRows(); i++) {
    for (int j = 0; j < m.NCols(); j++) {
      m2[j][i] = m[i][j];
    }
  }
  return m2;
}
Matrix LeftDiv(Matrix m1, Matrix m2){
    if(m1.NRows()!=m1.NCols()) throw std::domain_error("Divisor matrix is not square");
    if(m2.NRows()!=m1.NCols()) throw std::domain_error("Incompatible formats");
    auto eps=m1.GetEpsilon();
    for(int k=0; k<m1.NRows(); k++){
        int p=k;
        for(int i=k+1; i<m1.NCols(); i++ ){
            if(std::fabs(m1[i][k])>std::fabs(m1[p][k])) p=i;
        }
        if(std::fabs(m1[p][k])<eps) throw std::domain_error("Divisor matrix is singular");
        if(p!=k){
            std::swap(m1.mat[k],m1.mat[p]);
            std::swap(m2.mat[k], m2.mat[p]);
        }
        for(int i=k+1; i<m1.NCols(); i++){
            double temp=m1[i][k]/m1[k][k];
            for(int j=k+1; j<m1.NCols(); j++){
                m1[i][j]-=temp*m1[k][j];
            }
            for(int j=0; j<m2.NCols(); j++){
                m2[i][j]-=temp*m2[k][j];
            }
        }
    }
    //supstitucija unazad
    Matrix x(m1.NRows(),m2.NCols());
    for(int k=0; k<m2.NCols(); k++){
        for( int i=m1.NRows()-1; i>=0; i--){
            double temp=m2[i][k];
            for( int j=i+1; j<m1.NCols(); j++){
                temp-=m1[i][j]*x[j][k];
            }
            x[i][k]=temp/m1[i][i];
        }
    }
    return x;
}
Vector LeftDiv(Matrix m, Vector v){
    if(m.NRows()!=m.NCols()) throw std::domain_error("Divisor matrix is not square");
    if(v.NElems()!=m.NCols()) throw std::domain_error("Incompatible formats");
    double eps=m.GetEpsilon();
    for(int k=0; k<m.NRows(); k++){
        int p=k;
        for(int i=k+1; i<m.NRows(); i++){
            if(std::fabs(m[i][k])>std::fabs(m[p][k])) p=i;
        }
        if(std::fabs(m[p][k])<eps) throw std::domain_error("Divisor matrix is singular");
        
        if(p!=k){
            std::swap(m.mat[k],m.mat[p]);
            std::swap(v[k], v[p]);
        }
        for(int i=k+1; i<v.NElems(); i++){
            double temp=m[i][k]/m[k][k];
            for(int j=k+1; j<m.NCols(); j++){
                m[i][j]-=temp*m[k][j];
            }
            v[i]-=temp*v[k];
        }
    }

    Vector x(v.NElems());
    for(int i=v.NElems()-1; i>=0; i--){
        double temp=v[i];
        for(int j=i+1; j<v.NElems(); j++){
            temp-=m[i][j]*x[j];
        }
        x[i]=temp/m[i][i];
    }
    
    return x;
}
Matrix operator /(const Matrix &m, double s){
    if(s==0) throw std::domain_error("Division by zero");
    Matrix x(m.NRows(),m.NCols());
    
    for(int i=0; i<m.NRows(); i++){
        for(int j=0; j<m.NCols(); j++){
            x[i][j]=m[i][j]/s;
        }
    }
    return x;
}
Matrix operator /(Matrix m1, Matrix m2){
    m1/=m2;
    Matrix m3=m1;
    return m3;
}
double Det(Matrix m){
    return m.Det();
}
Matrix Inverse(Matrix m){
    m.Invert();
    return m;
}
Matrix RREF(Matrix m){
    m.ReduceToRREF();
    return m;
}
int Rank(Matrix m){
    return m.Rank();
}

class LUDecomposer{
    std::vector<int> w;
    Matrix LU;
    public:
    LUDecomposer(Matrix m):LU(m), w(m.NCols()){
        if(m.NRows()!=m.NCols()) throw std::domain_error("Matrix is not square");
        long double s, temp;
        double eps=m.GetEpsilon();
        int p;
        
        for(int j=0; j<LU.NRows(); j++){
            for(int i=0; i<=j; i++){
                s= LU[i][j];
                for(int k=0; k<i-1; k++){
                    s-=LU[i][k]*LU[k][j];
                }
                LU[i][j]=s;
            }
            p=j;
            for(int i=j+1; i<LU.NRows(); i++ ){
                s=LU[i][j];
                for(int k=0; k<=j-1; k++){
                    s-=LU[i][k]*LU[k][j];
                }
                LU[i][j]=s;
                if(std::fabs(s)>std::fabs(LU[p][j])){
                    p=i;
                }
            }
            if(std::fabs(LU[p][j])<LU.GetEpsilon()){
                throw std::domain_error("Matrix is singular");
            }
            if(p!=j){
                for(int x=0; x<LU.NCols(); x++){
                    std::swap(LU[p][x],LU[j][x]);
                }
                
            }
            w[j]=p;
            temp=LU[j][j];
            for(int i=j+1; i<LU.NRows(); i++){
                LU[i][j]/=temp;
            }
        }
        
    }
    void Solve(const Vector &b, Vector &x) const{
        if(LU.NCols()!=x.NElems() || LU.NRows()!=b.NElems()) throw std::domain_error("Incompatible formats");
        //supstitucija unaprijed
        std::vector<double>y(LU.NRows());
        double s;
        int p;
        x=b;
        for(int i=0; i<LU.NRows(); i++){
            p=w[i];
            s=x[p];
            x[p]=x[i];
            for(int j=0; j<i-1; j++){
                s-=LU[i][j]*y[j];
            }
            y[i]=s;
        }
        //supstitucija unazad
        for(int i=LU.NRows()-1; i>=0; i--){
            s=y[i];
            for(int j=i+1; j<LU.NRows(); j++){
                s-=LU[i][j]*x[j];
            }
            x[i]=s/LU[i][i];
        }
    }
    Vector Solve(Vector b) const{
        if(b.NElems()!=LU.NRows()) throw std::domain_error("Incompatible formats");
        //supstitucija unaprijed
        std::vector<double>y(LU.NRows());
        Vector x(LU.NCols());
        double s;
        int p;
        for(int i=0; i<LU.NRows(); i++){
            p=w[i];
            s=b[p];
            b[p]=b[i];
            for(int j=0; j<i-1; j++){
                s-=LU[i][j]*y[j];
            }
            y[i]=s;
        }
        //supstitucija unazad
        for(int i=LU.NRows()-1; i>=0; i--){
            s=y[i];
            for(int j=i+1; j<LU.NRows(); j++){
                s-=LU[i][j]*x[j];
            }
            x[i]=s/LU[i][i];
        }
        return x;
    }
    void Solve(const Matrix &b, Matrix &x) const{
        if(b.NRows()!=LU.NCols() || x.NRows()!=LU.NCols() || x.NCols()!=b.NCols()) throw std::domain_error("Incompatible formats");
        //supstitucija unaprijed
        Matrix y(b.NRows(), b.NCols());
        long double s;
        int p;
        y=b;
        
        for(int k=0; k<b.NCols(); k++){

        
        for(int i=0; i<b.NRows(); i++){
            p=w[i];
            s=y[p][k];
            y[p][k]=y[i][k];
            for(int j=0; j<i; j++){
                s-=LU[i][j]*y[j][k];
            }
            y[i][k]=s;
        }
        //supstitucija unazad
        for(int i=b.NRows()-1; i>=0; i--){
            s=y[i][k];
            for(int j=i+1; j<b.NRows(); j++){
                s-=LU[i][j]*x[j][k];
            }
            x[i][k]=s/LU[i][i];
        }
        }
    }
    Matrix Solve(Matrix b) const{
        if(b.NRows()!=LU.NCols()) throw std::domain_error("Incompatible formats");
        //supstitucija unaprijed
        Matrix y(b.NRows(),b.NCols());
        Matrix x(b.NRows(),b.NCols());
        long double s;
        int p;
        y=b;
        for(int k=0; k<b.NCols(); k++){
        
        for(int i=0; i<b.NRows(); i++){
            p=w[i];
            s=y[p][k];
            y[p][k]=y[i][k];
            for(int j=0; j<i; j++){
                s-=LU[i][j]*y[j][k];
            }
            y[i][k]=s;
        }
        //supstitucija unazad
        for(int i=b.NRows()-1; i>=0; i--){
            s=y[i][k];
            for(int j=i+1; j<b.NRows(); j++){
                s-=LU[i][j]*x[j][k];
            }
            x[i][k]=s/LU[i][i];
        }
        }
        return x;
    }
    Matrix GetCompactLU() const{
        return LU;
    }
    Matrix GetL() const{
        Matrix L(LU.NRows(), LU.NRows());
        for(int i=0; i<LU.NRows(); i++){
            for(int j=0; j<LU.NRows(); j++){
                if(i==j)L[i][j]=1;
                else if(i>j)L[i][j]=LU[i][j];
                else L[i][j]=0;
            }
        }
        return L;
    }
    Matrix GetU() const{
        Matrix U(LU.NRows(), LU.NRows());
        for(int i=0; i<LU.NRows(); i++){
            for(int j=0; j<LU.NRows(); j++){
                if(i<=j)U[i][j]=LU[i][j];
                else U[i][j]=0;
            }
        }
        return U;
    }
    Vector GetPermuation() const{
        Vector x(w.size());
        for(int i=0; i<w.size(); i++){
            x[i]=w[i]+1;
        }
        return x;
    }
};

class QRDecomposer{
    private:
    std::vector<double> d;
    Matrix QR;
    public:
    QRDecomposer(Matrix m):QR(m.NRows(),m.NCols()){
        if(m.NRows()<m.NCols()) throw std::domain_error("Invalid matrix format");
        double s;
        double temp;
        d.resize(m.NCols());
        for(int k=0; k<m.NCols(); k++){
            s=0;
            for(int i=k; i<m.NRows(); i++){
                s+=m[i][k]*m[i][k];
            }
            s=std::sqrt(s);
            temp=std::sqrt(s*(s+std::fabs(m[k][k])));
            if(m[k][k]<0){
                s=-s;
            }
            if(std::fabs(temp)<m.GetEpsilon()){
                throw std::domain_error("Matrix is singular");
            }
            QR[k][k]=(m[k][k]+s)/temp;
            for(int i=k+1; i<m.NRows(); i++){
                QR[i][k]=m[i][k]/temp;
            }
            d[k]=-s;
            for(int j=k+1; j<m.NCols(); j++){
                s=0;
                for(int i=k; i<m.NRows(); i++){
                    s+=QR[i][k]*m[i][j];
                }
                for(int i=k; i<m.NRows(); i++){
                    m[i][j]-=s*QR[i][k];
                }
            }
        }
        for(int i=0; i<m.NRows(); i++){
            for(int j=0; j<m.NCols(); j++){
                if(j>i){
                    QR[i][j]=m[i][j];
                }
            }
        }
    }
    void Solve(const Vector &b, Vector &x) const{
        if(QR.NCols()!=QR.NRows()) throw std::domain_error("Matrix is not square");
        if(x.NElems()!=QR.NCols() || b.NElems()!=QR.NRows()) throw std::domain_error("Incompatible formats");
        Vector y(b);
        double s;
        //Q'b
        for(int k=0; k<QR.NCols(); k++){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*y[i];
            }
            for(int i=k; i<QR.NRows();i++){
                y[i]-=s*QR[i][k];
            }
        }
        //supstitucija unazad
        for(int i=QR.NCols()-1; i>=0; i--){
            s=y[i];
            for(int j=i+1; j<QR.NCols(); j++){
                s-=QR[i][j]*x[j];
            }
            x[i]=s/d[i];
        }

    }
    Vector Solve(Vector b) const{
        if(QR.NCols()!=QR.NRows()) throw std::domain_error("Matrix is not square");
        if(b.NElems()!=QR.NRows()) throw std::domain_error("Incompatible formats");
        Vector y(b);
        Vector x(b.NElems());
        double s;
        //Q'b
        for(int k=0; k<QR.NCols(); k++){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[k][i]*y[i];
            }
            for(int i=k; i<QR.NRows();i++){
                y[i]-=s*QR[k][i];
            }
        }
        //supstitucija unazad
        for(int i=QR.NCols()-1; i>=0; i--){
            s=y[i];
            for(int j=i+1; j<QR.NCols(); j++){
                s-=QR[i][j]*x[j];
            }
            x[i]=s/d[i];
        }
        return x;
    }
    void Solve(Matrix &b, Matrix &x) const{
        if(QR.NCols()!=QR.NRows()) throw std::domain_error("Matrix is not square");
        if(b.NRows()!=QR.NRows() || b.NCols()!=x.NCols() || x.NRows()!=QR.NCols()) throw std::domain_error("Incompatible formats");
        Matrix y(b);
        double s;
        //Q'B
        for(int j=0; j<y.NCols(); j++){
            for(int k=0; k<QR.NCols(); k++){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*y[i][j];
            }
            for(int i=k; i<QR.NRows();i++){
                y[i][j]-=s*QR[i][k];
            }
        }
        }
        //supstitucija unazad
        for(int k=0; k<x.NCols(); k++){
            for(int i=QR.NCols()-1; i>=0; i--){
            s=y[i][k];
            for(int j=i+1; j<QR.NCols(); j++){
                s-=QR[i][j]*x[j][k];
            }
            x[i][k]=s/d[i];
        }
        }
    }
    Matrix Solve(Matrix b) const{
        if(QR.NCols()!=QR.NRows()) throw std::domain_error("Matrix is not square");
        if(b.NRows()!=QR.NRows()) throw std::domain_error("Incompatible formats");
        Matrix y(b);
        Matrix x(QR.NCols(), b.NCols());
        double s;
        //Q'B
        for(int j=0; j<y.NCols(); j++){
            for(int k=0; k<QR.NCols(); k++){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*y[i][j];
            }
            for(int i=k; i<QR.NRows();i++){
                y[i][j]-=s*QR[i][k];
            }
        }
        }
        //supstitucija unazad
        for(int k=0; k<x.NCols(); k++){
            for(int i=QR.NCols()-1; i>=0; i--){
            s=y[i][k];
            for(int j=i+1; j<QR.NCols(); j++){
                s-=QR[i][j]*x[j][k];
            }
            x[i][k]=s/d[i];
        }
        }
        return x;
    }
    Vector MulQWith(Vector v) const{
        if(QR.NRows()!=v.NElems()) throw std::domain_error("Incompatible formats");
        double s;
        for(int k=v.NElems()-1; k>=0; k--){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*v[i];
            }
            for(int i=k; i<QR.NRows(); i++){
                v[i]-=s*QR[i][k];
            }
        }
        return v;
    }
    Matrix MulQWith(Matrix m) const{
        if(QR.NRows()!=m.NRows()) throw std::domain_error("Incompatible formats");
        double s;
        for(int j=0; j<m.NCols(); j++){
            for(int k=m.NRows()-1; k>=0; k--){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*m[i][j];
            }
            for(int i=k; i<QR.NRows(); i++){
                m[i][j]-=s*QR[i][k];
            }
        }
        }
        return m;
    }
    Vector MulQTWith(Vector v) const{
        if(QR.NRows()!=v.NElems()) throw std::domain_error("Incompatible formats");
        //Q'b
        for(int k=0; k<QR.NCols(); k++){
            double s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*v[i];
            }
            for(int i=k; i<QR.NRows();i++){
                v[i]-=s*QR[i][k];
            }
        }
        return v;
    }
    Matrix MulQTWith(Matrix m) const{
        if(QR.NRows()!=m.NRows()) throw std::domain_error("Incompatible formats");
        //Q'm
        double s;
        for(int j=0; j<m.NCols(); j++){
            for(int k=0; k<m.NRows(); k++){
            s=0;
            for(int i=k; i<QR.NRows(); i++){
                s+=QR[i][k]*m[i][j];
            }
            for(int i=k; i<QR.NRows(); i++){
                m[i][j]-=s*QR[i][k];
            }
        }
        }
        return m;
    }
    Matrix GetQ() const{
        Matrix Q(QR.NRows(),QR.NRows());
        for(int j=0; j<Q.NRows(); j++){
            for(int i=0; i<Q.NCols(); i++){
                Q[i][j]=0;
            }
            Q[j][j]=1;
            for(int k=QR.NCols()-1; k>=0; k--){
                double s=0;
                for(int i=k; i<QR.NRows(); i++){
                    s+=QR[i][k]*Q[i][j];
                }
                for(int i=k; i<QR.NRows(); i++){
                    Q[i][j]-=s*QR[i][k];
                }
            }
        }
        return Q;
    }
    Matrix GetR() const{
        Matrix R(QR.NRows(), QR.NCols());
        for(int i=0; i<R.NRows(); i++){
            for(int j=0; j<R.NCols(); j++){
                if(i==j)R[i][j]=d[i];
                else if(j>i)R[i][j]=QR[i][j];
                else R[i][j]=0;
            }
        }
        return R;
    }
};
int main() {
    
   Matrix m1={{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
    Matrix m2={{8, -1, 4}, {6, 2, 2}, {7, -3, 6}};
    Vector v1={3, -2, -1};
    try{
        std::cout << "LeftDiv test:" << std::endl;
        Matrix rez1=LeftDiv(m1,m2);
        rez1.Print(5);
        std::cout <<std::endl;
        Matrix rez2=LeftDiv(m1,v1);
        rez2.Print(5);
        std::cout <<std::endl;
        Matrix rez3={{3,8,1},{-4,1,1}, {-4,1,1}};
        rez3=LeftDiv(rez3,m1);
        rez3.Print();
        std::cout <<std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "/= and / test" << std::endl;
        Matrix rez1=m1/m2;
        rez1.Print(5);
        std::cout <<std::endl;
        Matrix rez3={{3,8,1},{-4,1,1}, {-4,1,1}};
        rez3/=m1;
        rez3.Print(5);
        std::cout <<std::endl;
        rez3/=0.5;
        rez3.Print(5);
        std::cout <<std::endl;
        Matrix rez2=m1/v1;
        rez2.Print(5);
        std::cout <<std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "Det test:" << std::endl;
        std::cout << Det(m1) << std::endl;
        Matrix rez3={{3,8,1},{-4,1,1}, {-4,1,1}};
        rez3.Det();
        rez3.Print(5);
        std::cout <<std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "Invert test:" << std::endl;
        Matrix rez1=Inverse(m1);
        rez1.Print(5);
        std::cout <<std::endl;
        Matrix rez2=m2;
        rez2.Invert();
        rez2.Print(5);
        std::cout <<std::endl;
        Matrix rez3={{3,8,1},{-4,1,1}, {-4,1,1}};
        rez3.Invert();
        rez3.Print();
        std::cout <<std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "Rref test:" << std::endl;
        Matrix rez1=RREF(m1);
        rez1.Print(5);
        std::cout <<std::endl;
        Matrix rez2=m2;
        rez2.ReduceToRREF();
        rez2.Print(5);
        std::cout <<std::endl;
        Matrix rez3={{3,8,1},{-4,1,1}, {-4,1,1}};
        rez3.ReduceToRREF();
        rez3.Print();
        std::cout <<std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "Rank test:" << std::endl;
        int rez1=Rank(m1);
    
        std::cout << rez1<<std::endl;
        Matrix rez2=m2;
        std::cout << rez2.Rank() << std::endl;
        Matrix rez3={{3,8,1},{-4,1,1}, {-4,1,1}};
        
        std::cout << rez3.Rank() <<std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "EqualTo test:" << std::endl;
        Matrix rez1=Inverse(m1);
        rez1.Invert();
        
        std::cout << rez1.EqualTo(m1)<<std::endl;
        Matrix rez2=Transpose(m2);
        rez2.Transpose();
        
        std::cout << rez2.EqualTo(m2) << std::endl;
        Vector v2=v1/3;
        std::cout << v2.EqualTo(v1) << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LU solve vektor test:" << std::endl;
        LUDecomposer lu(m1);
        Vector rez1=lu.Solve(v1);
        rez1.Print();
        std::cout << std::endl;
        Vector pom1={2,4,6};
        lu.Solve(pom1,pom1);
        pom1.Print();
        std::cout << std::endl;
        Vector pom2={3,4,5,6,7};
        Vector x(3);
        lu.Solve(pom2,x);
        pom2.Print();
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LU Get______ test:" << std::endl;
        LUDecomposer lu(m1);
        Matrix rez1=lu.GetCompactLU();
        rez1.Print(3);
        std::cout << std::endl;
        Matrix rez2=lu.GetL();
        rez2.Print(3);
        std::cout << std::endl;
        Matrix rez3=lu.GetU();
        rez3.Print(3);
        std::cout << std::endl;
        Vector rez4=lu.GetPermuation();
        rez4.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LU konstruktor test:" << std::endl;
        LUDecomposer lu(m1);
        lu.GetL().Print(3);
        std::cout << std::endl;
        lu.GetU().Print();
        std::cout << std::endl;
        Matrix pom={{1,2,3},{2,4,6},{3,6,9}};
        LUDecomposer lu2(pom);
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR solve vektor test:" << std::endl;
        QRDecomposer qr(m1);
        Vector rez1=qr.Solve(v1);
        rez1.Print(3);
        std::cout << std::endl;
        Vector rez2={1,2,2};
        Vector rez3(rez2.NElems());
        qr.Solve(rez2,rez3);
        rez3.Print(3);
        std::cout << std::endl;
        rez2={5,7,9,1,4,2,5};
        qr.Solve(rez2,rez2);
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR solve matrica test:" << std::endl;
        QRDecomposer qr(m1);
        Matrix rez1=qr.Solve(m2);
        rez1.Print(3);
        std::cout << std::endl;
        Matrix rez2={{1,2,2},{3,4,9}, {12,5,17}};
        Matrix rez3(3,3);
        qr.Solve(rez2,rez3);
        rez3.Print(3);
        std::cout << std::endl;
        rez2={{5,2,6,8,4}, {5,7,8,8,5}, {5,7,3,7,8}, {8,9,1,3,4}};
        qr.Solve(rez2,rez2);
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR Get_____ test:" << std::endl;
        QRDecomposer qr(m1);
        qr.GetQ().Print(3);
        std::cout << std::endl;
        qr.GetR().Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR MulWith vektor test:" << std::endl;
        QRDecomposer qr(m1);
        Vector rez1=qr.MulQWith(v1);
        rez1.Print(3);
        std::cout << std::endl;
        Vector rez2=qr.MulQTWith(v1);
        rez2.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR MulWith matrica test:" << std::endl;
        QRDecomposer qr(m1);
        Matrix rez1=qr.MulQWith(m2);
        rez1.Print(3);
        std::cout << std::endl;
        Matrix rez2=qr.MulQTWith(m2);
        rez2.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LU pivotizacija test:" << std::endl;
        Matrix temp1={{2,1,1,0}, {4,3,3,1}, {8,7,9,5}, {6,7,9,8}};
        LUDecomposer lu(temp1);
        lu.GetL().Print(3);
        std::cout << std::endl;
        lu.GetU().Print(3);
        std::cout << std::endl;
        Matrix temp2={{0,1},{2,1}};
        LUDecomposer lu2(temp2);
        lu2.GetU().Print(3);
        std::cout << std::endl;
        lu2.GetU().Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR pivotizacija test:" << std::endl;
        Matrix temp1={{2,1,1,0}, {4,3,3,1}, {8,7,9,5}, {6,7,9,8}};
        QRDecomposer qr(temp1);
        qr.GetR().Print(3);
        std::cout << std::endl;
        qr.GetQ().Print(3);
        std::cout << std::endl;
        Matrix temp2={{0,1,2},{3,4,5}, {6,7,8}};
        QRDecomposer qr2(temp2);
        qr2.GetQ().Print(3);
        std::cout << std::endl;
        qr2.GetR().Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "Inverse pivotizacija test:" << std::endl;
        Matrix temp1={{2,3,5},{2,3,7},{4,1,8}};
        Matrix rez1=Inverse(temp1);
        rez1.Print(3);
        std::cout << std::endl;
        Matrix temp2={{0,3,2},{4,6,1},{3,1,7}};
        Matrix rez2=Inverse(temp2);
        rez2.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LeftDiv pivotizacija test:" << std::endl;
        Matrix temp1={{0,3,2},{4,6,1},{3,1,7}};
        Matrix temp2={{4,1,5},{1,2,1},{8,7,9}};
        Matrix rez1=temp1*temp2;
        Matrix temp3={{2,3,5},{2,3,7},{4,1,8}};
        Matrix rez2=LeftDiv(temp3, rez1);
        rez2.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LU singularna matrica test:" << std::endl;
        Matrix temp1={{2,2,2}, {4,4,4}, {6,6,6}};
        LUDecomposer lu(temp1);
        lu.GetL().Print(3);
        std::cout << std::endl;
        lu.GetU().Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "QR singularna matrica test:" << std::endl;
        Matrix temp1={{2,2,2}, {4,4,4}, {6,6,6}};
        QRDecomposer qr(temp1);
        qr.GetQ().Print(3);
        std::cout << std::endl;
        qr.GetR().Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "Inverse singularna matrica test:" << std::endl;
        Matrix temp1={{2,2,2}, {4,4,4}, {6,6,6}};
        Matrix rez=Inverse(temp1);
        rez.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
    try{
        std::cout << "LeftDiv singularna matrica test:" << std::endl;
        Matrix temp1={{1,2},{2,2},{3,3}};
        Matrix temp2= {{1,2},{1,2}};
        Matrix rez=temp1/temp2;
        rez.Print(3);
        std::cout << std::endl;
    }catch(std::domain_error e){
        std::cout << e.what() << std::endl;
    }
  return 0;
}