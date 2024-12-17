#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include<algorithm>

using std::vector;

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
      if (std::abs(v[i]) < eps)
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
      if(s==0) throw std::domain_error("Division by zero");
    for (auto &x : v) {
      x /= s;
    }
    return *this;
  }
  void Chop(double eps=-1){
      if(eps==-1) eps=GetEpsilon();
      for(auto &x:v){
          if(x<eps) x=0;
      }
  }
  bool EqualTo(const Vector &v, double eps=-1) const{
      if(eps==-1) eps=GetEpsilon();
      if(this->NElems()!=v.NElems()) return false;
      for(int i=0; i<this->NElems(); i++){
          if(std::abs(v[i]-v.v[i])>=eps) return false;
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

    for (int i = 0; i < mat.size(); i++) {
      for (int j = 0; j < mat[0].size(); j++) {
        if (std::abs(mat[i][j]) < eps) {
          std::cout << std::setw(width) << 0.0;
        } else {
          std::cout << std::setw(width) << mat[i][j];
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
      if(eps==-1) eps=GetEpsilon();
      for(auto &red:mat){
          for(auto &x:red){
              if(x<eps) x=0;
          }
      }
  }
  bool EqualTo(const Matrix& m, double eps=-1)const{
      if(this->NRows()!= m.NRows() || this->NCols()!=m.NCols()) return false;
      if(eps==-1) eps=GetEpsilon();
      for(int i=0; i<this->NRows(); i++){
          for(int j=0; j<this->NCols(); j++){
              if(std::abs(mat[i][j]-m.mat[i][j])>=eps) return false;
          }
      }
      return true;
  }
  friend Matrix LeftDiv(Matrix m1, Matrix m2);  
  friend Vector LeftDiv(Matrix m, Vector v);    
  friend Matrix operator /(const Matrix &m, double s);  
  Matrix &operator /=(double s){
      if(s==0) throw std::domain_error("Division by zero");
      for(int i=0; i<NRows(); i++){
          for(int j=0; j<NCols(); j++){
              mat[i][j]/=s;
          }
      }
      return *this;
  }  
  friend Matrix operator /(Matrix m1, Matrix m2);   
  Matrix &operator /=(Matrix m){
      if(NRows()!=NCols()) throw std::domain_error("Divisor matris is not square");
      if(m.NCols()!=NCols()) throw std::domain_error("Incompatible formats");
    auto eps=GetEpsilon();
      for(int k=0; k<NRows(); k++){
          int p=k;
          for(int i=k+1; i<NCols(); i++){
              if(std::abs(mat[k][i])>std::abs(mat[k][p])) p=i;
          }
          if(std::abs(mat[k][p])<eps) throw std::domain_error("Divisor matrix is singular");
          if(p!=k){
              for(int i=0; i<NRows(); i++){
                  std::swap(mat[i][k],mat[i][p]);
              }
              for(int i=0; i<m.NRows(); i++){
                  std::swap(m.mat[i][k],m.mat[i][p]);
              }
          }
          for(int i=k+1; i<NCols(); i++){
              double temp=mat[k][i]/mat[k][k];
              for(int j=k+1; j<NCols(); j++){
                  mat[j][i]-=temp*mat[j][k];
              }
              for(int j=k+1; j<NCols(); j++){
                  m[j][i]-=temp*m[j][k];
              }
          }
      }
      for(int k=0; k<m.NCols(); k++){
          for(int i=m.NRows()-1;i>=0; i--){
              double temp=m[k][i];
              for(int j=i+1; j<m.NRows(); j++){
                  temp-=mat[j][i]*m[k][j];
              }
              m[k][i]=temp/mat[i][i];
          }
      }
      return *this;
  }    
  double Det() const{
      if(NRows()!=NCols()) throw std::domain_error("Matrix is not square");
      double d=1;
      Matrix x=*this;
      auto eps=GetEpsilon();
      for(int k=0; k<NRows(); k++){
          int p=k;
          for(int i=0; i<NRows(); i++){
              if(std::abs(x.mat[i][k])>std::abs(x.mat[p][k])) p=i;
          }
          if(x.mat[p][k]<eps) return 0;
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
      if(NCols()!=NRows()) throw std::domain_error("Matrix is not square");
      std::vector<int> w(NRows(),0);
      double eps=GetEpsilon();
      for(int k=0; k<NRows(); k++){
          int p=k;
          for(int i=k+1; i<NRows(); i++){
              if(std::abs(mat[i][k])>std::abs(mat[p][k])) p=i;
          }
          if(std::abs(mat[p][k])<eps) throw std::domain_error("Matrix is singular");
          if(p!=k){
              std::swap(mat[k], mat[p]);
          }
          w[k]=p;
          double temp=mat[k][k];
          mat[k][k]=1;
          for(int j=0; j<NRows(); j++){
              mat[k][j]/=temp;
          }
          for(int i=0; i<NRows(); i++){
              if(i!=k){
                  temp=mat[i][k];
                  mat[i][k]=0;
                  for(int j=0; j<NRows(); j++){
                      mat[i][j]-=temp*mat[k][j];
                  }
              }
          }
      }
      for(int j=NRows()-1; j>=0; j--){
          int p=w[j];
          if(p!=j){
              for(int i=0; i<NRows(); i++){
                  std::swap(mat[i][p], mat[i][j]);
              }
          }
      }
  }    
  friend Matrix Inverse(Matrix m);  
  void ReduceToRREF(){
      int k=0;
      int l=0;
      double eps=GetEpsilon();
      std::vector<bool> w(NCols());
      for(int j=0; j<NCols(); j++){
          w[j]=false;
      }
      while(k<=NRows() && l<=NCols()){
          l++;
          k++;
          double v=0;
          int p=0;
          double temp=0;
          while(v<eps && l<=NCols()){
              p=k;
              for(int i=k; i<NRows(); i++){
                  if(std::abs(mat[i][l])>v){
                      v=std::abs(mat[i][l]);
                      p=i;
                  }
              }
              if(v<eps){
                  l++;
              }
          }
          if(l<=NCols()){
              w[l]=true;
              if(p!=k){
                  std::swap(mat[k], mat[p]);
              }
              temp=mat[k][l];
              for(int j=l; j<NCols(); j++){
                  mat[k][j]/=temp;
              }
              for(int i=0; i<NRows(); i++){
                  if(i!=k){
                      temp=mat[i][l];
                      for(int j=l; j<NRows(); j++){
                          mat[i][j]-=temp*mat[k][j];
                      }
                  }
              }
          }
      }
  }  
  friend Matrix RREF(Matrix m); 
  int Rank() const{
      Matrix m(*this);
      int k=0;
      int l=0;
      double eps=GetEpsilon();
      std::vector<bool> w(NCols());
      for(int j=0; j<NCols(); j++){
          w[j]=false;
      }
      while(k<=NRows() && l<=NCols()){
          l++;
          k++;
          double v=0;
          int p=0;
          double temp=0;
          while(v<eps && l<=NCols()){
              p=k;
              for(int i=k; i<NRows(); i++){
                  if(std::abs(m[i][l])>v){
                      v=std::abs(m[i][l]);
                      p=i;
                  }
              }
              if(v<eps){
                  l++;
              }
          }
          if(l<=NCols()){
              w[l]=true;
              if(p!=k){
                  std::swap(m.mat[k], m.mat[p]);
              }
              temp=m[k][l];
              for(int j=l; j<NCols(); j++){
                  m[k][j]/=temp;
              }
              for(int i=0; i<NRows(); i++){
                  if(i!=k){
                      temp=mat[i][l];
                      for(int j=l; j<NRows(); j++){
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
      if (std::abs(m[i][j]) < eps) {
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
    if(m2.NCols()!=m1.NRows()) throw std::domain_error("Incompatible formats");
    auto eps=m1.GetEpsilon();
    for(int k=0; k<m1.NRows(); k++){
        int p=k;
        for(int i=k+1; i<m1.NCols(); i++ ){
            if(std::abs(m1[i][k])>std::abs(m1[p][k])) p=i;
        }
        if(std::abs(m1[p][k])<eps) throw std::domain_error("Divisor matrix is singular");
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
    auto eps=m.GetEpsilon();
    for(int k=0; k<m.NRows(); k++){
        int p=k;
        for(int i=k+1; i<m.NRows(); i++){
            if(std::abs(m[i][k])>std::abs(m[p][k])) p=i;
        }
        if(std::abs(m[p][k])<eps) throw std::domain_error("Divisor matrix is singular");
        if(p!=k){
            std::swap(m.mat[k],m.mat[p]);
            std::swap(v[k], v[p]);
        }
        for(int i=k+1; i<m.NCols(); i++){
            double temp=m[i][k]/m[k][k];
            for(int j=k+1; m.NCols(); j++){
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
    return m1/=m2;
}
double Det(Matrix m){
    if(m.NRows()!=m.NCols()) throw std::domain_error("Matrix is not square");
      double d=1;
      auto eps=m.GetEpsilon();
      for(int k=0; k<m.NRows(); k++){
          int p=k;
          for(int i=0; i<m.NRows(); i++){
              if(std::abs(m.mat[i][k])>std::abs(m.mat[p][k])) p=i;
          }
          if(m.mat[p][k]<eps) return 0;
          if(p!=k){
              std::swap(m.mat[k],m.mat[p]);
              d=-d;
          }
          d*=m.mat[k][k];
          for(int i=k+1; i<m.NRows(); i++){
              double temp=m.mat[i][k]/m.mat[k][k];
              for(int j=k+1; j<m.NRows(); j++){
                  m.mat[i][j]-=temp*m.mat[k][j];
              }
          }
      }
      return d;
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
    Matrix LU=Matrix(0,0);
    public:
    LUDecomposer(Matrix m){
        if(m.NRows()!=m.NCols()) throw std::domain_error("Matrix is not square");
        double s, temp;
        int p;
        w.resize(m.NRows());
        for(int i=0; m.NRows(); i++){
            w[i]=i;
        }
        for(int j=0; j<m.NRows(); j++){
            for(int i=0; i<j; i++){
                s= m[i][j];
                for(int k=0; k<i-1; k++){
                    s-=m[i][k]*m[k][j];
                }
                m[i][j]=s;
            }
            p=j;
            for(int i=j+1; i<m.NRows(); i++ ){
                s=m[i][j];
                for(int k=1; k<j-1; k++){
                    s-=m[i][k]*m[k][j];
                }
                m[i][j]=s;
                if(std::abs(s)>std::abs(m[p][j])){
                    p=i;
                }
            }
            if(std::abs(m[p][j])<m.GetEpsilon()){
                throw std::domain_error("Matrix is singular");
            }
            if(p!=j){
                for(int x=0; x<m.NRows(); x++){
                    std::swap(m[p][x],m[j][x]);
                }
            }
            w[j]=p;
            temp=m[j][j];
            for(int i=j+1; i<m.NRows(); i++){
                m[i][j]/=temp;
            }
        }
        LU=m;
    }
    void Solve(const Vector &b, Vector &x) const{
        //supstitucija unaprijed
        std::vector<double>y(LU.NRows());
        double s,temp;
        int p;
        x=b;
        for(int i=0; i<LU.NRows(); i++){
            p=w[i];
            s=x[p];
            x[p]=x[i];
            for(int j=0; j<i-1; j++){
                s-=LU[i][j];
            }
            y[i]=s;
        }
        //supstitucija unazad
        for(int i=LU.NRows(); i>=0; i--){
            s=y[i];
            for(int j=i+1; j<LU.NRows(); j++){
                s-=LU[i][j]*x[j];
            }
            x[i]=s/LU[i][i];
        }
    }
    Vector Solve(Vector b) const{
        //supstitucija unaprijed
        std::vector<double>y(LU.NRows());
        Vector x(LU.NRows());
        double s,temp;
        int p;
        for(int i=0; i<LU.NRows(); i++){
            p=w[i];
            s=b[p];
            b[p]=b[i];
            for(int j=0; j<i-1; j++){
                s-=LU[i][j];
            }
            y[i]=s;
        }
        //supstitucija unazad
        for(int i=LU.NRows(); i>=0; i--){
            s=y[i];
            for(int j=i+1; j<LU.NRows(); j++){
                s-=LU[i][j]*x[j];
            }
            x[i]=s/LU[i][i];
        }
        return x;
    }
    void Solve(const Matrix &b, Matrix &x) const{
        //supstitucija unaprijed
        Matrix y(LU.NRows(), LU.NRows());
        double s,temp;
        int p;
        x=b;
        for(int k=0; k<x.NRows(); k++){
            for(int i=0; i<LU.NRows(); i++){
                p=w[i];
                s=x[k][p];
                x[k][p]=x[k][i];
                for(int j=0; j<i; j++){
                    s-=LU[i][j]*x[k][j];
                }
                y[k][i]=s;
            }
        }
        //supstitucija unazad
        for(int k=0; k<x.NRows(); k++){
            for(int i=LU.NRows()-1; i>=0; i--){
                s=y[k][i];
                for(int j=i+1; j<LU.NRows(); j++){
                    s-=LU[i][j]*x[k][j];
                }
               x[k][i]=s/LU[i][i]; 
            }
        }
    }
    Matrix Solve(Matrix b) const{
        //supstitucija unaprijed
        Matrix y(LU.NRows(), LU.NRows());
        Matrix x(LU.NRows(), LU.NRows());
        double s,temp;
        int p;
        x=b;
        for(int k=0; k<x.NRows(); k++){
            for(int i=0; i<LU.NRows(); i++){
                p=w[i];
                s=x[k][p];
                x[k][p]=x[k][i];
                for(int j=0; j<i; j++){
                    s-=LU[i][j]*x[k][j];
                }
                y[k][i]=s;
            }
        }
        //supstitucija unazad
        for(int k=0; k<x.NRows(); k++){
            for(int i=LU.NRows()-1; i>=0; i--){
                s=y[k][i];
                for(int j=i+1; j<LU.NRows(); j++){
                    s-=LU[i][j]*x[k][j];
                }
               x[k][i]=s/LU[i][i]; 
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
    Matrix R=Matrix(0,0);
    public:
    QRDecomposer(Matrix m){
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
            temp=std::sqrt(s*(s+std::abs(m[k][k])));
            if(m[k][k]<0){
                s=-s;
            }
            if(std::abs(temp)<m.GetEpsilon()){
                throw std::domain_error("Matrix is singular");
            }
            R[k][k]=(m[k][k]+s)/temp;
            for(int i=k+1; i<m.NRows(); i++){
                R[i][k]=m[i][k]/temp;
            }
            d[k]=-s;
            for(int j=k+1; m.NCols(); j++){
                s=0;
                for(int i=k; i<m.NRows(); i++){
                    s+=R[i][k]*m[i][j];
                }
                for(int i=k; i<m.NRows(); i++){
                    m[i][j]-=s*R[i][k];
                }
            }
        }
        for(int i=0; i<m.NRows(); i++){
            for(int j=0; j<m.NCols(); j++){
                if(j>i){
                    R[i][j]=m[i][j];
                }
            }
        }
    }
    void Solve(const Vector &b, Vector &x) const{
        if(R.NCols()!=R.NRows()) throw std::domain_error("Matrix is not square");
        Vector y=b;
        double s, temp;
        //Q'b
        for(int k=0; k<R.NCols(); k++){
            s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[k][i]*y[i];
            }
            for(int i=k; i<R.NRows();i++){
                y[i]-=s*R[k][i];
            }
        }
        //supstitucija unazad
        for(int i=R.NCols()-1; i>=0; i--){
            s=y[i];
            for(int j=i+1; j<R.NCols(); j++){
                s-=R[i][j]*x[j];
            }
            x[i]=s/d[i];
        }

    }
    Vector Solve(Vector b) const{
        if(R.NCols()!=R.NRows()) throw std::domain_error("Matrix is not square");
        Vector y=b;
        Vector x=b;
        double s, temp;
        //Q'b
        for(int k=0; k<R.NCols(); k++){
            s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[k][i]*y[i];
            }
            for(int i=k; i<R.NRows();i++){
                y[i]-=s*R[k][i];
            }
        }
        //supstitucija unazad
        for(int i=R.NCols()-1; i>=0; i--){
            s=y[i];
            for(int j=i+1; j<R.NCols(); j++){
                s-=R[i][j]*x[j];
            }
            x[i]=s/d[i];
        }
        return x;
    }
    void Solve(Matrix &b, Matrix &x) const{
        if(R.NCols()!=R.NRows()) throw std::domain_error("Matrix is not square");
        Matrix y=b;
        double s, temp;
        //Q'B
        for(int j=0; j<y.NRows(); j++){
            for(int k=0; k<R.NCols(); k++){
            s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[k][i]*y[j][i];
            }
            for(int i=k; i<R.NRows();i++){
                y[j][i]-=s*R[k][i];
            }
        }
        }
        //supstitucija unazad
        for(int k=0; k<x.NRows(); k++){
            for(int i=R.NCols()-1; i>=0; i--){
            s=y[k][i];
            for(int j=i+1; j<R.NCols(); j++){
                s-=R[i][j]*x[k][j];
            }
            x[k][i]=s/d[i];
        }
        }
    }
    Matrix Solve(Matrix b) const{
        if(R.NCols()!=R.NRows()) throw std::domain_error("Matrix is not square");
        Matrix y=b;
        Matrix x=b;
        double s, temp;
        //Q'B
        for(int j=0; j<y.NRows(); j++){
            for(int k=0; k<R.NCols(); k++){
            s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[k][i]*y[j][i];
            }
            for(int i=k; i<R.NRows();i++){
                y[j][i]-=s*R[k][i];
            }
        }
        }
        //supstitucija unazad
        for(int k=0; k<x.NRows(); k++){
            for(int i=R.NCols()-1; i>=0; i--){
            s=y[k][i];
            for(int j=i+1; j<R.NCols(); j++){
                s-=R[i][j]*x[k][j];
            }
            x[k][i]=s/d[i];
        }
        }
        return x;
    }
    Vector MulQWith(Vector v) const{
        if(R.NCols()!=v.NElems()) throw std::domain_error("Incompatible formats");
        double s;
        for(int k=v.NElems()-1; k>=0; k--){
            s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[i][k]*v[i];
            }
            for(int i=k; i<R.NRows(); i++){
                v[i]-=s*R[i][k];
            }
        }
        return v;
    }
    Matrix MulQWith(Matrix m) const{
        if(R.NCols()!=m.NRows()) throw std::domain_error("Incompatible formats");
        for(int j=0; j<m.NRows(); j++){
            for(int k=R.NCols()-1; k>=0; k--){
                double s=0;
                for(int i=k; i<R.NRows(); i++){
                    s+=R[i][k]*m[j][i];
                }
                for(int i=k; i<R.NRows(); i++){
                    m[j][i]-=s*R[i][k];
                }
            }
        }
        return m;
    }
    Vector MulQTWith(Vector v) const{
        //Q'b
        for(int k=0; k<R.NCols(); k++){
            double s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[k][i]*v[i];
            }
            for(int i=k; i<R.NRows();i++){
                v[i]-=s*R[k][i];
            }
        }
        return v;
    }
    Matrix MulQTWith(Matrix m) const{
        //Q'B
        for(int j=0; j<m.NRows(); j++){
            for(int k=0; k<R.NCols(); k++){
            double s=0;
            for(int i=k; i<R.NRows(); i++){
                s+=R[k][i]*m[j][i];
            }
            for(int i=k; i<R.NRows();i++){
                m[j][i]-=s*R[k][i];
            }
        }
        
        }
        return m;
    }
    Matrix GetQ() const{
        Matrix Q(R.NRows(),R.NRows());
        for(int j=0; j<Q.NRows(); j++){
            for(int i=0; i<Q.NCols(); i++){
                Q[i][j]=0;
            }
            Q[j][j]=1;
            for(int k=R.NCols()-1; k>=0; k--){
                double s=0;
                for(int i=k; i<R.NRows(); i++){
                    s+=R[i][k]*Q[i][j];
                }
                for(int i=k; i<R.NRows(); i++){
                    Q[i][j]-=s*R[i][k];
                }
            }
        }
        return Q;
    }
    Matrix GetR() const{
        Matrix R1(R.NRows(), R.NCols());
        for(int i=0; i<R.NRows(); i++){
            for(int j=0; j<R.NCols(); j++){
                if(i==j)R1[i][j]=d[i];
                else if(j>i)R1[i][j]=R[i][j];
                else R1[i][j]=0;
            }
        }
        return R1;
    }
};
int main() {
  try {

    Vector v1(3);
    v1[0] = 1.0;
    v1[1] = 2.0;
    v1(3) = 3.0;
    Vector v2{2.1, 4.3, 8.7};
    Vector v3 = v1 + v2;
    Vector v4 = v1 - v2;
    Vector v5 = 6.2 * v1;
    Vector v6 = v1 * 7.3;
    Vector v7 = v2 / 2;
    double dot_product = v1 * v2;
    int v1_num_of_elements = v1.NElems();
    std::cout << "Vector v1: ";
    PrintVector(v1);
    std::cout << "Vector v1, number of elements: ";
    std::cout << v1_num_of_elements << std::endl;
    std::cout << "Vector v2: ";
    PrintVector(v2);
    std::cout << "Vector v3 (v1 + v2): ";
    PrintVector(v3);
    std::cout << "Vector v4 (v1 - v2): ";
    PrintVector(v4);
    std::cout << "Vector v5 (2.0 * v1): ";
    PrintVector(v5);
    std::cout << "Vector v6 (v1 * 2.0): ";
    PrintVector(v6);
    std::cout << "Dot product (v1 * v2): " << dot_product << std::endl;
    std::cout << "Vector v1 Norm: " << VectorNorm(v1) << std::endl;
    std::cout << "Epsilon for v1: " << v1.GetEpsilon() << std::endl;

    v3 += v1;
    v4 -= v1;
    v5 *= 19;
    v6 /= 0.53;

    std::cout << "Vector v3 (+=v1): ";
    v3.Print();
    std::cout << "Vector v4 (-=v1): ";
    v4.Print();
    std::cout << "Vector v5 (*=19): ";
    v5.Print('#');
    std::cout << "Vector v6 (/=0.53): ";
    v6.Print('_', 0.012);

    Matrix m1{{1.0, 2.0}, {3.0, 4.0}};
    Matrix m2(2, 2);
    m2[0][0] = 6.2;
    m2[0][1] = 5.3;
    m2(2, 1) = 9.1;
    m2(2, 2) = 3.6;
    Matrix m3 = m1 + m2;
    Matrix m4 = m1 - m2;
    Matrix m5 = 2.0 * m1;
    Matrix m6 = m1 * 2.0;
    Matrix m7 = m1 * m2;
    Vector v8{1.0, 2.0};
    Vector v9 = m1 * v8;
    Matrix m1Transposed = Transpose(m1);
    m1.Transpose();

    std::cout << "Matrix m1: " << std::endl;
    PrintMatrix(m1);
    std::cout << "Matrix m2: " << std::endl;
    PrintMatrix(m2);
    std::cout << "Matrix m3 (m1 + m2): " << std::endl;
    PrintMatrix(m3);
    std::cout << "Matrix m4 (m1 - m2): " << std::endl;
    PrintMatrix(m4);
    std::cout << "Matrix m5 (2.0 * m1): " << std::endl;
    PrintMatrix(m5);
    std::cout << "Matrix m6 (m1 * 2.0): " << std::endl;
    PrintMatrix(m6);
    std::cout << "Matrix m7 (m1 * m2): " << std::endl;
    PrintMatrix(m7);
    std::cout << "Matrix m1 Transposed: " << std::endl;
    PrintMatrix(m1Transposed);
    std::cout << "Vector v7: ";
    PrintVector(v8);
    std::cout << "Vector v8 (m1 * v7): ";
    PrintVector(v9);

  } catch (const std::exception &ex) {
    std::cout << "Exception: " << ex.what() << std::endl;
  }
  return 0;
}