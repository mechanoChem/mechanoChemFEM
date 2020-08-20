#ifndef EXTRAFUNCTIONS
#define EXTRAFUNCTIONS

#include <string>
#include <deal.II/base/table.h>
#include <Sacado.hpp>

/************************************
 * Convert comma-delimited string to vector of doubles
***********************************/
std::vector<double> str2vec(std::string str){

  while (str.find(",") != -1){
    str.replace(str.find(","),1," ");
  }
  std::istringstream iss(str);
  std::vector<double> vec;
  double num;
  while(iss >> num){
    vec.push_back(num);
  }

  return vec;
}

/************************************
 * Tensor operations on deal.II Tables
 ***********************************/
template<typename T>
void copy_to_q(dealii::Table<3,T> &A, const dealii::Table<2,T> &B, unsigned int q){
  for (unsigned int i=0; i<B.n_rows(); ++i){
    for (unsigned int j=0; j<B.n_cols(); ++j){
      A[q][i][j] = B[i][j];
    }
  }
}

template<typename T>
void copy_to_q(dealii::Table<2,T> &A, const dealii::Table<1,T> &B, unsigned int q){
  for (unsigned int i=0; i<B.size(0); ++i){
    A[q][i] = B[i];
  }
}

template<typename T> 
dealii::Table<2,T> copy_from_q(const dealii::Table<3,T> &A, unsigned int q){
  dealii::Table<2,T> B(A.size(1),A.size(2));
  for (unsigned int i=0; i<B.n_rows(); ++i){
    for (unsigned int j=0; j<B.n_cols(); ++j){
      B[i][j] = A[q][i][j];
    }
  }
  return B;
}

template<typename T>
dealii::Table<1,T> copy_from_q(const dealii::Table<2,T> &A, unsigned int q){
  dealii::Table<1,T> B(A.size(1));
  for (unsigned int i=0; i<B.size(0); ++i){
    B[i] = A[q][i];
  }
  return B;
}

template<typename T>
dealii::Table<2,T> operator*(const dealii::Table<2,T> &A, const dealii::Table<2,T> &B){
  unsigned int m1, n1, m2, n2;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  n2 = B.size(1);
  //if (n1 != m2){
  //  throw std::invalid_argument( "incompatible dimensions for multiplication" );
  //}
  //else{
    dealii::Table<2,T> C(m1,n2);
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n2; ++j){
	C[i][j] = 0;
	for (unsigned int k=0; k<m2; ++k){
	  C[i][j] += A[i][k]*B[k][j];
	}
      }
    }
    return C;
    //}
  
}

template<typename T>
dealii::Table<1,T> operator*(const dealii::Table<2,T> &A, const dealii::Table<1,T> &B){
  unsigned int m1, n1, m2;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  //if (n1 != m2){
  //  throw std::invalid_argument( "incompatible dimensions for multiplication" );
  //}
  //else{
    dealii::Table<1,T> C(m1);
    for (unsigned int i=0; i<m1; ++i){
      C[i] = 0;
      for (unsigned int k=0; k<m2; ++k){
	C[i] += A[i][k]*B[k];
      }
    }
    return C;
    //}
  
}

template<typename T>
dealii::Table<1,T> operator*(T c, const dealii::Table<1,T> &A){

  unsigned int m = A.size(0);
  dealii::Table<1,T> B(m);
  for (unsigned int i=0; i<m; ++i){
    B[i] = c*A[i];
  }
  return B;
  
}

template<typename T>
dealii::Table<1,T> operator*(double c, const dealii::Table<1,T> &A){

  unsigned int m = A.size(0);
  dealii::Table<1,T> B(m);
  for (unsigned int i=0; i<m; ++i){
    B[i] = c*A[i];
  }
  return B;
  
}

template<typename T>
dealii::Table<1,T> operator*(const dealii::Table<1,T> &A, double c){

  unsigned int m = A.size(0);
  dealii::Table<1,T> B(m);
  for (unsigned int i=0; i<m; ++i){
    B[i] = c*A[i];
  }
  return B;
  
}

template<typename T>
dealii::Table<1,T> operator+(const dealii::Table<1,T> &A, const dealii::Table<1,T> &B){

  unsigned int m = A.size(0);
  dealii::Table<1,T> C(m);
  for (unsigned int i=0; i<m; ++i){
    C[i] = A[i] + B[i];
  }
  return C;
  
}

template<typename T>
dealii::Table<2,T> operator*(double c, const dealii::Table<2,T> &A){

  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,T> B(m,n);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[i][j] = c*A[i][j];
    }
  }
  return B;
  
}

template<typename U>
dealii::Table<2,U> operator*(U c, const dealii::Table<2,double> &A){

  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,U> B(m,n);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[i][j] = c*A[i][j];
    }
  }
  return B;
  
}

template<typename T>
dealii::Table<2,T> operator*(T c,const dealii::Table<2,T> &A){

  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,T> B(m,n);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[i][j] = c*A[i][j];
    }
  }
  return B;
  
}

template<typename T>
dealii::Table<2,T> operator*(const dealii::Table<2,T> &A, T c){

  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,T> B(m,n);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[i][j] = c*A[i][j];
    }
  }
  return B;
  
}

template<typename T>
dealii::Table<2,T> operator*(const dealii::Table<2,T> &A, double c){

  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,T> B(m,n);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[i][j] = c*A[i][j];
    }
  }
  return B;
  
}


template<typename U>
dealii::Table<2,U> operator*(const dealii::Table<2,double> &A, U c){

  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,U> B(m,n);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[i][j] = c*A[i][j];
    }
  }
  return B;
  
}


template<typename T, typename U>
  dealii::Table<4,T> operator*(U c, const dealii::Table<4,T> &A){
  unsigned int I = A.size(0);
  unsigned int J = A.size(1);
  unsigned int K = A.size(2);
  unsigned int L = A.size(3);
  dealii::Table<4,T> B(I,J,K,L);
  for (unsigned int i=0; i<I; ++i){
    for (unsigned int j=0; j<J; ++j){
      for (unsigned int k=0; k<K; ++k){
	for (unsigned int l=0; l<L; ++l){
	  B[i][j][k][l] = c*A[i][j][k][l];
	}
      }
    }
  }
  return B;
  
}

template<typename T, typename U>
  dealii::Table<4,T> operator*(const dealii::Table<4,T> &A, U c){
  unsigned int I = A.size(0);
  unsigned int J = A.size(1);
  unsigned int K = A.size(2);
  unsigned int L = A.size(3);
  dealii::Table<4,T> B(I,J,K,L);
  for (unsigned int i=0; i<I; ++i){
    for (unsigned int j=0; j<J; ++j){
      for (unsigned int k=0; k<K; ++k){
	for (unsigned int l=0; l<L; ++l){
	  B[i][j][k][l] = c*A[i][j][k][l];
	}
      }
    }
  }
  return B;
  
}

template<typename T>
dealii::Table<2,T> operator+(const dealii::Table<2,T> &A, const dealii::Table<2,T> &B){
  unsigned int m1, n1, m2, n2;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  n2 = B.size(1);
  //if (n1 != n2 || m1 != m2){
  //  throw std::invalid_argument( "incompatible dimensions for addition" );
  //}
  //else{
    dealii::Table<2,T> C(m1,n2);
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n2; ++j){
	C[i][j] = A[i][j] + B[i][j];
      }
    }
    return C;
    //}
  
}

template<typename T>
dealii::Table<2,T> operator-(const dealii::Table<2,T> &A, const dealii::Table<2,T> &B){

  unsigned int m1, n1, m2, n2;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  n2 = B.size(1);
  //if (n1 != n2 || m1 != m2){
  //  throw std::invalid_argument( "incompatible dimensions for addition" );
  //}
  //else{
    dealii::Table<2,T> C(m1,n2);
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n2; ++j){
	C[i][j] = A[i][j] - B[i][j];
      }
    }
    return C;
    //}
  
}

template<typename T>
dealii::Table<4,T> operator+(const dealii::Table<4,T> &A, const dealii::Table<4,T> &B){
  unsigned int m1, n1, o1, p1, m2, n2, o2, p2;
  m1 = A.size(0);
  n1 = A.size(1);
  o1 = A.size(2);
  p1 = A.size(3);
  m2 = B.size(0);
  n2 = B.size(1);
  o2 = B.size(2);
  p2 = B.size(3);
  //if (n1 != n2 || m1 != m2 || o1 != o2 || p1 != p2){
  //  throw std::invalid_argument( "incompatible dimensions for addition" );
  //}
  //else{
    dealii::Table<4,T> C(m1,n1,o1,p1);
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n1; ++j){
	for (unsigned int k=0; k<o1; ++k){
	  for (unsigned int l=0; l<p1; ++l){
	    C[i][j][k][l] = A[i][j][k][l] + B[i][j][k][l];
	  }
	}
      }
    }
    return C;
    //}
  
}

template<typename T>
dealii::Table<4,T> operator-(const dealii::Table<4,T> &A, const dealii::Table<4,T> &B){
  unsigned int m1, n1, o1, p1, m2, n2, o2, p2;
  m1 = A.size(0);
  n1 = A.size(1);
  o1 = A.size(2);
  p1 = A.size(3);
  m2 = B.size(0);
  n2 = B.size(1);
  o2 = B.size(2);
  p2 = B.size(3);
  //if (n1 != n2 || m1 != m2 || o1 != o2 || p1 != p2){
  //  throw std::invalid_argument( "incompatible dimensions for subtraction" );
  //}
  //else{
    dealii::Table<4,T> C(m1,n1,o1,p1);
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n1; ++j){
	for (unsigned int k=0; k<o1; ++k){
	  for (unsigned int l=0; l<p1; ++l){
	    C[i][j][k][l] = A[i][j][k][l] - B[i][j][k][l];
	  }
	}
      }
    }
    return C;
    //}
  
}

template<typename T>
T trace(const dealii::Table<2,T> &A){
  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  T B = 0.;
  for (unsigned int i=0; i<m && i<n; ++i){
      B += A[i][i];
  }
  return B;
  
}

template<typename T>
dealii::Table<2,T> transpose(const dealii::Table<2,T> &A){
  unsigned int m = A.size(0);
  unsigned int n = A.size(1);
  dealii::Table<2,T> B(n,m);
  for (unsigned int i=0; i<m; ++i){
    for (unsigned int j=0; j<n; ++j){
      B[j][i] = A[i][j];
    }
  }
  return B;
  
}

template<typename T>
dealii::Table<2,T> double_contract(const dealii::Table<4,T> &A, const dealii::Table<2,T> &B){
  unsigned int o1, p1, m1, n1, m2, n2;
  m1 = A.size(0);
  n1 = A.size(1);
  o1 = A.size(2);
  p1 = A.size(3);
  m2 = B.size(0);
  n2 = B.size(1);
  //if (o1 != m2 || p1 != n2){
  //  throw std::invalid_argument( "incompatible dimensions for double contraction" );
  //}
  //else{
    dealii::Table<2,T> C(o1,p1);
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n1; ++j){
	C[i][j] = 0;;
	for (unsigned int k=0; k<o1; ++k){
	  for (unsigned int l=0; l<p1; ++l){
	    C[i][j] += A[i][j][k][l]*B[k][l];
	  }
	}
      }
    }
    return C;
    //}
}

template<typename T>
dealii::Table<2,T> double_contract(const dealii::Table<2,T> &A, const dealii::Table<4,T> &B){
  unsigned int m1, n1, m2, n2, o2, p2;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  n2 = B.size(1);
  o2 = B.size(2);
  p2 = B.size(3);
  //if (m1 != m2 || n1 != n2){
  //  throw std::invalid_argument( "incompatible dimensions for double contraction" );
  //}
  //else{
    dealii::Table<2,T> C(o2,p2);
    for (unsigned int i=0; i<o2; ++i){
      for (unsigned int j=0; j<p2; ++j){
	C[i][j] = 0;
	for (unsigned int k=0; k<m2; ++k){
	  for (unsigned int l=0; l<n2; ++l){
	    C[i][j] += A[k][l]*B[k][l][i][j];
	  }
	}
      }
    }
    return C;
    //}
}

template<typename T>
T double_contract(const dealii::Table<2,T> &A, const dealii::Table<2,T> &B){
  unsigned int m1, n1, m2, n2;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  n2 = B.size(1);
  //if (m1 != m2 || n1 != n2){
  //  throw std::invalid_argument( "incompatible dimensions for double contraction" );
  //}
  //else{
    T c = 0;
    for (unsigned int i=0; i<m1; ++i){
      for (unsigned int j=0; j<n1; ++j){
	c += A[i][j]*B[i][j];
      }
    }
    return c;
    //}
}

template<typename T>
T full_contract(const dealii::Table<2,T> &A, const dealii::Table<4,T> &B, const dealii::Table<2,T> &C){
  unsigned int m1, n1, m2, n2, o2, p2, o3, p3;
  m1 = A.size(0);
  n1 = A.size(1);
  m2 = B.size(0);
  n2 = B.size(1);
  o2 = B.size(2);
  p2 = B.size(3);
  o3 = C.size(0);
  p3 = C.size(1);
  //if (m1 != m2 || n1 != n2 || o2 != o3 || p2 != p3){
  //  throw std::invalid_argument( "incompatible dimensions for full contraction" );
  //}
  //else{
    T c=0;
    for (unsigned int i=0; i<m2; ++i){
      for (unsigned int j=0; j<n2; ++j){
	for (unsigned int k=0; k<o2; ++k){
	  for (unsigned int l=0; l<p2; ++l){
	    c += A[i][j]*B[i][j][k][l]*C[k][l];
	  }
	}
      }
    }
    return c;
    //}
}

#endif
