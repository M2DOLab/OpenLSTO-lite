template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::fill(const _Scalar &val) {
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           data[ii][jj] = val;
       }
   }
}

template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::setIdentity() {
   assert(_Rows == _Cols);
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           if (ii == jj) {
               data[ii][jj] = 1.0;
           }
           else {
               data[ii][jj] = 0.0;
           }
       }
   }
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar,_Rows,_Cols> Matrix<_Scalar,_Rows,_Cols>::Zero() {
   Matrix<_Scalar,_Rows,_Cols> outMat;
       for (uint ii = 0; ii < _Rows; ++ii) {
           for (uint jj = 0; jj < _Cols; ++jj) {
               outMat.data[ii][jj] = 1;
           }
       }
   return outMat;
}
template<typename _Scalar, int _Rows, int _Cols>
uint Matrix<_Scalar,_Rows,_Cols>::rows() const {
   return nRows;
}

template<typename _Scalar, int _Rows, int _Cols>
uint Matrix<_Scalar,_Rows,_Cols>::cols() const {
   return nCols;
}

template<typename _Scalar, int _Rows, int _Cols>
uint Matrix<_Scalar,_Rows,_Cols>::size() const {
   return nCols*nRows;
}

template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::print() const {
   std::cout << "\n [ ";
   for (uint jj = 0; jj < _Cols-1; ++jj) {
       std::cout << data[0][jj] << " , ";
   }
   std::cout << data[0][_Cols-1] << "\n";
   uint ii = 1;
   while (ii < _Rows-1) {
       std::cout << "   ";
       for (uint jj = 0; jj < _Cols-1; ++jj) {
           std::cout << data[ii][jj] << " , ";
       }
       std::cout << data[ii][_Cols-1] << "\n";
       ii++;
   }
   std::cout << "   ";
   for (uint jj = 0; jj < _Cols-1; ++jj) {
       std::cout << data[_Rows-1][jj] << " , ";
   }
   std::cout << data[_Rows-1][_Cols-1] << " ]\n";
}

// matrix operations
template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar,_Cols,_Rows> Matrix<_Scalar,_Rows,_Cols>::transpose() {
   Matrix<_Scalar,_Cols,_Rows> Matrix_T;
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           Matrix_T(jj,ii) = this->data[ii][jj];
       }
   }
   return Matrix_T;
}

template<typename _Scalar, int _Rows, int _Cols>
double Matrix<_Scalar,_Rows,_Cols>::determinant() {
   assert(this->nRows == this->nCols);
   double detVal; // determinant
   switch (this->nRows) {
   case 2:
       detVal = data[0][0]*data[1][1] - data[0][1]*data[1][0];
       return detVal;
       break;
   case 3:
       detVal = data[0][0]*(data[1][1]*data[2][2] - data[1][2]*data[2][1])
               -data[0][1]*(data[1][0]*data[2][2] - data[1][2]*data[2][0])
               +data[0][2]*(data[1][0]*data[2][1] - data[1][1]*data[2][0]);
       return detVal;
       break;
   default:
       std::cout << "Currently 2 or 3 dim Matrix are supported\n";
       return detVal;
       break;
   }
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar,_Rows,_Cols> Matrix<_Scalar,_Rows,_Cols>::inverse() {
   assert(this->nRows == this->nCols);
   double detVal = this->determinant(); // determinant
   Matrix<_Scalar,_Rows,_Cols> invMat;
   switch (this->nRows) {
       case 2:
           invMat.data[0][0] = +data[1][1]/detVal;
           invMat.data[0][1] = -data[0][1]/detVal;
           invMat.data[1][0] = -data[1][0]/detVal;
           invMat.data[1][1] = +data[0][0]/detVal;
           break;
       case 3:
           double a,b,c,d,e,f,g,h,i;
           a = data[0][0]; b = data[0][1]; c = data[0][2];
           d = data[1][0]; e = data[1][1]; f = data[1][2];
           g = data[2][0]; h = data[2][1]; i = data[2][2];

           invMat.data[0][0] = +(e*i - f*h)/detVal;
           invMat.data[0][1] = -(b*i - c*h)/detVal;
           invMat.data[0][2] = +(b*f - c*e)/detVal;
           invMat.data[1][0] = -(d*i - f*g)/detVal;
           invMat.data[1][1] = +(a*i - c*g)/detVal;
           invMat.data[1][2] = -(a*f - c*d)/detVal;
           invMat.data[2][0] = +(d*h - e*g)/detVal;
           invMat.data[2][1] = -(a*h - b*g)/detVal;
           invMat.data[2][2] = +(a*e - b*d)/detVal;
           break;
       default:
           std::cout << "Currently 2 or 3 dim Matrix are supported\n";
           break;
   }
   return invMat;
}

template<typename _Scalar, int _Rows, int _Cols>
template<typename _OtherScalar, int _OtherRows, int _OtherCols>
Matrix<_Scalar,_Rows,_OtherCols> Matrix<_Scalar,_Rows,_Cols>::dot(Matrix<_OtherScalar,_OtherRows,_OtherCols> X) {
   assert(_Cols == _OtherRows);
   Matrix<_Scalar,_Rows,_OtherCols> MatrixOut;
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _OtherCols; ++jj) {
           for (uint pp = 0; pp < _Cols; ++pp) {
               MatrixOut.data[ii][jj] += this->data[ii][pp]*X.data[pp][jj];
           }
       }
   }
   return MatrixOut;
}

template<typename _Scalar, int _Rows, int _Cols>
template<typename _OtherScalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,_Rows,_Cols>::dot(Matrix<_OtherScalar,-1,-1> X) {
   assert(_Cols == this->nRows);
   Matrix<_Scalar,-1,-1> MatrixOut;
   MatrixOut.resize(_Rows,this->nCols);
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           for (uint pp = 0; pp < _Cols; ++pp) {
               MatrixOut.data[ii][jj] += this->data[ii][pp]*X.data[pp][jj];
           }
       }
   }
   return MatrixOut;
}

template<typename _Scalar, int _Rows, int _Cols>
template<typename _OtherScalar, int _OtherRows>
Vector<_Scalar,_Rows> Matrix<_Scalar,_Rows,_Cols>::dot(Vector<_OtherScalar,_OtherRows> V) {
   Vector<_Scalar,_Rows> VectOut;
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint pp = 0; pp < _Cols; ++pp) {
           VectOut.data[ii] += this->data[ii][pp]*V.data[pp];
       }
   }
   return VectOut;

}

template<typename _Scalar, int _Rows, int _Cols>
template<typename _OtherScalar>
Vector<_Scalar,-1> Matrix<_Scalar,_Rows,_Cols>::dot(Vector<_OtherScalar,-1> V) {
   Vector<_Scalar,-1> VectOut(_Rows);
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint pp = 0; pp < _Cols; ++pp) {
           VectOut.data[ii] += this->data[ii][pp]*V.data[pp];
       }
   }
   return VectOut;

}

// Overloading operators
template<typename _Scalar, int _Rows, int _Cols>
_Scalar &Matrix<_Scalar,_Rows,_Cols>::operator () (int row, int col) {
   assert(col >= 0 && col < nCols);
   assert(row >= 0 && row < nRows);
   return data[row][col];
}

template<typename _Scalar, int _Rows, int _Cols>
const bool Matrix<_Scalar,_Rows,_Cols>::operator == (const Matrix<_Scalar,_Rows,_Cols> &rhs) {
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           if (this->data[ii][jj] != rhs(ii,jj)) {
               return false;
           }
       }
   }
   return true;
}

template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::operator = (const Matrix<_Scalar,_Rows,_Cols> &rhs) {
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           this->data[ii][jj] = rhs.data[ii][jj];
       }
   }
}

template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::operator += (const _Scalar &val) {
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           this->data[ii][jj] += val;
       }
   }
}

template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::operator *= (const _Scalar &val) {
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           this->data[ii][jj] *= val;
       }
   }
}

template<typename _Scalar, int _Rows, int _Cols>
void Matrix<_Scalar,_Rows,_Cols>::operator -= (const _Scalar &val) {
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           this->data[ii][jj] -= val;
       }
   }
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar,_Rows,_Cols> Matrix<_Scalar,_Rows,_Cols>::operator + (Matrix<_Scalar,_Rows,_Cols> rhs) {
   Matrix<_Scalar,_Rows,_Cols> MatrixOut;
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] + rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar, int _Rows, int _Cols>
Matrix<_Scalar,_Rows,_Cols> Matrix<_Scalar,_Rows,_Cols>::operator - (Matrix<_Scalar,_Rows,_Cols> rhs) {
   Matrix<_Scalar,_Rows,_Cols> MatrixOut;
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] - rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar, int _Rows, int _Cols>
template<typename _OtherScalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,_Rows,_Cols>::operator + (Matrix<_OtherScalar,-1,-1> rhs) {
   Matrix<_Scalar,-1,-1> MatrixOut;
   MatrixOut.resize(this->nRows, this->nCols);
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] + rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar, int _Rows, int _Cols>
template<typename _OtherScalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,_Rows,_Cols>::operator - (Matrix<_OtherScalar,-1,-1> rhs) {
   Matrix<_Scalar,-1,-1> MatrixOut;
   MatrixOut.resize(this->nRows,this->nCols);
   for (uint ii = 0; ii < _Rows; ++ii) {
       for (uint jj = 0; jj < _Cols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] - rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::resize(int _nRows, int _nCols) {
   data.clear();
   data.resize(_nRows,std::vector<_Scalar>(_nCols,0.0));
   nRows = _nRows; nCols = _nCols;
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::fill(_Scalar m) {
   for (uint ii = 0; ii < nRows; ++ii) {
       std::fill(data[ii].begin(),data[ii].end(),m);
   }
}

template<typename _Scalar>
uint Matrix<_Scalar,-1,-1>::rows() const {
   return this->nRows;
}

template<typename _Scalar>
uint Matrix<_Scalar,-1,-1>::cols() const {
   return this->nCols;
}

template<typename _Scalar>
uint Matrix<_Scalar,-1,-1>::size() const {
   return this->nRows*this->nCols;
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::print() const {
   std::cout << "\n [ ";
   for (uint jj = 0; jj < nCols-1; ++jj) {
       std::cout << data[0][jj] << " , ";
   }
   std::cout << data[0][nCols-1]<< "\n";
   uint ii = 1;
   while (ii < nRows-1) {
       std::cout << "   ";
       for (uint jj = 0; jj < nCols-1; ++jj) {
           std::cout << data[ii][jj]<< " , ";
       }
       std::cout << data[ii][nCols-1] << "\n";
       ii++;
   }
   std::cout <<"   ";
   for (uint jj = 0; jj < nCols-1; ++jj) {
       std::cout << data[nRows-1][jj] << " , ";
   }

   std::cout << data[nRows-1][nCols-1] << " ]\n";
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::setIdentity() {
   assert(rows() == cols());
   for (uint ii = 0; ii < rows(); ++ii) {
       for (uint jj = 0; jj < cols(); ++jj) {
           if (ii == jj) {
               data[ii][jj] = 1.0;
           }
           else {
               data[ii][jj] = 0.0;
           }
       }
   }
}

template<typename _Scalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::transpose() {
   Matrix<_Scalar,-1,-1 > outMatrix;
   outMatrix.resize(this->nCols,this->nRows);
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           outMatrix.data[jj][ii] = this->data[ii][jj];
       }
   }
   return outMatrix;
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::operator = (const Matrix<_Scalar,-1,-1> &rhs) {
   this->resize(rhs.rows(),rhs.cols());
   for (uint ii = 0; ii < rhs.rows(); ++ii) {
       for (uint jj = 0; jj < rhs.cols(); ++jj) {
           this->data[ii][jj] = rhs.data[ii][jj];
       }
   }
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::operator += (const _Scalar &val) {
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           this->data[ii][jj] += val;
       }
   }
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::operator -= (const _Scalar &val) {
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           this->data[ii][jj] -= val;
       }
   }
}

template<typename _Scalar>
void Matrix<_Scalar,-1,-1>::operator *= (const _Scalar &val) {
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           this->data[ii][jj] *= val;
       }
   }
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherRows, int _OtherCols>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::operator + (Matrix<_OtherScalar,_OtherRows,_OtherCols> rhs) {
   assert(this->nRows == rhs.nRows);
   assert(this->nCols == rhs.nCols);
   Matrix<_Scalar,-1,-1> MatrixOut;
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] + rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherRows, int _OtherCols>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::operator - (Matrix<_OtherScalar,_OtherRows,_OtherCols> rhs) {
   assert(this->nRows == rhs.nRows);
   assert(this->nCols == rhs.nCols);
   Matrix<_Scalar,-1,-1> MatrixOut;
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] - rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar>
template<typename _OtherScalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::operator + (Matrix<_OtherScalar,-1,-1> rhs) {
   assert(this->nRows == rhs.nRows);
   assert(this->nCols == rhs.nCols);
   Matrix<_Scalar,-1,-1> MatrixOut;
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] + rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar>
template<typename _OtherScalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::operator - (Matrix<_OtherScalar,-1,-1> rhs) {
   assert(this->nRows == rhs.nRows);
   assert(this->nCols == rhs.nCols);
   Matrix<_Scalar,-1,-1> MatrixOut;
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < this->nCols; ++jj) {
           MatrixOut.data[ii][jj] = this->data[ii][jj] - rhs.data[ii][jj];
       }
   }
   return MatrixOut;
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherRows, int _OtherCols >
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::dot(Matrix<_OtherScalar,_OtherRows,_OtherCols > X) {
   Matrix<_Scalar,-1,-1> MatrixOut(this->nRows,X.cols());
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < X.cols(); ++jj) {
           for (uint pp = 0; pp < this->nCols; ++pp) {
               MatrixOut.data[ii][jj] += this->data[ii][pp]*X.data[pp][jj];
           }
       }
   }
   return MatrixOut;
}

// Overloading operators
template<typename _Scalar>
_Scalar &Matrix<_Scalar,-1,-1>::operator () (int row, int col) {
   assert(col >= 0 && col < nCols);
   assert(row >= 0 && row < nRows);
   return data[row][col];
}

template<typename _Scalar>
template<typename _OtherScalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::dot(Matrix<_OtherScalar,-1,-1> X) {
   Matrix<_Scalar,-1,-1> MatrixOut(this->nRows,X.cols());
   for (uint ii = 0; ii < this->nRows; ++ii) {
       for (uint jj = 0; jj < X.cols(); ++jj) {
           for (uint pp = 0; pp < this->nCols; ++pp) {
               MatrixOut.data[ii][jj] += this->data[ii][pp]*X.data[pp][jj];
           }
       }
   }
   return MatrixOut;
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherRows>
Vector<_Scalar,-1> Matrix<_Scalar,-1,-1>::dot(Vector<_OtherScalar,_OtherRows> V) {
   Vector<_Scalar,-1> VectOut(this->rows());
   for (uint ii = 0; ii < this->rows(); ++ii) {
       for (uint pp = 0; pp < this->cols(); ++pp) {
           VectOut.data[ii] += this->data[ii][pp]*V.data[pp];
       }
   }
   return VectOut;
}

template<typename _Scalar>
template<typename _OtherScalar>
Vector<_Scalar,-1> Matrix<_Scalar,-1,-1>::dot(Vector<_OtherScalar,-1> V) {
   Vector<_Scalar,-1> VectOut(this->rows());
   for (uint ii = 0; ii < this->rows(); ++ii) {
       for (uint pp = 0; pp < this->cols(); ++pp) {
           VectOut.data[ii] += this->data[ii][pp]*V.data[pp];
       }
   }
   return VectOut;
}


template<typename _Scalar>
double Matrix<_Scalar,-1,-1>::determinant() {
   assert(this->nRows == this->nCols);
   double detVal; // determinant
   switch (this->nRows) {
   case 2:
       detVal = data[0][0]*data[1][1] - data[0][1]*data[1][0];
       return detVal;
       break;
   case 3:
       detVal = data[0][0]*(data[1][1]*data[2][2] - data[1][2]*data[2][1])
               -data[0][1]*(data[1][0]*data[2][2] - data[1][2]*data[2][0])
               +data[0][2]*(data[1][0]*data[2][1] - data[1][1]*data[2][0]);
       return detVal;
       break;
   default:
       std::cout << "Currently 2 or 3 dim Matrix are supported\n";
       return detVal;
       break;
   }
}

template<typename _Scalar>
Matrix<_Scalar,-1,-1> Matrix<_Scalar,-1,-1>::inverse() {
   assert(this->nRows == this->nCols);
   double detVal = this->determinant(); // determinant
   Matrix<_Scalar,-1,-1> invMat;
   invMat.resize(this->nRows,this->nCols);
   switch (this->nRows) {
       case 2:
           invMat.data[0][0] = +data[1][1]/detVal;
           invMat.data[0][1] = -data[0][1]/detVal;
           invMat.data[1][0] = -data[1][0]/detVal;
           invMat.data[1][1] = +data[0][0]/detVal;
           break;
       case 3:
           double a,b,c,d,e,f,g,h,i;
           a = data[0][0]; b = data[0][1]; c = data[0][2];
           d = data[1][0]; e = data[1][1]; f = data[1][2];
           g = data[2][0]; h = data[2][1]; i = data[2][2];

           invMat.data[0][0] = +(e*i - f*h)/detVal;
           invMat.data[0][1] = -(b*i - c*h)/detVal;
           invMat.data[0][2] = +(b*f - c*e)/detVal;
           invMat.data[1][0] = -(d*i - f*g)/detVal;
           invMat.data[1][1] = +(a*i - c*g)/detVal;
           invMat.data[1][2] = -(a*f - c*d)/detVal;
           invMat.data[2][0] = +(d*h - e*g)/detVal;
           invMat.data[2][1] = -(a*h - b*g)/detVal;
           invMat.data[2][2] = +(a*e - b*d)/detVal;
           break;
       default:
           std::cout << "Currently 2 or 3 dim Matrix are supported\n";
           break;
   }
   return invMat;
}
