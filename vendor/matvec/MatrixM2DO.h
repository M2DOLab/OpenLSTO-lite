#ifndef MATRIX_M2DO_H
#define MATRIX_M2DO_H

#include <iostream>
#include <cmath>
#include <vector>
#include <assert.h>
#include "VectorM2DO.h"
typedef unsigned int uint;

template<typename _Scalar, int _Rows, int _Cols>
class Matrix { // Static-sized

   public:
   Matrix() {
       nRows = _Rows; nCols = _Cols;
       data.resize(nRows,std::vector<_Scalar> (nCols,0.0));
       this->fill(0.0);
   };

   // unary member functions
   void fill(const _Scalar &val); // fill with val
   uint rows() const; // get nRows
   uint cols() const; // get nCols
   uint size() const; // nRows * nCols
   void print() const;

   Matrix<_Scalar,_Cols,_Rows> transpose(); // transpose of the self
   double determinant(); // determinant, 2x2, 3x3 only
   Matrix<_Scalar,_Rows,_Cols> inverse(); // 3x3, 2x2 only

   // operator override
   _Scalar &operator () (int row, int col); // element Accessor
   const bool operator == (const Matrix &rhs); // check equality
   void operator = (const Matrix &rhs); // bulk Assignment

   void setIdentity();
   static Matrix<_Scalar,_Rows,_Cols> Zero();

   template<typename _OtherScalar, int _OtherRows, int _OtherCols> // dot with static
   Matrix<_Scalar,_Rows,_OtherCols> dot(Matrix<_OtherScalar,_OtherRows,_OtherCols> X);

   template<typename _OtherScalar> // dot product with dynamic
   Matrix<_Scalar,-1,-1> dot(Matrix<_OtherScalar,-1,-1> X);

   template<typename _OtherScalar, int _OtherRows>
   Vector<_Scalar,_Rows> dot(Vector<_OtherScalar,_OtherRows> V);

   template<typename _OtherScalar>
   Vector<_Scalar,-1> dot(Vector<_OtherScalar,-1> V);

   // Operator override for elem-wise
   void operator += (const _Scalar &val); // element-wise addition
   void operator -= (const _Scalar &val); // element-wise subtraction
   void operator *= (const _Scalar &val); // scalar multiplication

   Matrix<_Scalar,_Rows,_Cols> operator + (Matrix<_Scalar,_Rows,_Cols> rhs);
   Matrix<_Scalar,_Rows,_Cols> operator - (Matrix<_Scalar,_Rows,_Cols> rhs);

   template<typename _OtherScalar>
   Matrix<_Scalar,-1,-1> operator + (Matrix<_OtherScalar,-1,-1> rhs);
   template<typename _OtherScalar>
   Matrix<_Scalar,-1,-1> operator - (Matrix<_OtherScalar,-1,-1> rhs);

   // protected:
   std::vector<std::vector<_Scalar> > data;
   private:
   int nRows,nCols;
};


template<typename _Scalar>
class Matrix<_Scalar,-1,-1> { // Dynamic-sized
   public:
   Matrix() {
       nRows = 0; nCols = 0;
   };
   Matrix(int nx, int ny) {
       this->resize(nx,ny);
   }

   // unary member functions
   void resize(int _nRows, int _nCols); // resize matrix
   void fill(_Scalar);
   uint rows() const; // get nRows
   uint cols() const; // get nCols
   uint size() const; // nRows * nCols
   void print() const;

   Matrix<_Scalar,-1,-1> transpose(); // transpose of the self
   double determinant(); // determinant, 2x2, 3x3 only
   Matrix<_Scalar,-1,-1> inverse(); // 3x3, 2x2 only

   void setIdentity();

   // accessor
   _Scalar &operator () (int row, int col); // element Accessor
   void operator = (const Matrix<_Scalar,-1,-1> &rhs); // assignment

   // binary operator (dot)
   template<typename _OtherScalar, int _OtherRows, int _OtherCols>
   const bool operator == (const Matrix &rhs);

   template<typename _OtherScalar, int _OtherRows, int _OtherCols> // dot product with static matrix
   Matrix<_Scalar,-1,-1> dot(Matrix<_OtherScalar,_OtherRows,_OtherCols> X);

   template<typename _OtherScalar> // dot product with dynamic matrix
   Matrix<_Scalar,-1,-1> dot(Matrix<_OtherScalar,-1,-1> X);

   template<typename _OtherScalar, int _OtherRows>
   Vector<_Scalar,-1> dot(Vector<_OtherScalar,_OtherRows> V);

   template<typename _OtherScalar>
   Vector<_Scalar,-1> dot(Vector<_OtherScalar,-1> V);

   // operator override
   void operator += (const _Scalar &val); // element-wise addition
   void operator -= (const _Scalar &val); // element-wise subtraction
   void operator *= (const _Scalar &val); // scalar multiplication

   template<typename _OtherScalar, int _OtherRows, int _OtherCols>
   Matrix<_Scalar,-1,-1> operator + (Matrix<_OtherScalar, _OtherRows, _OtherCols> rhs);
   template<typename _OtherScalar, int _OtherRows, int _OtherCols>
   Matrix<_Scalar,-1,-1> operator - (Matrix<_OtherScalar, _OtherRows, _OtherCols> rhs);

   template<typename _OtherScalar>
   Matrix<_Scalar,-1,-1> operator + (Matrix<_OtherScalar,-1,-1> rhs);
   template<typename _OtherScalar>
   Matrix<_Scalar,-1,-1> operator - (Matrix<_OtherScalar,-1,-1> rhs);

   // private:
   typedef std::vector<std::vector<_Scalar> > Matrix2D;
   Matrix2D data;
   uint nRows,nCols;
};
#include "MatrixM2DO.cpp"
#endif
