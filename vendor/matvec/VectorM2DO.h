#ifndef VECTOR_M2DO_H
#define VECTOR_M2DO_H

#include <cmath>
#include <vector>
#include <assert.h>
#include <iostream>
typedef unsigned int uint;

using namespace std;

template<typename _Scalar, const int _nElem>
class Vector { // Static-sized

   public:
   // Constructor
   Vector() {
       data.resize(_nElem);
       this->fill(0.0); // Initializes "vect" into _nElem elements of zeros
   };

   // Class Function Initializations
   uint size() const;
   void fill(_Scalar x = 0.0);
   double norm() const;
   void print() const;

   double dot(Vector<_Scalar,_nElem> &rhs); // Static-Static
   double dot(Vector<_Scalar,-1> &rhs); // Static-Dynamic

   // Operator Override Initializations
   void operator = (const Vector<_Scalar,_nElem> &rhs); // Vector Assignment; Static-Static
   void operator = (const Vector<_Scalar,-1> &rhs); // Vector Assignment; Static-Dynamic

   _Scalar &operator () (const double indx); // Value Retrieval

   void operator += (const double alpha); // Constant Addition
   void operator -= (const double alpha); // Constant Subtraction
   void operator *= (const double alpha); // Constant Multiplication
   void operator /= (const double alpha); // Constant Division

   Vector<_Scalar,_nElem> operator + (const Vector<_Scalar,_nElem> &rhs); // Vector Addition; Static-Static
   Vector<_Scalar,_nElem> operator - (const Vector<_Scalar,_nElem> &rhs); // Vector Subtraction; Static-Static
   Vector<_Scalar,_nElem> operator * (const Vector<_Scalar,_nElem> &rhs); // Vector Multiplication; Static-Static
   // Vector<_Scalar,_nElem> operator / (const Vector<_Scalar,_nElem> &rhs); // Vector Division; Static-Static

   Vector<_Scalar,-1> operator + (const Vector<_Scalar,-1> &rhs); // Vector Addition; Static-Dynamic
   Vector<_Scalar,-1> operator - (const Vector<_Scalar,-1> &rhs); // Vector Subtraction; Static-Dynamic
   Vector<_Scalar,-1> operator * (const Vector<_Scalar,-1> &rhs); // Vector Multiplication; Static-Dynamic
   // Vector<_Scalar,-1> operator / (const Vector<_Scalar,-1> &rhs); // Vector Division; Static-Dynamic

   bool operator == (const Vector<_Scalar,_nElem> &rhs); // Equality Check; Static-Static
   bool operator == (const Vector<_Scalar,-1> &rhs); // Equality Check; Static-Dynamic

   // Variable Declaration
   std::vector<_Scalar> data;
};


template<typename _Scalar>
class Vector<_Scalar,-1> { // Dynamic-sized
   public:
   Vector()
   {
     resize(0);
   };
   Vector(int nx){
       resize(nx);
       fill(0.0);
   }
   // Class Function Initializations
   void resize(int);
   void fill(_Scalar);
   uint size() const;
   double norm() const;
   void print() const;

   template<typename _OtherScalar, int _OtherElem>
   double dot(const Vector<_OtherScalar,_OtherElem> rhs);

   // Operator Override Initializations
   template<typename _OtherScalar, int _OtherElem>
   void operator = (Vector<_OtherScalar,_OtherElem> rhs); // Vector Assignment
   template<typename _OtherScalar>
   void operator = (Vector<_OtherScalar,-1> rhs); // Vector Assignment


   _Scalar &operator () (const double indx); // Value Retrieval

   void operator += (const double alpha); // Constant Addition
   void operator -= (const double alpha); // Constant Subtraction
   void operator *= (const double alpha); // Constant Multiplication
   void operator /= (const double alpha); // Constant Division

   template<typename _OtherScalar, int _OtherElem>
   Vector<_Scalar,-1> operator + (const Vector<_OtherScalar,_OtherElem> &rhs); // Vector Addition
   template<typename _OtherScalar, int _OtherElem>
   Vector<_Scalar,-1> operator - (const Vector<_OtherScalar,_OtherElem> &rhs); // Vector Subtraction
   template<typename _OtherScalar, int _OtherElem>
   Vector<_Scalar,-1> operator * (const Vector<_OtherScalar,_OtherElem> &rhs); // Vector Multiplication
   // template<typename _OtherScalar, int _OtherElem>
   // Vector<_Scalar,-1> operator / (const Vector<_OtherScalar,_OtherElem> &rhs); // Vector Division

   Vector<_Scalar,-1> operator + (const Vector<_Scalar,-1> &rhs); // Vector Addition
   Vector<_Scalar,-1> operator - (const Vector<_Scalar,-1> &rhs); // Vector Subtraction
   Vector<_Scalar,-1> operator * (const Vector<_Scalar,-1> &rhs); // Vector Multiplication

   template<typename _OtherScalar, int _OtherElem>
   bool operator == (const Vector<_OtherScalar,_OtherElem> &rhs); // Equality Check

   // Variable Declaration
   std::vector<_Scalar> data;
};
#include "VectorM2DO.cpp"
#endif
