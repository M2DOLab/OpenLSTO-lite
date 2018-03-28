template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::operator = (const Vector<_Scalar,_nElem> &rhs) { // Vector Assignment; Static-Static
   assert(_nElem == rhs.data.size());
   for (uint ii = 0; ii < _nElem; ++ii) {
       data[ii] = rhs.data[ii];
   }
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::operator = (const Vector<_Scalar,-1> &rhs) { // Vector Assignment; Static-Dynamic
   assert(_nElem == rhs.data.size());
   for (uint ii = 0; ii < _nElem; ++ii) {
       data[ii] = rhs.data[ii];
   }
}

template<typename _Scalar, int _nElem>
_Scalar &Vector<_Scalar,_nElem>::operator () (const double indx) { // Value Retrieval
   if (indx > _nElem-1 || indx != round(indx) || indx < 0) {
       std::cout << "\nERROR: The Vector index should be an Integer from 0 to _nElem, inclusive.\n\n";
       assert(false);
   }
   else {
       return data[indx];
   }
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::operator += (const double alpha) { // Constant Addition
   uint flag = 0,ii;

   for (ii = 0; ii < _nElem; ++ii) {
       if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
           flag++;
       }
       data[ii] += alpha;
   }
   if (flag > 0) {
       std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of addition!\n";
   }
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::operator -= (const double alpha) { // Constant Subtraction
   uint flag = 0,ii;
   for (ii = 0; ii < _nElem; ++ii) {
       if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
           flag++;
       }
       data[ii] -= alpha;
   }
   if (flag > 0) {
       std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of subtraction!\n";
   }
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::operator *= (const double alpha) { // Constant Multiplication
   uint flag = 0,ii;
   for (ii = 0; ii < _nElem; ++ii) {
       if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
           flag++;
       }
       data[ii] *= alpha;
   }
   if (flag > 0) {
       std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of multiplication!\n";
   }
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::operator /= (const double alpha) { // Constant Division
   if (alpha == 0) {
       std::cout << "ERROR: The divisor cannot be zero!\n\n";
       assert(false);
   }
   else {
       uint flag = 0,ii;
       for (ii = 0; ii < _nElem; ++ii) {
           if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
               flag++;
           }
           data[ii] /= alpha;
       }
       if (flag > 0) {
           std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of division!\n";
       }
   }
}

// binary
template<typename _Scalar, int _nElem>
Vector<_Scalar,_nElem> Vector<_Scalar,_nElem>::operator + (const Vector<_Scalar,_nElem> &rhs) { // Vector Addition; Static-Static
   Vector<_Scalar,_nElem> vect;
   for (uint ii = 0; ii < _nElem; ++ii) {
       vect[ii] = this->data[ii] + rhs.data[ii];
   }
   return vect;
}

template<typename _Scalar, int _nElem>
Vector<_Scalar,_nElem> Vector<_Scalar,_nElem>::operator - (const Vector<_Scalar,_nElem> &rhs) { // Vector Subtraction; Static-Static
   Vector<_Scalar,_nElem> vect;
   for (uint ii = 0; ii < _nElem; ++ii) {
       vect[ii] = this->data[ii] - rhs.data[ii];
   }
   return vect;
}

template<typename _Scalar, int _nElem>
Vector<_Scalar,_nElem> Vector<_Scalar,_nElem>::operator * (const Vector<_Scalar,_nElem> &rhs) { // Vector Multiplication; Static-Static
   Vector<_Scalar,_nElem> vect;
   for (uint ii = 0; ii < _nElem; ++ii) {
       vect[ii] = this->data[ii] * rhs.data[ii];
   }
   return vect;
}

// template<typename _Scalar, int _nElem>
// Vector<_Scalar,_nElem> &Vector<_Scalar,_nElem>::operator /= (const Vector<_Scalar,_nElem> &rhs) { // Vector Division; Static-Static
//     uint flag = 0,jj = 0;
//     while (flag == 0 && jj < _nElem) {
//         if (rhs.data[jj] == 0) {
//             flag = 1;
//         }
//         jj++;
//     }
//     if (flag == 1) {
//         std::cout << "ERROR: The divisor vector contains a zero value element!\n\n";
//         assert(false);
//     }
//     else {
//         for (uint ii = 0; ii < _nElem; ++ii) {
//             data[ii] /= rhs.data[ii];
//         }
//     return data;
//     }
// }

template<typename _Scalar, int _nElem>
Vector<_Scalar,-1> Vector<_Scalar,_nElem>::operator + (const Vector<_Scalar,-1> &rhs) { // Vector Addition; Static-Dynamic
   Vector<_Scalar,-1> vect(_nElem);
   for (uint ii = 0; ii < _nElem; ++ii) {
       vect[ii] = data[ii] + rhs.data[ii];
   }
   return vect;
}

template<typename _Scalar, int _nElem>
Vector<_Scalar,-1> Vector<_Scalar,_nElem>::operator - (const Vector<_Scalar,-1> &rhs) { // Vector Subtraction; Static-Dynamic
   Vector<_Scalar,-1> vect(_nElem);
   for (uint ii = 0; ii < _nElem; ++ii) {
       vect[ii] = data[ii] - rhs.data[ii];
   }
   return vect;
}

template<typename _Scalar, int _nElem>
Vector<_Scalar,-1> Vector<_Scalar,_nElem>::operator * (const Vector<_Scalar,-1> &rhs) { // Vector Multiplication; Static-Dynamic
   Vector<_Scalar,-1> vect(_nElem);
   for (uint ii = 0; ii < _nElem; ++ii) {
       vect[ii] = data[ii] * rhs.data[ii];
   }
   return vect;
}

// template<typename _Scalar, int _nElem>
// Vector<_Scalar,_nElem> &Vector<_Scalar,_nElem>::operator / (const Vector<_Scalar,-1> &rhs) { // Vector Division; Static-Dynamic
//     uint flag = 0,jj = 0;
//     while (flag == 0 && jj < _nElem) {
//         if (rhs.data[jj] == 0) {
//             flag = 1;
//         }
//         jj++;
//     }
//     if (flag == 1) {
//         std::cout << "ERROR: The divisor vector contains a zero value element!\n\n";
//         assert(false);
//     }
//     else {
//         for (uint ii = 0; ii < _nElem; ++ii) {
//             data[ii] /= rhs.data[ii];
//         }
//     return data;
//     }
// }

template<typename _Scalar, int _nElem>
bool Vector<_Scalar,_nElem>::operator == (const Vector<_Scalar,_nElem> &rhs) { // Equality Check; Static-Static
   if (_nElem != rhs.data.size()) {
       return false;
   }
   else {
       uint ii = 0;
       bool equality = true;
       while (ii < _nElem) {
           if (data[ii] != rhs.data[ii]) {
               ii = _nElem;
               equality = false;
           }
           ii++;
       }
       return equality;
   }
}

template<typename _Scalar, int _nElem>
bool Vector<_Scalar,_nElem>::operator == (const Vector<_Scalar,-1> &rhs) { // Equality Check; Static-Dynamic
   if (_nElem != rhs.data.size()) {
       return false;
   }
   else {
       uint ii = 0;
       bool equality = true;
       while (ii < _nElem) {
           if (data[ii] != rhs.data[ii]) {
               ii = _nElem;
               equality = false;
           }
           ii++;
       }
       return equality;
   }
}

template<typename _Scalar, int _nElem>
uint Vector<_Scalar,_nElem>::size() const{
   return _nElem;
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::fill(_Scalar _value) {
   for (uint ii = 0; ii < _nElem; ++ii) {
       data[ii] = _value;
   }
}

template<typename _Scalar, int _nElem>
double Vector<_Scalar,_nElem>::norm() const {
   double length_ = 0.0;
   for (uint ii = 0; ii < _nElem; ++ii) {
       length_ += data[ii]*data[ii];
   }
   length_ = std::sqrt(length_);
   return length_;
}

template<typename _Scalar, int _nElem>
void Vector<_Scalar,_nElem>::print() const{
   std::cout << "\n [ ";
   for (uint ii = 0; ii < _nElem-1; ++ii) {
       std::cout << data[ii] << " , ";
   }
   std::cout << data[_nElem-1] << " ]\n";
}

template<typename _Scalar, int _nElem>
double Vector<_Scalar,_nElem>::dot(Vector<_Scalar,_nElem> &rhs) { // Static-Static
   double dot_ = 0.0;
   for (uint ii = 0; ii < _nElem; ++ii) {
       dot_ += data[ii]*rhs(ii);
   }
   return dot_;
}

template<typename _Scalar, int _nElem>
double Vector<_Scalar,_nElem>::dot(Vector<_Scalar,-1> &rhs) { // Static-Dynamic
   double dot_ = 0.0;
   for (uint ii = 0; ii < _nElem; ++ii) {
       dot_ += data[ii]*rhs(ii);
   }
   return dot_;
}
template<typename _Scalar>
uint Vector<_Scalar,-1>::size() const{
   return data.size();
}

template<typename _Scalar>
void Vector<_Scalar,-1>::fill(_Scalar _value) {
   for (uint ii = 0; ii < this->size(); ++ii) {
       data[ii] = _value;
   }
}

template<typename _Scalar>
double Vector<_Scalar,-1>::norm() const {
   double length_ = 0.0;
   for (uint ii = 0; ii < data.size(); ++ii) {
       length_ += data[ii]*data[ii];
   }
   length_ = std::sqrt(length_);
   return length_;
}

template<typename _Scalar>
void Vector<_Scalar,-1>::resize(int nx) {
   data.clear();
   data.resize(nx);
}

template<typename _Scalar>
void Vector<_Scalar,-1>::print() const {
   std::cout << "\n [ ";
   for (uint ii = 0; ii < data.size()-1; ++ii) {
       std::cout << data[ii] << " , ";
   }
   std::cout << data[size()-1] << " ]\n";
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherElem>
double Vector<_Scalar,-1>::dot(Vector<_OtherScalar,_OtherElem> rhs) {
   double dot_ = 0.0;
   for (uint ii = 0; ii < data.size(); ++ii) {
       dot_ += data[ii]*rhs.data[ii];
   }
   return dot_;
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherElem>
void Vector<_Scalar,-1>::operator = (Vector<_OtherScalar,_OtherElem> rhs) {
   data.resize(_OtherElem);
   for (uint ii = 0; ii < _OtherElem; ++ii) {
       data[ii] = rhs.data[ii];
   }
}

template<typename _Scalar>
template<typename _OtherScalar>
void Vector<_Scalar,-1>::operator = (Vector<_OtherScalar,-1> rhs) {
   data.resize(rhs.size());
   for (uint ii = 0; ii < rhs.size(); ++ii) {
       data[ii] = rhs.data[ii];
   }
}

template<typename _Scalar>
_Scalar &Vector<_Scalar,-1>::operator () (const double indx) { // Value Retrieval
   if (indx > data.size()-1 || indx != round(indx) || indx < 0) {
       std::cout << "\nERROR: The Vector index should be an Integer from 0 to _nElem, inclusive.\n\n";
       assert(false);
   }
   else {
       return data[indx];
   }
}

template<typename _Scalar>
void Vector<_Scalar,-1>::operator += (const double alpha) { // Constant Addition
   uint flag = 0,ii;
   for (ii = 0; ii < data.size(); ++ii) {
       if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
           flag++;
       }
       data[ii] += alpha;
   }
   if (flag > 0) {
       std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of addition!\n";
   }
}

template<typename _Scalar>
void Vector<_Scalar,-1>::operator -= (const double alpha) { // Constant Subtraction
   uint flag = 0,ii;
   for (ii = 0; ii < data.size(); ++ii) {
       if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
           flag++;
       }
       data[ii] -= alpha;
   }
   if (flag > 0) {
       std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of subtraction!\n";
   }
}

template<typename _Scalar>
void Vector<_Scalar,-1>::operator *= (const double alpha) { // Constant Multiplication
   uint flag = 0,ii;
   for (ii = 0; ii < data.size(); ++ii) {
       if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
           flag++;
       }
       data[ii] *= alpha;
   }
   if (flag > 0) {
       std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of multiplication!\n";
   }
}

template<typename _Scalar>
void Vector<_Scalar,-1>::operator /= (const double alpha) { // Constant Division
   if (alpha == 0) {
       std::cout << "ERROR: The divisor cannot be zero!\n\n";
       assert(false);
   }
   else {
       uint flag = 0,ii;
       for (ii = 0; ii < data.size(); ++ii) {
           if (data[ii] == round(data[ii]) && alpha != round(alpha)) {
               flag++;
           }
           data[ii] /= alpha;
       }
       if (flag > 0) {
           std::cout << "WARNING: The final vector may have been rounded due to the DOUBLE constant of division!\n";
       }
   }
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherElem>
Vector<_Scalar,-1> Vector<_Scalar,-1>::operator + (const Vector<_OtherScalar,_OtherElem> &rhs) { // Vector Addition
   Vector<_Scalar,-1> vect(this->size());
   for (uint ii = 0; ii < data.size(); ++ii) {
       vect[ii] = this->data[ii] + rhs.data[ii];
   }
   return vect;
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherElem>
Vector<_Scalar,-1> Vector<_Scalar,-1>::operator - (const Vector<_OtherScalar,_OtherElem> &rhs) { // Vector Subtraction
   Vector<_Scalar,-1> vect(this->size());
   for (uint ii = 0; ii < data.size(); ++ii) {
       vect[ii] = this->data[ii] - rhs.data[ii];
   }
   return vect;
}

template<typename _Scalar>
template<typename _OtherScalar, int _OtherElem>
Vector<_Scalar,-1> Vector<_Scalar,-1>::operator * (const Vector<_OtherScalar,_OtherElem> &rhs) { // Vector Multiplication
   Vector<_Scalar,-1> vect(this->size());
   for (uint ii = 0; ii < data.size(); ++ii) {
       vect[ii] = this->data[ii] * rhs.data[ii];
   }
   return vect;
}

// template<typename _Scalar>
// template<typename _OtherScalar, int _OtherElem>
// Vector<_Scalar,-1> &Vector<_Scalar,-1>::operator /= (const Vector<_OtherScalar,_OtherElem> &rhs) { // Vector Division
//     uint flag = 0,jj = 0;
//     while (flag == 0 && jj < data.size()) {
//         if (rhs.data[jj] == 0) {
//             flag = 1;
//         }
//         jj++;
//     }
//     if (flag == 1) {
//         std::cout << "ERROR: The divisor vector contains a zero value element!\n\n";
//         assert(false);
//     }
//     else {
//         for (uint ii = 0; ii < data.size(); ++ii) {
//             data[ii] /= rhs.data[ii];
//         }
//     return data;
//     }
// }

template<typename _Scalar>
template<typename _OtherScalar, int _OtherElem>
bool Vector<_Scalar,-1>::operator == (const Vector<_OtherScalar,_OtherElem> &rhs) { // Equality Check
   if (data.size() != rhs.data.size()) {
       return false;
   }
   else {
       uint ii = 0;
       bool equality = true;
       while (ii < data.size()) {
           if (data[ii] != rhs.data[ii]) {
               ii = data.size();
               equality = false;
           }
           ii++;
       }
       return equality;
   }
}
