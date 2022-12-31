//  Inertial Measurement Unit Maths Library
//
//  Copyright 2013-2021 Sam Cowen <samuel.cowen@camelsoftware.com>
//  Bug fixes and cleanups by Gé Vissers (gvissers@gmail.com)
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.

#ifndef IMUMATH_MATRIX_HPP
#define IMUMATH_MATRIX_HPP

#include <stdint.h>
#include <string.h>

#include "vector.h"

namespace imu {

template <uint8_t N> class Matrix {
public:
  Matrix() { memset(_cell_data, 0, N * N * sizeof(double)); }

  Matrix(const Matrix &m) {
    for (int ij = 0; ij < N * N; ++ij) {
      _cell_data[ij] = m._cell_data[ij];
    }
  }

  ~Matrix() {}

  Matrix &operator=(const Matrix &m) {
    for (int ij = 0; ij < N * N; ++ij) {
      _cell_data[ij] = m._cell_data[ij];
    }
    return *this;
  }

  Vector<N> row_to_vector(int i) const {
    Vector<N> ret;
    for (int j = 0; j < N; j++) {
      ret[j] = cell(i, j);
    }
    return ret;
  }

  Vector<N> col_to_vector(int j) const {
    Vector<N> ret;
    for (int i = 0; i < N; i++) {
      ret[i] = cell(i, j);
    }
    return ret;
  }

  void vector_to_row(const Vector<N> &v, int i) {
    for (int j = 0; j < N; j++) {
      cell(i, j) = v[j];
    }
  }

  void vector_to_col(const Vector<N> &v, int j) {
    for (int i = 0; i < N; i++) {
      cell(i, j) = v[i];
    }
  }

  double operator()(int i, int j) const { return cell(i, j); }
  double &operator()(int i, int j) { return cell(i, j); }

  double cell(int i, int j) const { return _cell_data[i * N + j]; }
  double &cell(int i, int j) { return _cell_data[i * N + j]; }

  Matrix operator+(const Matrix &m) const {
    Matrix ret;
    for (int ij = 0; ij < N * N; ++ij) {
      ret._cell_data[ij] = _cell_data[ij] + m._cell_data[ij];
    }
    return ret;
  }

  Matrix operator-(const Matrix &m) const {
    Matrix ret;
    for (int ij = 0; ij < N * N; ++ij) {
      ret._cell_data[ij] = _cell_data[ij] - m._cell_data[ij];
    }
    return ret;
  }

  Matrix operator*(double scalar) const {
    Matrix ret;
    for (int ij = 0; ij < N * N; ++ij) {
      ret._cell_data[ij] = _cell_data[ij] * scalar;
    }
    return ret;
  }

  Matrix operator*(const Matrix &m) const {
    Matrix ret;
    for (int i = 0; i < N; i++) {
      Vector<N> row = row_to_vector(i);
      for (int j = 0; j < N; j++) {
        ret(i, j) = row.dot(m.col_to_vector(j));
      }
    }
    return ret;
  }

  Vector<4> operator*(const Vector<4> &v) const {
    
    float w = v.quat_w();
    float x = v.quat_x();
    float y = v.quat_y();
    float z = v.quat_z();

    float ret_w = cell(0,0) * w + cell(0,1) * y + cell(0,2) * z + cell(0,3) * z;
    float ret_x = cell(1,0) * w + cell(1,1) * y + cell(1,2) * z + cell(1,3) * z;
    float ret_y = cell(2,0) * w + cell(2,1) * y + cell(2,2) * z + cell(2,3) * z;
    float ret_z = cell(3,0) * w + cell(3,1) * y + cell(3,2) * z + cell(3,3) * z;
    
    return Vector<4>(ret_w, ret_x, ret_y, ret_z);
  }

  Matrix transpose() const {
    Matrix ret;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        ret(j, i) = cell(i, j);
      }
    }
    return ret;
  }

  Matrix<N - 1> minor_matrix(int row, int col) const {
    Matrix<N - 1> ret;
    for (int i = 0, im = 0; i < N; i++) {
      if (i == row)
        continue;

      for (int j = 0, jm = 0; j < N; j++) {
        if (j != col) {
          ret(im, jm++) = cell(i, j);
        }
      }
      im++;
    }
    return ret;
  }

  double determinant() const {
    // specialization for N == 1 given below this class
    double det = 0.0, sign = 1.0;
    for (int i = 0; i < N; ++i, sign = -sign)
      det += sign * cell(0, i) * minor_matrix(0, i).determinant();
    return det;
  }

  Matrix invert() const {
    Matrix ret;
    double det = determinant();

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        ret(i, j) = minor_matrix(j, i).determinant() / det;
        if ((i + j) % 2 == 1)
          ret(i, j) = -ret(i, j);
      }
    }
    return ret;
  }

  double trace() const {
    double tr = 0.0;
    for (int i = 0; i < N; ++i)
      tr += cell(i, i);
    return tr;
  }
  
  Matrix<4> createOmegaOperator(imu::Vector<3> w){
    Matrix<4> mat;
    mat.cell(0,0) = 0;
    mat.cell(0,1) = -w.x();
    mat.cell(0,2) = -w.y();
    mat.cell(0,3) = -w.z();

    mat.cell(1,0) = w.x();
    mat.cell(1,1) = 0;
    mat.cell(1,2) = w.z();
    mat.cell(1,3) = -w.y();
    
    mat.cell(2,0) = w.y();
    mat.cell(2,1) = -w.z();
    mat.cell(2,2) = 0;
    mat.cell(2,3) = w.x();

    mat.cell(3,0) = w.z();
    mat.cell(3,1) = w.y();
    mat.cell(3,2) = -w.x();
    mat.cell(3,3) = 0;

    return mat;
  }

  

  Matrix<4> createI4Matrix(){
    Matrix<4> I4;
    for(int i=0; i<4; i++ ){
      for(int j=0; j<4; j++ ){
        if(i==j){
          I4.cell(i,j) = 1;
        }else{
          I4.cell(i,j) = 0;
        }
      }
    }
    return I4;
  }

private:
  double _cell_data[N * N];
};

template <> inline double Matrix<1>::determinant() const { return cell(0, 0); }

}; // namespace imu

#endif
