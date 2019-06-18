#include <cmath>
#include <sstream>
#include <stdexcept>

#include "matrix.hpp"

namespace ekumath {
  Matrix::Matrix(size_t rows, size_t columns, double value) :
    rows_(rows), columns_(columns), data_(columns)
  {
    if (!rows || !columns) {
      throw std::runtime_error("Can't create matrix of zero rows or columns");
    }
    for (auto & column: data_) {
      column.resize(rows, value);
    }
  }

  Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list)
  {
    rows_ = list.size();
    columns_ = list.begin()->size();
    for (auto it = list.begin() + 1; it != list.end(); it++) {
      if (it->size() != columns_) {
        throw std::runtime_error("All the rows should have the same number of elements");
      }
    }

    data_.resize(columns_);
    for (auto & column: data_) {
      column.resize(rows_);
    }

    size_t row = 0;
    for (auto & subl: list) {
      size_t col = 0;
      for (double elem: subl) {
        this->unsafe_at(row, col) = elem;
        col++;
      }
      row++;
    }
  }

  Matrix::Matrix(std::initializer_list<double> list)
  {
    rows_ = list.size();
    columns_ = 1U;
    data_.resize(1);
    data_[0] = list;
  }

  Matrix
  Matrix::Identity(size_t size)
  {
    Matrix A(size, size);
    for (size_t i; i < size; i++) {
      A.unsafe_at(i, i) = 1.;
    }
    return A;
  }

  Matrix
  Matrix::Random(size_t rows, size_t columns)
  {
    Matrix A(rows, columns);
    for (auto & column: A.data_) {
      for (auto & element: column) {
        element = static_cast<double>(std::rand()) + static_cast<double>(std::rand()) / 1000;
      }
    }
    return A;
  }

  double &
  Matrix::operator()(size_t row, size_t column)
  {
    if (row >= rows_ || column >= columns_) {
      std::ostringstream ss;
      ss << "Tried to access member (" << row << ", " << column <<
        "\n Matrix size was (" << rows_ << ", " << columns_ << ")";
      throw std::out_of_range(ss.str().c_str());
    }
    return unsafe_at(row, column);
  }

  double
  Matrix::operator()(size_t row, size_t column) const
  {
    return const_cast<Matrix *>(this)->operator()(row, column);
  }

  double &
  Matrix::unsafe_at(size_t row, size_t column)
  {
    return data_[column][row];
  }

  double
  Matrix::unsafe_at(size_t row, size_t column) const
  {
    return const_cast<Matrix *>(this)->unsafe_at(row, column);
  }

  size_t
  Matrix::rows() const
  {
    return rows_;
  }

  size_t
  Matrix::cols() const
  {
    return columns_;
  }

  size_t
  Matrix::size() const
  {
    return data_.size();
  }

  bool
  operator==(const Matrix & lhs, const Matrix & rhs)
  {
    if (lhs.rows_ != rhs.rows_ || lhs.columns_ != rhs.columns_) {
      return false;
    }
    return lhs.data_ == rhs.data_;
  }

  std::ostream &
  operator<<(std::ostream &stream, const Matrix & A)
  {
    size_t row = 0U;

    if (A.cols() > 1) {
      stream << "[";
      for (; row < (A.rows()-1U); row++) {
        A.print_row(stream, row);
        stream << ", ";
      }
      A.print_row(stream, row);
      stream << "]";
    } else {
      A.print_col(stream, 0);
    }
    return stream;
  }

  void
  Matrix::print_row(std::ostream & stream, size_t row) const
  {
    size_t col = 0;

    stream << "[";
    for (; col < (columns_ - 1U); col++) {
      stream << unsafe_at(row, col) << ", ";
    }
    stream << unsafe_at(row, col);
    stream << "]";
  }

  void
  Matrix::print_col(std::ostream & stream, size_t col) const
  {
    size_t row = 0;

    stream << "[";
    for (; row < (rows_ - 1U); row++) {
      stream << unsafe_at(row, col) << ", ";
    }
    stream << unsafe_at(row, col);
    stream << "]";
  }

#define ELEMENT_WISE_OPERATION_WITH_DOUBLE(A, x, op) \
  __ELEMENT_WISE_OPERATION_WITH_DOUBLE((A), (x), op)
#define __ELEMENT_WISE_OPERATION_WITH_DOUBLE(A, x, op) \
  { \
    Matrix result(A.rows(), A.cols()); \
    for (size_t col; col < A.cols(); col++) { \
      for (size_t row; row < A.rows(); row++) { \
        result.unsafe_at(row, col) = A.unsafe_at(row, col) op x; \
      } \
    } \
    return result; \
  }

  Matrix
  operator+(double x, const Matrix & A)
  {
    return A + x;
  }

  Matrix
  operator+(const Matrix & A, double x)
  ELEMENT_WISE_OPERATION_WITH_DOUBLE(A, x, +)

  Matrix
  operator*(double x, const Matrix & A)
  {
    return A * x;
  }

  Matrix
  operator*(const Matrix & A, double x)
  ELEMENT_WISE_OPERATION_WITH_DOUBLE(A, x, *)

  Matrix
  operator/(const Matrix & A, double x)
  ELEMENT_WISE_OPERATION_WITH_DOUBLE(A, x, /)

#define ELEMENT_WISE_OPERATION_WITH_ANOTHER_MATRIX(A, B, op) \
  __ELEMENT_WISE_OPERATION_WITH_ANOTHER_MATRIX((A), (B), op)
#define __ELEMENT_WISE_OPERATION_WITH_ANOTHER_MATRIX(A, B, op) \
  { \
    if (A.rows() != B.rows() || A.cols() != B.cols()) { \
      throw std::runtime_error("Matrices must have the same dimension"); \
    } \
    Matrix result(A.rows(), A.cols()); \
    for (size_t col; col < A.cols(); col++) { \
      for (size_t row; row < A.rows(); row++) { \
        result.unsafe_at(row, col) = A.unsafe_at(row, col) op B.unsafe_at(row, col); \
      } \
    } \
    return result; \
  }

  Matrix
  operator+(const Matrix & A, const Matrix & B)
  ELEMENT_WISE_OPERATION_WITH_ANOTHER_MATRIX(A, B, +)

  Matrix
  operator*(const Matrix & A, const Matrix & B)
  ELEMENT_WISE_OPERATION_WITH_ANOTHER_MATRIX(A, B, *)

  Matrix
  operator^(const Matrix & A, double x)
  {
    Matrix result(A.rows(), A.cols());
    for (size_t col; col < A.cols(); col++) {
      for (size_t row; row < A.rows(); row++) {
        result.unsafe_at(row, col) = std::pow(A.unsafe_at(row, col), x);
      }
    }
    return result;
  }

  Matrix
  pow(const Matrix & A, double x)
  {
    return A ^ x;
  }

  Matrix
  sqrt(const Matrix & A)
  {
    return A ^ .5;
  }

  Matrix
  Matrix::get_transpose() const
  {
    Matrix T(this->rows(), this->cols());
    for (size_t col; col < this->cols(); col++) {
      for (size_t row; row < this->rows(); row++) {
        T.unsafe_at(row, col) = this->unsafe_at(col, row);
      }
    }
    return T;
  }

  Matrix &
  Matrix::transpose()
  {
    *this = get_transpose();
    return *this;
  }

  Matrix &
  Matrix::concatenate(const Matrix & B) {
    if (B.rows() != this->rows()) {
      throw std::runtime_error("Trying to concatenate matrices with different number of rows");
    }
    size_t actual_cols = this->cols();
    size_t extra_cols = B.cols();
    this->data_.resize(actual_cols + extra_cols);
    for (size_t col = 0; col < extra_cols; col++) {
      this->data_[actual_cols + col] = B.data_[col];
    }
    return *this;
  }

  Matrix
  operator|(const Matrix & A, const Matrix & B)
  {
    Matrix result = A;
    return result.concatenate(B);
  }
}
