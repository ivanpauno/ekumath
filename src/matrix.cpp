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
}
