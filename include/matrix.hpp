#include <initializer_list>
#include <vector>

namespace ekumath {

  class Matrix {
  friend
  std::ostream &
  operator<<(std::ostream &stream, const Matrix & A);

  friend
  bool
  operator==(const Matrix & lhs, const Matrix & rhs);

  public:
    /**
      * Constructs a `rows` x `columns` matrix.
      * All the elements are initialized to `value`.
      */
    Matrix(size_t rows, size_t columns, double value=0.);
    /**
      * Constructs a matrix from a nested initializer list.
      * e.g.: {{1.0, 2.0, 3.0}, {3.0, 1.0, 2.0}}
      *   has two rows and three columns.
      */
    Matrix(std::initializer_list<std::initializer_list<double>> list);
    /**
      * Constructs a column vector from an initializer list.
      * e.g.: {1.0, 2.0, 3.0}
      *   Is a 1x3 matrix.
      */
    Matrix(std::initializer_list<double> list);

    /* Get element in the specified row and column. */
    double &
    operator()(size_t row, size_t column);

    /* Get element in the specified row and column. */
    double
    operator()(size_t row, size_t column) const;

    /* Get number of rows. */
    size_t
    rows() const;

    /* Get number of columns. */
    size_t
    cols() const;

    /* Get size. */
    size_t
    size() const;
  private:
    double &
    unsafe_at(size_t row, size_t col);
    double
    unsafe_at(size_t row, size_t col) const;
    void
    print_row(std::ostream & stream, size_t row) const;
    void
    print_col(std::ostream & stream, size_t col) const;

    size_t rows_;
    size_t columns_;
    std::vector<double> data_;
  };

  std::ostream &
  operator<<(std::ostream &stream, const Matrix & A);

  bool
  operator==(const Matrix & lhs, const Matrix & rhs);
}