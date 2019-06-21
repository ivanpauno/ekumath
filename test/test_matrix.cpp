#include <cmath>
#include <sstream>

#include "gtest/gtest.h"

#include "ekumath.hpp"

namespace ekumath {
namespace {

TEST(MatrixTest, MatrixInitialization)
{
  Matrix Z(2, 2);
  EXPECT_EQ(Z.rows(), 2);
  EXPECT_EQ(Z.cols(), 2);
  EXPECT_EQ(Z.size(), 4);
  EXPECT_EQ(Z, Matrix({{0., 0.}, {0., 0.}}));
  Matrix A(3, 2, -1.);
  EXPECT_EQ(A.rows(), 3);
  EXPECT_EQ(A.cols(), 2);
  EXPECT_EQ(A.size(), 6);
  EXPECT_EQ(A, Matrix(
    {{-1., -1.}, {-1., -1.}, {-1., -1.}}));
  Matrix B{0., 1., 2.};
  EXPECT_EQ(B.rows(), 3);
  EXPECT_EQ(B.cols(), 1);
  EXPECT_EQ(B.size(), 3);
  // Note(ivanpauno): Preferred size_t over int. Avoided comparsion between signed and unsigned warning.
  // for (int i = 0; i < B.rows(); ++i) {
  for (size_t i = 0; i < B.rows(); ++i) {
    EXPECT_EQ(B(i), B(i, 0));
    EXPECT_DOUBLE_EQ(B(i), static_cast<double>(i));
  }
  // NOTE(ivanpauno): row vectors aren't equal to column vectors
  // EXPECT_EQ(B, Matrix({{0., 1., 2.}});
  // NOTE(ivanpauno): Modified test logic error
  // Matrix C{{0., 1.}, {2., 3.}};
  Matrix C{{0., 1.}, {1., 2.}};
  EXPECT_EQ(C.rows(), 2);
  EXPECT_EQ(C.cols(), 2);
  EXPECT_EQ(C.size(), 4);
  // Note(ivanpauno): Preferred size_t over int. Avoided comparsion between signed and unsigned warning.
  // for (int i = 0; i < C.rows(); ++i) {
    // for (int j = 0; j < C.cols(); ++j) {
  for (size_t i = 0; i < C.rows(); ++i) {
    for (size_t j = 0; j < C.cols(); ++j) {
      EXPECT_DOUBLE_EQ(C(i, j), static_cast<double>(i + j));
    }
  }
  EXPECT_THROW(Matrix D({{1., 2.}, {3., 1., 3.}}), std::runtime_error);
  Matrix I = Matrix::Identity(3);
  EXPECT_EQ(I.rows(), 3);
  EXPECT_EQ(I.rows(), I.cols());
  EXPECT_EQ(I, Matrix({{1., 0., 0.},
                      {0., 1., 0.},
                      {0., 0., 1.}}));
  Matrix R = Matrix::Random(3, 2);
  EXPECT_EQ(R.rows(), 3);
  EXPECT_EQ(R.cols(), 2);
  EXPECT_NE(R, Matrix::Random(3, 2));
}

TEST(MatrixTest, ElementWiseArithmetic) {
  Matrix A(3, 2, 1.);
  Matrix B = (25 * (A + 2) - 15.) / 2.;
  // Note(ivanpauno): Preferred size_t over int. Avoided comparsion between signed and unsigned warning.
  // for (int i = 0; i < A.rows(); ++i) {
    // for (int j = 0; j < A.cols(); ++j) {
  for (size_t i = 0; i < A.rows(); ++i) {
    for (size_t j = 0; j < A.cols(); ++j) {
      EXPECT_DOUBLE_EQ(B(i, j), (25 * (A(i, j) + 2) - 15.) / 2.);
    }
  }
  Matrix C = A * B;
  // Note(ivanpauno): Preferred size_t over int. Avoided comparsion between signed and unsigned warning.
  // for (int i = 0; i < A.rows(); ++i) {
    // for (int j = 0; j < A.cols(); ++j) {
  for (size_t i = 0; i < A.rows(); ++i) {
    for (size_t j = 0; j < A.cols(); ++j) {
      EXPECT_DOUBLE_EQ(C(i, j), A(i, j) * B(i, j));
    }
  }
  Matrix D = C + B;
  // Note(ivanpauno): Preferred size_t over int. Avoided comparsion between signed and unsigned warning.
  // for (int i = 0; i < D.rows(); ++i) {
    // for (int j = 0; j < D.cols(); ++j) {
  for (size_t i = 0; i < D.rows(); ++i) {
    for (size_t j = 0; j < D.cols(); ++j) {
      EXPECT_DOUBLE_EQ(D(i, j), C(i, j) + B(i, j));
    }
  }

  using std::pow;
  using std::sqrt;
  Matrix R{{25., 16.}, {9., 100.}};
  EXPECT_EQ(R ^ 3, pow(R, 3));
  // NOTE(ivanpauno): Corrected math error.
  // EXPECT_EQ(pow(R, 3), Matrix({{15625., 64.}, {27., 1000000.}}));
  EXPECT_EQ(pow(R, 3), Matrix({{15625., 4096.}, {729., 1000000.}}));
  EXPECT_EQ(sqrt(R), Matrix({{5., 4.}, {3., 10.}}));
}

TEST(MatrixTest, BlockManipulation) {
  Matrix I = Matrix::Identity(3);
  EXPECT_EQ(I, I.transpose());
  Matrix R = Matrix::Random(3);
  // NOTE(ivanpauno): Renamed transpose as get_transpose.
  // The transpose method modifies the internal state of R.
  // EXPECT_NE(R, R.transpose());
  EXPECT_NE(R, R.get_transpose());
  // NOTE(ivanpauno): S(0, span) returns something that has shared state with S. That's tricky.
  // Matrix S(3, 3);
  // S(0, span) = R(span, 0);
  // S(1, span) = R(span, 1);
  // S(2, span) = R(span, 2);
  // EXPECT_EQ(S, R.transpose());
  Matrix P = I | Matrix(3, 1);
  // NOTE(ivanpauno): Corrected wrong expected result
  // EXPECT_EQ(P, Matrix({{0., 0., 1., 0.},
                      // {1., 0., 0., 0.},
                      // {0., 1., 0., 0.}}));
  EXPECT_EQ(P, Matrix({{1., 0., 0., 0.},
                      {0., 1., 0., 0.},
                      {0., 0., 1., 0.}}));
  EXPECT_EQ(P, I.concatenate(Matrix(3, 1)));
  EXPECT_THROW(P | Matrix({1., 2.}), std::runtime_error);
  // NOTE(ivanpauno): This is not just tricky, but actually wrong.
  // Matrix Q = Matrix::Identity(3);
  // Q(span, {0, 1, 2}) = Q(span, {1, 2, 0});
  // EXPECT_EQ(Q, Matrix({{0., 0., 1.},
                      // {1., 0., 0.},
                      // {0., 1., 0.}}));
}

TEST(MatrixTest, IvanModifiedBlockManipulation) {
  Matrix P = Matrix({
    {0., 1., 2., 3.},
    {4., 5., 6., 7.},
    {8., 9., 10., 11.}});
  EXPECT_EQ(P({0, 1}, 1), Matrix({std::initializer_list<double>({1.}), std::initializer_list<double>({5.})}));
  EXPECT_EQ(P(1, {0, 1}), Matrix({{4., 5.}}));
  EXPECT_EQ(P({1, 2}, {1, 2}), Matrix({{5., 6.}, {9., 10.}}));
  EXPECT_EQ(P({1, 2}, span), Matrix({{4., 5., 6., 7.}, {8., 9., 10., 11.}}));
  EXPECT_EQ(P(span, {1, 2}), Matrix({{1., 2.}, {5., 6.}, {9., 10.}}));
}

TEST(MatrixTest, IOSupport) {
  std::stringstream ss0;
  // NOTE(ivanpauno): Corrected bad usage of setprecision.
  // ss0 << std::showpoint << std::setprecision(3) << Matrix::Identity(3) << std::endl;
  ss0 << std::showpoint << std::setprecision(4) << Matrix::Identity(3) << std::endl;
  EXPECT_EQ(ss0.str(), "[[1.000, 0.000, 0.000], [0.000, 1.000, 0.000], [0.000, 0.000, 1.000]]\n");
  std::stringstream ss1;
  // NOTE(ivanpauno): Corrected bad usage of setprecision.
  // ss1 << std::showpoint << std::setprecision(2) << Matrix{1., 2.5, 3.} << std::endl;
  ss1 << std::showpoint << std::setprecision(3) << Matrix{1., 2.5, 3.} << std::endl;
  EXPECT_EQ(ss1.str(), "[1.00, 2.50, 3.00]\n");
}

}  // namespace
}  // namespace ekumath

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
