#include "spatial_organization/exact_vector.h"

#include "spatial_organization/rational.h"

namespace cx3d {
namespace spatial_organization {

using std::array;
using std::shared_ptr;
using cx3d::spatial_organization::ExactVector;

Rational* ExactVector::det(const std::array<ExactVector*, 3>& c) {
  return c[0]->elements_[0]->multiply(c[1]->elements_[1])->multiply(c[2]->elements_[2])->add(
      c[0]->elements_[1]->multiply(c[1]->elements_[2])->multiply(c[2]->elements_[0]))->add(
      c[0]->elements_[2]->multiply(c[1]->elements_[0])->multiply(c[2]->elements_[1]))->subtract(
      c[0]->elements_[0]->multiply(c[1]->elements_[2])->multiply(c[2]->elements_[1]))->subtract(
      c[0]->elements_[1]->multiply(c[1]->elements_[0])->multiply(c[2]->elements_[2]))->subtract(
      c[0]->elements_[2]->multiply(c[1]->elements_[1])->multiply(c[2]->elements_[0]));
}

ExactVector::ExactVector(const array<Rational*, 3>& values)
    : elements_(values) {
}

ExactVector::ExactVector(const array<double, 3>& values)
    : elements_(
        array<Rational*, 3>({ Rational::create(values[0]), Rational::create(values[1]),
            Rational::create(values[2]) })) {
}

#ifndef EXACTVECTOR_NATIVE
ExactVector::ExactVector() : elements_(
      array<Rational*, 3>({ Rational::create(0.0), Rational::create(0.0),
          Rational::create(0.0) })) {
}
#endif

ExactVector::~ExactVector() {
}

Rational* ExactVector::squaredLength() const {
  auto rational = Rational::create(0L, 1L);
  for (auto element : elements_) {
    rational->add(element->multiply(element));
  }
  return rational;
}

ExactVector* ExactVector::add(const ExactVector* other) const {
  array<Rational*, 3> vector;
  vector[0] = elements_[0]->add(other->elements_[0]);
  vector[1] = elements_[1]->add(other->elements_[1]);
  vector[2] = elements_[2]->add(other->elements_[2]);

  return ExactVector::create(vector);
}

ExactVector* ExactVector::increaseBy(const ExactVector* other) {
  for (int i = 0; i < 3; i++) {
    elements_[i]->increaseBy(other->elements_[i]);
  }
  return this;
}

ExactVector* ExactVector::subtract(const ExactVector* other) {
  array<Rational*, 3> vector;
  vector[0] = elements_[0]->subtract(other->elements_[0]);
  vector[1] = elements_[1]->subtract(other->elements_[1]);
  vector[2] = elements_[2]->subtract(other->elements_[2]);

  return ExactVector::create(vector);
}

ExactVector* ExactVector::decreaseBy(const ExactVector* other) {
  for (int i = 0; i < 3; i++) {
    elements_[i]->decreaseBy(other->elements_[i]);
  }
  return this;
}

ExactVector* ExactVector::multiply(const Rational* factor) {
  array<Rational*, 3> vector;
  vector[0] = elements_[0]->multiply(factor);
  vector[1] = elements_[1]->multiply(factor);
  vector[2] = elements_[2]->multiply(factor);

  return ExactVector::create(vector);
}

ExactVector* ExactVector::multiplyBy(const Rational* factor) {
  for (auto element : elements_) {
    element->multiplyBy(factor);
  }
  return this;
}

ExactVector* ExactVector::divide(const Rational* factor) {
  array<Rational*, 3> vector;
  vector[0] = elements_[0]->divide(factor);
  vector[1] = elements_[1]->divide(factor);
  vector[2] = elements_[2]->divide(factor);

  return ExactVector::create(vector);
}

ExactVector* ExactVector::divideBy(const Rational* factor) {
  for (auto element : elements_) {
    element->divideBy(factor);
  }
  return this;
}

Rational* ExactVector::dotProduct(const ExactVector* other) {
  auto rational = Rational::create(0L, 1L);
  for (int i = 0; i < 3; i++) {
    rational = rational->add(other->elements_[i]->multiply(elements_[i]));
  }
  return rational;
}

ExactVector* ExactVector::negate() {
  for (auto element : elements_) {
    element->negate();
  }
  return this;
}

ExactVector* ExactVector::crossProduct(const ExactVector* other) {
  array<Rational*, 3> vector;
  for (int i = 0; i < 3; i++) {
    vector[i] = elements_[((i + 1) % 3)]->multiply(other->elements_[((i + 2) % 3)])->subtract(
        elements_[((i + 2) % 3)]->multiply(other->elements_[((i + 1) % 3)]));
  }
  return ExactVector::create(vector);
}

std::string ExactVector::toString() {
  std::stringstream ret;
  ret << "(";
  for (size_t i = 0; i < elements_.size(); i++) {
    if (i != 0)
      ret << ", ";
    ret << elements_[i]->toString();
  }
  ret << ")";
  return ret.str();
}

bool ExactVector::equalTo(const ExactVector* other) {
  for (size_t i = 0; i < elements_.size(); i++) {
    if (elements_[i]->compareTo(other->elements_[i]) != 0) {
      return false;
    }
  }
  return true;
}

}  // namespace spatial_organization
}  // namespace cx3d
