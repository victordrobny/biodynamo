#include "spatial_organization/edge_hash_key.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "physics/physical_node.h"
#include "matrix.h"
#include "spatial_organization/space_node.h"

namespace cx3d {
namespace spatial_organization {


EdgeHashKey::EdgeHashKey(const std::shared_ptr<SpaceNode>& a,
                            const std::shared_ptr<SpaceNode>& b,
                            const std::shared_ptr<SpaceNode>& opposite_node)
    : a_(a),
      b_(b),
      hash_code_(0),
      ab_ { 0.0, 0.0, 0.0 },
      last_normal_vector_ { 0.0, 0.0, 0.0 } {
  ab_ = Matrix::subtract(b_->getPosition(), a_->getPosition());
  auto subtraction = Matrix::subtract(opposite_node->getPosition(), a_->getPosition());
  last_normal_vector_ = Matrix::normalize(Matrix::crossProduct(ab_, subtraction));
  hash_code_ = std::max(a_->getId(), b_->getId()) * 11 + std::min(a_->getId(), b_->getId()) * 31;
}


EdgeHashKey::EdgeHashKey(const EdgeHashKey& other){
  hash_code_ = other.hash_code_;
}


EdgeHashKey& EdgeHashKey::operator=(const EdgeHashKey& rhs) {
  if(this == &rhs){
    return *this;
  }
  a_ = rhs.a_;
  b_ = rhs.b_;
  hash_code_ = rhs.hash_code_;
  ab_ = rhs.ab_;
  last_normal_vector_ = rhs.last_normal_vector_;
  return *this;
}


std::string EdgeHashKey::toString() const {
  std::ostringstream str_stream;
  str_stream << "(";
  str_stream << a_->toString();
  str_stream << ", ";
  str_stream << b_->toString();
  str_stream << ")";
  return str_stream.str();
}


int EdgeHashKey::hashCode() const {
  return hash_code_;
}


bool EdgeHashKey::equalTo(const std::shared_ptr<EdgeHashKey>& other) const {
  return other.get() == this;
}


double EdgeHashKey::getCosine(const std::array<double, 3>& fourth_point) const {
  auto difference = Matrix::subtract(fourth_point, a_->getPosition());
  auto cross_product = Matrix::crossProduct(ab_, difference);
  auto normal = Matrix::normalize(cross_product);
  double cosine = Matrix::dot(normal, last_normal_vector_);
  if (cosine > 0.999999999) {
    return 1;
  } else if (cosine < -0.99999999) {
    return -1;
  }
  return cosine;
}


std::shared_ptr<SpaceNode> EdgeHashKey::getEndpointA() const {
  return a_;
}


std::shared_ptr<SpaceNode> EdgeHashKey::getEndpointB() const {
  return b_;
}


std::shared_ptr<SpaceNode> EdgeHashKey::oppositeNode(
    const std::shared_ptr<SpaceNode>& node) const {
  if (node == a_) {
    return b_;
  } else if (node == b_) {
    return a_;
  } else {
    //fnoexceptionthrow std::invalid_argument("Could not find an opposite node for" + node->toString());
  }
}


std::size_t EdgeHashKeyHash::operator()( const EdgeHashKey& element) const {
  return element.hash_code_;
}


bool EdgeHashKeyEqual::operator()(const EdgeHashKey& lhs, const EdgeHashKey& rhs) const {
  return lhs.hash_code_ == rhs.hash_code_;
}


}  // namespace spatial_organization
}  // namespace cx3d
