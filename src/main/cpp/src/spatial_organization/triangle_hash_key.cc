#include "spatial_organization/triangle_hash_key.h"

#include <cmath>

#include "physics/physical_node.h"
#include "spatial_organization/space_node.h"

namespace cx3d {
namespace spatial_organization {


TriangleHashKey::TriangleHashKey(SpaceNode* a,
                                    SpaceNode* b,
                                    SpaceNode* c)
    : a_(a),
      b_(b),
      c_(c),
      hash_code_(0) {
  int a_id = a_ != nullptr ? a_->getId() : -1;
  int b_id = b_ != nullptr ? b_->getId() : -1;
  int c_id = c_ != nullptr ? c_->getId() : -1;
  createHashCode(a_id, b_id, c_id);
}


TriangleHashKey::TriangleHashKey(const TriangleHashKey& other) {
    hash_code_ = other.hash_code_;
    a_ = other.a_;
    b_ = other.b_;
    c_ = other.c_;
  }


int TriangleHashKey::hashCode() const {
  return hash_code_;
}


bool TriangleHashKey::equalTo(const TriangleHashKey* other) const {
  return (a_ == other->a_
      && ((b_ == other->b_ && c_ == other->c_) || (b_ == other->c_ && c_ == other->b_)))
      || (a_ == other->b_
          && ((b_ == other->a_ && c_ == other->c_) || (b_ == other->c_ && c_ == other->a_)))
      || (a_ == other->c_
          && ((b_ == other->a_ && c_ == other->b_) || (b_ == other->b_ && c_ == other->a_)));
}


void TriangleHashKey::createHashCode(int a_id, int b_id, int c_id) {
  int min = std::min(a_id, std::min(b_id, c_id));
  int max = std::max(a_id, std::max(b_id, c_id));
  hash_code_ = (min * 31 + max * 11 + a_id + b_id + c_id) % 2000000001;
}


std::size_t TriangleHashKeyHash::operator()(const TriangleHashKey& key) const {
  return key.hash_code_;
}


bool TriangleHashKeyEqual::operator()(const TriangleHashKey& lhs,
                                         const TriangleHashKey& rhs) const {
  return (lhs.a_ == rhs.a_
      && ((lhs.b_ == rhs.b_ && lhs.c_ == rhs.c_)
          || (lhs.b_ == rhs.c_ && lhs.c_ == rhs.b_)))
      || (lhs.a_ == rhs.b_
          && ((lhs.b_ == rhs.a_ && lhs.c_ == rhs.c_)
              || (lhs.b_ == rhs.c_ && lhs.c_ == rhs.a_)))
      || (lhs.a_ == rhs.c_
          && ((lhs.b_ == rhs.a_ && lhs.c_ == rhs.b_)
              || (lhs.b_ == rhs.b_ && lhs.c_ == rhs.a_)));
}


}  // namespace spatial_organization
}  // namespace cx3d

