#include "spatial_organization/edge.h"

#include <sstream>
#include <stdexcept>
#include "stl_util.h"

#include "physics/physical_node.h"
#include "string_util.h"
#include "spatial_organization/space_node.h"
#include "spatial_organization/tetrahedron.h"

namespace cx3d {
namespace spatial_organization {


Edge::Edge(SpaceNode* a, SpaceNode* b)
    : a_(a),
      b_(b),
      cross_section_area_(0.0),
      adjacent_tetrahedra_() {
}


SpaceNode* Edge::getOpposite(
    const SpaceNode* node) const {
  if (node == a_) {
    return b_;
  } else if (node == b_) {
    return a_;
  } else {
    //fnoexceptionthrow std::invalid_argument(
//        "The edge " + toString() + " is not adjacent to the node " + node->toString());
  }
}


physics::PhysicalNode* Edge::getOppositeElement(const physics::PhysicalNode* element) const {
  if (a_ != nullptr && b_ != nullptr) {
    if (element == a_->getUserObject()) {
      return b_->getUserObject();
    } else {
      return a_->getUserObject();
    }
  }
  return nullptr;
}


physics::PhysicalNode* Edge::getFirstElement() const {
  return a_->getUserObject();
}


physics::PhysicalNode* Edge::getSecondElement() const {
  return b_->getUserObject();
}


double Edge::getCrossSection() const {
  return cross_section_area_;
}


const std::string Edge::toString() const {
  std::ostringstream str_stream;
  str_stream << "(";
  str_stream << "Edge";
  str_stream << StringUtil::toStr(cross_section_area_);
//  str_stream << " - ";
//  str_stream << StringUtil::toStr(a_);
//  str_stream << " - ";
//  str_stream << StringUtil::toStr(b_);
//  str_stream << " - ";
//  str_stream << StringUtil::toStr(adjacent_tetrahedra_);
  str_stream << ")";
  return str_stream.str();
}


bool Edge::equalTo(const Edge* other) {
  return other == this;
}


bool Edge::equals(const SpaceNode* a,
                     const SpaceNode* b) const {
  return ((a_ == a) && (b_ == b)) || ((b_ == a) && (a_ == b));
}


void Edge::removeTetrahedron(Tetrahedron* tetrahedron) {
  STLUtil::vectorRemove(adjacent_tetrahedra_, tetrahedron);
  if (adjacent_tetrahedra_.empty()) {
    remove();
  }
}


void Edge::addTetrahedron(Tetrahedron* tetrahedron) {
  adjacent_tetrahedra_.push_back(tetrahedron);
}


void Edge::remove() {
  if (a_ != nullptr) {
    a_->removeEdge(this);
  }
  if (b_ != nullptr) {
    b_->removeEdge(this);
  }
}


std::vector<Tetrahedron* > Edge::getAdjacentTetrahedra() const {
  return adjacent_tetrahedra_;
}


void Edge::changeCrossSectionArea(double change) {
  cross_section_area_ += change;
}

void Edge::initializationHelper() {
  if (a_ != nullptr) {
    a_->addEdge(this);
  }
  if (b_ != nullptr) {
    b_->addEdge(this);
  }
}

}  // namespace spatial_organization
}  // namespace cx3d
