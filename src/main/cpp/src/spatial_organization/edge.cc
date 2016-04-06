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


Edge::Edge(const std::shared_ptr<SpaceNode>& a, const std::shared_ptr<SpaceNode>& b)
    : a_(a),
      b_(b),
      cross_section_area_(0.0),
      adjacent_tetrahedra_() {
}


std::shared_ptr<SpaceNode> Edge::getOpposite(
    const std::shared_ptr<const SpaceNode>& node) const {
  if (node == a_) {
    return b_;
  } else if (node == b_) {
    return a_;
  } else {
    //fnoexceptionthrow std::invalid_argument(
//        "The edge " + toString() + " is not adjacent to the node " + node->toString());
  }
}


std::shared_ptr<physics::PhysicalNode> Edge::getOppositeElement(const std::shared_ptr<physics::PhysicalNode>& element) const {
  if (a_.get() != nullptr && b_.get() != nullptr) {
    if (element == a_->getUserObject()) {
      return b_->getUserObject();
    } else {
      return a_->getUserObject();
    }
  }
  return std::shared_ptr<physics::PhysicalNode>(nullptr);
}


std::shared_ptr<physics::PhysicalNode> Edge::getFirstElement() const {
  return a_->getUserObject();
}


std::shared_ptr<physics::PhysicalNode> Edge::getSecondElement() const {
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


bool Edge::equalTo(const std::shared_ptr<Edge>& other) {
  return other.get() == this;
}


bool Edge::equals(const std::shared_ptr<SpaceNode>& a,
                     const std::shared_ptr<SpaceNode>& b) const {
  return ((a_ == a) && (b_ == b)) || ((b_ == a) && (a_ == b));
}


void Edge::removeTetrahedron(const std::shared_ptr<Tetrahedron>& tetrahedron) {
  STLUtil::vectorRemove(adjacent_tetrahedra_, tetrahedron);
  if (adjacent_tetrahedra_.empty()) {
    remove();
  }
}


void Edge::addTetrahedron(const std::shared_ptr<Tetrahedron>& tetrahedron) {
  adjacent_tetrahedra_.push_back(tetrahedron);
}


void Edge::remove() {
  if (a_.get() != nullptr) {
    a_->removeEdge(this->shared_from_this());
  }
  if (b_.get() != nullptr) {
    b_->removeEdge(this->shared_from_this());
  }
}


std::vector<std::shared_ptr<Tetrahedron> > Edge::getAdjacentTetrahedra() const {
  return adjacent_tetrahedra_;
}


void Edge::changeCrossSectionArea(double change) {
  cross_section_area_ += change;
}

void Edge::initializationHelper() {
  if (a_.get() != nullptr) {
    a_->addEdge(this->shared_from_this());
  }
  if (b_.get() != nullptr) {
    b_->addEdge(this->shared_from_this());
  }
}

}  // namespace spatial_organization
}  // namespace cx3d
