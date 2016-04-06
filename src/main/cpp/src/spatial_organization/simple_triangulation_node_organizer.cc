#include "spatial_organization/simple_triangulation_node_organizer.h"

#include <sstream>

#include "string_util.h"
#include "physics/physical_node.h"
#include "spatial_organization/space_node.h"
#include "spatial_organization/triangle_3d.h"
#include "spatial_organization/binary_tree_element.h"

namespace cx3d {
namespace spatial_organization {


SimpleTriangulationNodeOrganizer::SimpleTriangulationNodeOrganizer()
    : tree_head_ { BinaryTreeElement::generateTreeHead() } {
}


SimpleTriangulationNodeOrganizer::~SimpleTriangulationNodeOrganizer() {
  delete tree_head_;
}


std::vector<std::shared_ptr<SpaceNode>>SimpleTriangulationNodeOrganizer::getNodes(const std::shared_ptr<SpaceNode>& reference_point) {
  return tree_head_->inOrderTraversal();
}


void SimpleTriangulationNodeOrganizer::removeNode(
    const std::shared_ptr<SpaceNode>& node) {
  tree_head_->remove(node, nullptr);
}


void SimpleTriangulationNodeOrganizer::addNode(
    const std::shared_ptr<SpaceNode>& node) {
  tree_head_->insert(node);
}


std::shared_ptr<SpaceNode> SimpleTriangulationNodeOrganizer::getFirstNode() const {
  return tree_head_->bigger_->content_;
}


std::string SimpleTriangulationNodeOrganizer::toString() const {
  std::stringstream str_stream;
  str_stream << "[";
  str_stream << StringUtil::toStr(tree_head_);
  str_stream << "]";
  return str_stream.str();
}

// TODO move to AbstractTriangulationNodeOrganizer once porting has been finished
void SimpleTriangulationNodeOrganizer::addTriangleNodes(
    const std::shared_ptr<Triangle3D>& triangle) {
  auto nodes = triangle->getNodes();
  addNode(nodes[1]);
  addNode(nodes[2]);
  addNode(nodes[0]);
}

bool SimpleTriangulationNodeOrganizer::equalTo(
    const std::shared_ptr<SimpleTriangulationNodeOrganizer>& other) {
  return this == other.get();
}


}  // namespace spatial_organization
}  // namespace cx3d
