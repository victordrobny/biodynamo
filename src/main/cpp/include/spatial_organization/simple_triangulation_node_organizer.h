#ifndef SPATIAL_ORGANIZATION_SIMPLE_TRIANGULATION_NODE_ORGANIZER_H_
#define SPATIAL_ORGANIZATION_SIMPLE_TRIANGULATION_NODE_ORGANIZER_H_

#include <spatial_organization/abstract_triangulation_node_organizer.h>

namespace cx3d {
namespace spatial_organization {

 class BinaryTreeElement;

/**
 * This class is a very simple implementation of {@link AbstractTriangulationNodeOrganizer}.
 * All nodes are stored in a binary tree in order to obtain a good performance when
 * checking whether a certain node is already added to this organizer.
 *
 * @param  The type of user objects associated with the nodes in the current triangulation.
 */

class SimpleTriangulationNodeOrganizer :
    public AbstractTriangulationNodeOrganizer {
 public:
  static std::shared_ptr<SimpleTriangulationNodeOrganizer> create() {
    return std::shared_ptr<SimpleTriangulationNodeOrganizer>(
        new SimpleTriangulationNodeOrganizer());
  }

  SimpleTriangulationNodeOrganizer();

  virtual ~SimpleTriangulationNodeOrganizer();

  virtual void removeNode(const std::shared_ptr<SpaceNode>& node) override;

  virtual void addNode(const std::shared_ptr<SpaceNode>& node) override;

  virtual std::shared_ptr<SpaceNode> getFirstNode() const override;

  virtual std::string toString() const override;

  // TODO should be implemented in AbstractTriangulationNodeOrganizer, but SWIG does not generate
  // this function on the java side
  virtual void addTriangleNodes(const std::shared_ptr<Triangle3D>& triangle) override;

  std::vector<std::shared_ptr<SpaceNode>>getNodes(const std::shared_ptr<SpaceNode>& reference_point) override;

  bool equalTo(const std::shared_ptr<SimpleTriangulationNodeOrganizer>& other);

private:
  BinaryTreeElement* tree_head_;

  SimpleTriangulationNodeOrganizer(const SimpleTriangulationNodeOrganizer&) = delete;
  SimpleTriangulationNodeOrganizer& operator=(const SimpleTriangulationNodeOrganizer&) = delete;
};

}  // namespace spatial_organization
}  // namespace cx3d

#endif // SPATIAL_ORGANIZATION_SIMPLE_TRIANGULATION_NODE_ORGANIZER_H_
