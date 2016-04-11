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
  static SimpleTriangulationNodeOrganizer* create() {
    return new SimpleTriangulationNodeOrganizer();
  }

  SimpleTriangulationNodeOrganizer();

  virtual ~SimpleTriangulationNodeOrganizer();

  virtual void removeNode(const SpaceNode* node) override;

  virtual void addNode(SpaceNode* node) override;

  virtual SpaceNode* getFirstNode() const override;

  virtual std::string toString() const override;

  // TODO should be implemented in AbstractTriangulationNodeOrganizer, but SWIG does not generate
  // this function on the java side
  virtual void addTriangleNodes(const Triangle3D* triangle) override;

  std::vector<SpaceNode*>getNodes(const SpaceNode* reference_point) override;

  bool equalTo(const SimpleTriangulationNodeOrganizer* other);

private:
  BinaryTreeElement* tree_head_;

  SimpleTriangulationNodeOrganizer(const SimpleTriangulationNodeOrganizer&) = delete;
  SimpleTriangulationNodeOrganizer& operator=(const SimpleTriangulationNodeOrganizer&) = delete;
};

}  // namespace spatial_organization
}  // namespace cx3d

#endif // SPATIAL_ORGANIZATION_SIMPLE_TRIANGULATION_NODE_ORGANIZER_H_
