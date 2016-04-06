#ifndef SPATIAL_ORGANIZATION_ABSTRACT_TRIANGULATION_NODE_ORGANIZER_H_
#define SPATIAL_ORGANIZATION_ABSTRACT_TRIANGULATION_NODE_ORGANIZER_H_

#include <vector> //former list
#include <memory>
#include <string>

namespace cx3d {
namespace spatial_organization {

 class SpaceNode;
 class Triangle3D;

/**
 * Instances of child classes of this class are used to keep track of
 * nodes in an incomplete triangulation which might possibly become neighbors of open triangles.
 *
 * @param  The type of user objects associated with nodes in the current triangulation.
 */
class AbstractTriangulationNodeOrganizer {
 public:
  AbstractTriangulationNodeOrganizer() {
  }

  virtual ~AbstractTriangulationNodeOrganizer() {
  }

  virtual std::vector<std::shared_ptr<SpaceNode> > getNodes(const std::shared_ptr<SpaceNode>& reference_point) = 0;

  virtual void addTriangleNodes(const std::shared_ptr<Triangle3D>& triangle) = 0;

  virtual void removeNode(const std::shared_ptr<SpaceNode>& node) = 0;

  virtual void addNode(const std::shared_ptr<SpaceNode>& node) = 0;

  virtual std::shared_ptr<SpaceNode> getFirstNode() const = 0;

  virtual std::string toString() const = 0;

 private:
  AbstractTriangulationNodeOrganizer(const AbstractTriangulationNodeOrganizer&) = delete;
  AbstractTriangulationNodeOrganizer& operator=(const AbstractTriangulationNodeOrganizer& ) = delete;
};

}  // namespace spatial_organization
}  // namespace cx3d

#endif // SPATIAL_ORGANIZATION_ABSTRACT_TRIANGULATION_NODE_ORGANIZER_H_
