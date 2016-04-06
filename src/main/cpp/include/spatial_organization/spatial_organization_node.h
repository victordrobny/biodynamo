#ifndef SPATIAL_ORGANIZATION_SPATIAL_ORGANIZATION_NODE_H_
#define SPATIAL_ORGANIZATION_SPATIAL_ORGANIZATION_NODE_H_

#include <string>
#include <vector> //former list
#include <array>
#include <memory>

#include "string_builder.h"
#include "physics/physical_node.h"

namespace cx3d {
namespace spatial_organization {

 class Edge;
 class SpaceNode;
 class SpatialOrganizationEdge;
 class SpatialOrganizationNodeMovementListener;

/**
 * Interface to define the basic properties of a node in the triangulation.
 *
 * @param  The type of user objects associated with each node in the triangulation.
 */

class SpatialOrganizationNode {
 public:
  virtual ~SpatialOrganizationNode() {
  }

  virtual void addSpatialOrganizationNodeMovementListener(
      const std::shared_ptr<SpatialOrganizationNodeMovementListener>& listener) = 0;

  /**
   * Returns a list that allows to iterate over all edges
   * incident to this node.
   */
  virtual std::vector<std::shared_ptr<Edge> > getEdges() const = 0;  //TODO change to SpatialOrganizationEdge once porting has been finished

  virtual std::vector<std::shared_ptr<cx3d::physics::PhysicalNode>> getNeighbors() const = 0;

  // todo change to interface type
  virtual std::shared_ptr<SpaceNode> getNewInstance(
      const std::array<double, 3>& position,
      const std::shared_ptr<physics::PhysicalNode>& user_object) = 0;

  virtual std::vector<std::shared_ptr<cx3d::physics::PhysicalNode>> getPermanentListOfNeighbors() const = 0;

  virtual std::array<double, 3> getPosition() const = 0;

  virtual std::shared_ptr<physics::PhysicalNode> getUserObject() const = 0;

  virtual std::array<std::shared_ptr<cx3d::physics::PhysicalNode>, 4> getVerticesOfTheTetrahedronContaining(
      const std::array<double, 3>& position, std::array<int, 1>& returned_null) const = 0;

  virtual double getVolume() const = 0;

  virtual void moveFrom(const std::array<double, 3>& delta) = 0;

  virtual void remove() = 0;

  virtual std::string toString() const =0;
};

}  // namespace spatial_organization
}  // namespace cx3d

#endif  // SPATIAL_ORGANIZATION_SPATIAL_ORGANIZATION_NODE_H_
