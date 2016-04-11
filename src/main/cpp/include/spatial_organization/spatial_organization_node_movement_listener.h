#ifndef SPATIAL_ORGANIZATION_SPATIAL_ORGANIZATION_NODE_MOVEMENT_LISTENER_H_
#define SPATIAL_ORGANIZATION_SPATIAL_ORGANIZATION_NODE_MOVEMENT_LISTENER_H_

#include <vector> //former list
#include <array>
#include <memory>
#include "physics/physical_node.h"
namespace cx3d {
namespace spatial_organization {

// forward declaration
 class SpaceNode;

//TODO change all occurences to SpatialOrganizationNode once porting has been finished

class SpatialOrganizationNodeMovementListener {
 public:
  SpatialOrganizationNodeMovementListener(){

  }
  virtual ~SpatialOrganizationNodeMovementListener() {
  }

  virtual void nodeAboutToMove(const SpaceNode* node,
                               const std::array<double, 3>& planned_movement) = 0;

  virtual void nodeMoved(const SpaceNode* node) = 0;

  virtual void nodeAboutToBeRemoved(const SpaceNode* node) = 0;

  virtual void nodeRemoved(const SpaceNode* node) = 0;

  virtual void nodeAboutToBeAdded(
      const SpaceNode* node, const std::array<double, 3>& planned_position,
      const std::array<physics::PhysicalNode*, 4>& vertices_of_the_tetrahedron_containing_the_position) = 0;

  virtual void nodeAdded(const SpaceNode* node) = 0;

  /**
   * Returns a String representation of this SpatialOrganizationNodeMovementListener
   */
  virtual std::string toString() const = 0;

  virtual bool equalTo(const SpatialOrganizationNodeMovementListener* other) const {
    return this == other;
  }
};

}  // namespace spatial_organization
}  // namespace cx3d

#endif  // SPATIAL_ORGANIZATION_SPATIAL_ORGANIZATION_NODE_MOVEMENT_LISTENER_H_
