#include "spatial_organization/flat_tetrahedron.h"

#include "physics/physical_node.h"
#include "spatial_organization/edge.h"
#include "spatial_organization/triangle_3d.h"
#include "spatial_organization/space_node.h"
#include "spatial_organization/open_triangle_organizer.h"

namespace cx3d {
namespace spatial_organization {



void FlatTetrahedron::updateCirumSphereAfterNodeMovement(
    const std::shared_ptr<SpaceNode>& moved_node) {
  for (size_t i = 0; i < 4; i++) {
    if (this->adjacent_nodes_[i] != moved_node) {
      this->adjacent_triangles_[i]->informAboutNodeMovement();
    }
  }
}


void FlatTetrahedron::calculateVolume() {
  this->volume_ = 0.0;
}


void FlatTetrahedron::updateCrossSectionAreas() {
  for (size_t i = 0; i < 6; i++) {
    this->changeCrossSection(i, 0.0);
  }
}


void FlatTetrahedron::calculateCircumSphere() {
}


bool FlatTetrahedron::isFlat() const {
  return true;
}


int FlatTetrahedron::orientation(const std::array<double, 3>& point) {
  this->adjacent_triangles_[0]->updatePlaneEquationIfNecessary();
  int orientation = this->adjacent_triangles_[0]->orientation(point, point);
  if (orientation == 0) {
    int memory = -1;
    for (size_t i = 0; i < 4; i++) {
      if (this->adjacent_triangles_[i] != nullptr) {
        int dummy = this->adjacent_triangles_[i]->circleOrientation(point);
        if (dummy == 1) {
          return 1;
        } else if (dummy == 0) {
          memory = 0;
        }
      }
    }
    return memory;
  } else {
    return orientation;
  }
}


bool FlatTetrahedron::isTrulyInsideSphere(const std::array<double, 3>& point) {
  return orientation(point) > 0;
}


bool FlatTetrahedron::isInsideSphere(const std::array<double, 3>& point) {
  return orientation(point) >= 0;
}


bool FlatTetrahedron::isPointInConvexPosition(const std::array<double, 3>& point,
                                                 int connecting_triangle_number) const {
  this->adjacent_triangles_[0]->updatePlaneEquationIfNecessary();
  return this->adjacent_triangles_[0]->orientation(point, point) == 0;
}


int FlatTetrahedron::isInConvexPosition(const std::array<double, 3>& point,
                                           int connecting_triangle_number) const {
  this->adjacent_triangles_[0]->updatePlaneEquationIfNecessary();
  if (this->adjacent_triangles_[0]->orientation(point, point) == 0) {
    return 0;
  } else {
    return -1;
  }
}


FlatTetrahedron::FlatTetrahedron()
    : Tetrahedron() {
}


}  // namespace spatial_organization
}  // namespace cx3d
