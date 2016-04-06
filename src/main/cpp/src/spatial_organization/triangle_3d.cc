#include "spatial_organization/triangle_3d.h"

#include <stdint.h>
#include <cmath>
#include <array>
#include <limits>
#include <memory>

#include "matrix.h"
#include "string_util.h"
#include "physics/physical_node.h"
#include "spatial_organization/tetrahedron.h"
#include "spatial_organization/space_node.h"

namespace cx3d {
namespace spatial_organization {

using std::shared_ptr;


std::array<double, 3> Triangle3D::calculate3PlaneXPoint(
    const std::array<std::array<double, 3>, 3>& normals, const std::array<double, 3>& offsets,
    double normal_det) {
  std::array<double, 3> result;
  if (normal_det != 0.0) {
    auto vec_1 = Matrix::scalarMult(offsets[0], Matrix::crossProduct(normals[1], normals[2]));
    auto vec_2 = Matrix::scalarMult(offsets[1], Matrix::crossProduct(normals[2], normals[0]));
    auto vec_3 = Matrix::scalarMult(offsets[2], Matrix::crossProduct(normals[0], normals[1]));

    auto sum_vec = Matrix::add(vec_1, Matrix::add(vec_2, vec_3));
    result = Matrix::scalarMult(1 / normal_det, sum_vec);
  } else {
    double max_value = std::numeric_limits<double>::max();
    for (int i = 0; i < 3; i++) {
      result[i] = max_value;
    }
  }
  return result;
}


std::array<double, 3> Triangle3D::calculate3PlaneXPoint(
    const std::array<std::array<double, 3>, 3>& normals, const std::array<double, 3>& offsets) {
  return calculate3PlaneXPoint(normals, offsets, Matrix::det(normals));
}


std::shared_ptr<ExactVector> Triangle3D::calculate3PlaneXPoint(
    const std::array<std::shared_ptr<ExactVector>, 3>& normals,
    const std::array<std::shared_ptr<Rational>, 3>& offsets,
    const std::shared_ptr<Rational>& normal_det) {
  if (!normal_det->isZero()) {
    auto ret = normals[1]->crossProduct(normals[2])->multiplyBy(offsets[0])->increaseBy(
        normals[2]->crossProduct(normals[0])->multiplyBy(offsets[1])->increaseBy(
            normals[0]->crossProduct(normals[1])->multiplyBy(offsets[2])))->divideBy(normal_det);
#ifdef TRIANGLE3D_DEBUG
    logReturnStatic(ret);
#endif
    return ret;
  } else {
    double max_value = std::numeric_limits<double>::max();
    std::array<std::shared_ptr<Rational>, 3> rationals;
    for (int i = 0; i < 3; i++) {
      rationals[i] = Rational::create(max_value, 1);
    }
    auto ret = ExactVector::create(rationals);
#ifdef TRIANGLE3D_DEBUG
    logReturnStatic(ret);
#endif
    return ret;
  }
}


Triangle3D::Triangle3D(const std::shared_ptr<SpaceNode>& sn_1,
                          const std::shared_ptr<SpaceNode>& sn_2,
                          const std::shared_ptr<SpaceNode>& sn_3,
                          const std::shared_ptr<Tetrahedron>& tetrahedron_1,
                          const std::shared_ptr<Tetrahedron>& tetrahedron_2)
    : Plane3D(),
      adjacent_tetrahedra_({ tetrahedron_1, tetrahedron_2 }),
      nodes_({ sn_1, sn_2, sn_3 }),
      circum_center_ { 0.0, 0.0, 0.0 },
      plane_updated_(false),
      circum_center_updated_(false),
      upper_side_positive_(true),
      connection_checked_(-1) {

  if (sn_2.get() == nullptr) {
    nodes_[1] = sn_1;
    nodes_[0] = std::shared_ptr<SpaceNode>(nullptr);
  }
  if (sn_3.get() == nullptr) {
    nodes_[2] = sn_1;
    nodes_[0] = std::shared_ptr<SpaceNode>(nullptr);
  }
}


bool Triangle3D::isSimilarTo(const std::shared_ptr<Triangle3D>& other_triangle) const {
  auto other_nodes = other_triangle->getNodes();
  return isAdjacentTo(other_nodes[0]) && isAdjacentTo(other_nodes[1])
      && isAdjacentTo(other_nodes[2]);
}


double Triangle3D::getSDDistance(const std::array<double, 3>& fourth_point) const {
  if (!isInfinite() && onUpperSide(fourth_point)) {
    double sd_distance = calculateSDDistance(fourth_point);
    if (sd_distance != std::numeric_limits<double>::max()) {
      return (upper_side_positive_) ? sd_distance : -sd_distance;
    } else {
      return std::numeric_limits<double>::max();
    }
  } else {
    return std::numeric_limits<double>::max();
  }
}


std::shared_ptr<Rational> Triangle3D::getSDDistanceExact(
    const std::array<double, 3>& fourth_point) const {
  if (!isInfinite() && onUpperSide(fourth_point)) {
    std::array<std::shared_ptr<ExactVector>, 4> points;
    std::array<std::shared_ptr<ExactVector>, 3> points_3;
    for (int i = 0; i < 3; i++) {
      points[i] = ExactVector::create(nodes_[i]->getPosition());
      points_3[i] = points[i];
    }
    points[3] = ExactVector::create(fourth_point);
    auto normal_vector = calculateExactNormalVector(points_3);
    if (normal_vector->dotProduct(ExactVector::create(this->normal_vector_))->compareTo(
        Rational::create(0, 1)) < 0) {
      normal_vector->negate();
    }
    if (upper_side_positive_) {
      return calculateSDDistanceExact(points, normal_vector);
    } else {
      return calculateSDDistanceExact(points, normal_vector)->negate();
    }
  } else {
    return Rational::create(std::numeric_limits<int64_t>::max(), 1);
  }
}


std::array<double, 3> Triangle3D::calculateCircumSphereCenter(
    const std::array<double, 3>& fourth_point) const {
  if (!isInfinite()) {
    double sd = calculateSDDistance(fourth_point);
    return Matrix::add(circum_center_, Matrix::scalarMult(sd, this->normal_vector_));
  }
  //fnoexceptionthrow std::logic_error("could not calculate circum sphere, because triangle is infinite");
}


std::array<double, 3> Triangle3D::calculateCircumSphereCenterIfEasy(
    const std::array<double, 3>& fourth_point) const {
  if (circum_center_updated_) {
    return calculateCircumSphereCenter(fourth_point);
  }
  //fnoexceptionthrow std::logic_error("could not calculate circum sphere, because triangle is infinite");
}


void Triangle3D::informAboutNodeMovement() {
  circum_center_updated_ = false;
  plane_updated_ = false;
  this->normal_vector_updated_ = false;
}


void Triangle3D::updatePlaneEquationIfNecessary() {
  if (!plane_updated_ && !isInfinite()) {
    auto node_0_position = nodes_[0]->getPosition();
    auto diff_1 = Matrix::subtract(nodes_[1]->getPosition(), node_0_position);
    auto diff_2 = Matrix::subtract(nodes_[2]->getPosition(), node_0_position);
    this->initPlane(diff_1, diff_2, node_0_position, false);
    plane_updated_ = true;
  }
}


void Triangle3D::update() {
  updateCircumCenterIfNecessary();
  updatePlaneEquationIfNecessary();
}


int Triangle3D::orientationExact(const std::array<double, 3>& point_1,
                                    const std::array<double, 3>& point_2) const {
  auto points = getExactPositionVectors();
  auto normal_vector = points[1]->subtract(points[0])->crossProduct(points[2]->subtract(points[0]));
  auto offset = normal_vector->dotProduct(points[0]);
  return normal_vector->dotProduct(ExactVector::create(point_1))->compareTo(offset)
      * normal_vector->dotProduct(ExactVector::create(point_2))->compareTo(offset);
}


int Triangle3D::circleOrientation(const std::array<double, 3>& point) {
  updateCircumCenterIfNecessary();
  auto dummy = Matrix::subtract(point, circum_center_);
  double squared_distance = Matrix::dot(dummy, dummy);
  auto radial = Matrix::subtract(nodes_[0]->getPosition(), circum_center_);
  double squared_radius = Matrix::dot(radial, radial);
  double tolerance = squared_radius * 0.000000001;
  if (squared_distance < squared_radius + tolerance) {
    if (squared_distance > squared_radius - tolerance) {
      auto points = getExactPositionVectors();
      auto circum_center = calculateCircumCenterExact(points, calculateExactNormalVector(points));
      auto point_distance = circum_center->subtract(ExactVector::create(point))->squaredLength();
      auto squared_radius_x = circum_center->subtract(points[0])->squaredLength();
      return squared_radius_x->compareTo(point_distance);
    } else {
      return 1;
    }
  } else {
    return -1;
  }
}


std::shared_ptr<Tetrahedron> Triangle3D::getOppositeTetrahedron(
    const std::shared_ptr<Tetrahedron>& incident_tetrahedron) const {
  if (adjacent_tetrahedra_[0] == incident_tetrahedron) {
    return adjacent_tetrahedra_[1];
  } else if (adjacent_tetrahedra_[1] == incident_tetrahedron) {
    return adjacent_tetrahedra_[0];
  } else {
    std::cout << __FUNCTION__ << std::endl; //fnoexceptionthrow std::invalid_argument("Tetrahedron not known!");
  }
}


void Triangle3D::removeTetrahedron(const std::shared_ptr<Tetrahedron>& tetrahedron) {
  if (adjacent_tetrahedra_[0] == tetrahedron) {
    adjacent_tetrahedra_[0] = std::shared_ptr<Tetrahedron>(nullptr);
  } else {
    adjacent_tetrahedra_[1] = std::shared_ptr<Tetrahedron>(nullptr);
  }
}


bool Triangle3D::isOpenToSide(const std::array<double, 3>& point) {
  if (adjacent_tetrahedra_[0].get() == nullptr) {
    if (adjacent_tetrahedra_[1].get() == nullptr) {
      return true;
    } else {
      if (adjacent_tetrahedra_[1]->isInfinite()) {
        return true;
      }
      auto position = adjacent_tetrahedra_[1]->getOppositeNode(this->shared_from_this())
          ->getPosition();
      return !(this->onSameSide(position, point));
    }
  } else if (adjacent_tetrahedra_[1].get() == nullptr) {
    if (adjacent_tetrahedra_[0]->isInfinite()) {
      return true;
    }
    auto position =
        adjacent_tetrahedra_[0]->getOppositeNode(this->shared_from_this())->getPosition();
    return !(this->onSameSide(position, point));
  } else {
    return false;
  }
}


void Triangle3D::orientToSide(const std::array<double, 3>& position) {
  if (!isInfinite()) {
    updatePlaneEquationIfNecessary();
    double dot = Matrix::dot(position, this->normal_vector_);
    if (dot > this->offset_ + this->tolerance_) {
      upper_side_positive_ = true;
    } else if (dot < this->offset_ - this->tolerance_) {
      upper_side_positive_ = false;
    } else {
      auto points = getExactPositionVectors();
      auto normal_vector = calculateExactNormalVector(points);
      auto dot_1 = normal_vector->dotProduct(points[0]);
      auto dot_2 = normal_vector->dotProduct(ExactVector::create(position));
      int comparison = dot_1->compareTo(dot_2);
      if (comparison == 0) {
        //fnoexceptionthrow std::logic_error(
//            "The triangle cannot be oriented to because that point lies in the plane!");
      }
      upper_side_positive_ = comparison < 0;
    }
  }
}


void Triangle3D::orientToOpenSide() {
  if (!isInfinite()) {
    if (adjacent_tetrahedra_[0].get() == nullptr) {
      if (adjacent_tetrahedra_[1].get() == nullptr) {
        //fnoexceptionthrow std::logic_error("The triangle has two open sides!");
      }
      if (!adjacent_tetrahedra_[1]->isInfinite()) {
        orientToSide(
            adjacent_tetrahedra_[1]->getOppositeNode(this->shared_from_this())->getPosition());
        upper_side_positive_ ^= true;
      }
    } else if (adjacent_tetrahedra_[1].get() == nullptr) {
      if (!adjacent_tetrahedra_[0]->isInfinite()) {
        orientToSide(
            adjacent_tetrahedra_[0]->getOppositeNode(this->shared_from_this())->getPosition());
        upper_side_positive_ ^= true;
      }
    } else {
      //fnoexceptionthrow std::logic_error("The triangle has no open side!");
    }
  }
}


int Triangle3D::orientationToUpperSide(const std::array<double, 3>& point) const {
  double dot = Matrix::dot(point, this->normal_vector_);
  if (dot > this->offset_ + this->tolerance_) {
    return upper_side_positive_ ? 1 : -1;
  } else if (dot < this->offset_ - this->tolerance_) {
    return upper_side_positive_ ? -1 : 1;
  } else {
    auto points = getExactPositionVectors();
    auto normal_vector = calculateExactNormalVector(points);
    auto dot_1 = normal_vector->dotProduct(points[0]);
    auto dot_2 = normal_vector->dotProduct(ExactVector::create(point));
    if (dot_1->compareTo(dot_2) == 0) {
      return 0;
    } else {
      return ((dot_1->compareTo(dot_2) > 0) ^ upper_side_positive_) ? 1 : -1;
    }
  }
}


bool Triangle3D::onUpperSide(const std::array<double, 3>& point) const {
  return orientationToUpperSide(point) >= 0;
}


bool Triangle3D::trulyOnUpperSide(const std::array<double, 3>& point) const {
  return orientationToUpperSide(point) > 0;
}


double Triangle3D::getTypicalSDDistance() const {
  if (isInfinite()) {
    return std::numeric_limits<double>::max();
  } else {
    auto dummy = Matrix::subtract(nodes_[0]->getPosition(), circum_center_);
    return Matrix::norm(dummy) / Matrix::norm(this->normal_vector_);
  }
}


std::string Triangle3D::toString() const {
  return "T3D";
//  return "{(" + StringUtil::toStr(nodes_[0]) + "," + StringUtil::toStr(nodes_[1]) + ","
//      + StringUtil::toStr(nodes_[2]) + "), " + "(" + StringUtil::toStr(adjacent_tetrahedra_[0])
//      + "," + StringUtil::toStr(adjacent_tetrahedra_[1]) + "), " + StringUtil::toStr(circum_center_)
//      + ", " + StringUtil::toStr(plane_updated_) + ", " + StringUtil::toStr(circum_center_updated_)
//      + ", " + StringUtil::toStr(upper_side_positive_) + ", "
//      + StringUtil::toStr(connection_checked_) + ", " + StringUtil::toStr(this->normal_vector_)
//      + ", " + StringUtil::toStr(this->offset_) + ", " + StringUtil::toStr(this->tolerance_) + ", "
//      + StringUtil::toStr(this->normal_vector_updated_) + "}";
}


bool Triangle3D::isInfinite() const {
  return nodes_[0].get() == nullptr;
}


std::array<std::shared_ptr<SpaceNode >, 3> Triangle3D::getNodes() const {
  return nodes_;
}


void Triangle3D::addTetrahedron(const std::shared_ptr<Tetrahedron>& tetrahedron) {
  if (adjacent_tetrahedra_[0].get() == nullptr) {
    adjacent_tetrahedra_[0] = tetrahedron;
  } else {
    adjacent_tetrahedra_[1] = tetrahedron;
  }
  connection_checked_ = -1;
}


bool Triangle3D::wasCheckedAlready(int checking_index) {
  if (checking_index == connection_checked_) {
    return true;
  } else {
    connection_checked_ = checking_index;
    return false;
  }
}


bool Triangle3D::isAdjacentTo(const std::shared_ptr<Tetrahedron>& tetrahedron) const {
  return (adjacent_tetrahedra_[0] == tetrahedron) || (adjacent_tetrahedra_[1] == tetrahedron);
}


bool Triangle3D::isAdjacentTo(const std::shared_ptr<SpaceNode>& node) const {
  return (nodes_[0] == node) || (nodes_[1] == node) || (nodes_[2] == node);
}


bool Triangle3D::isCompletelyOpen() const {
  return (adjacent_tetrahedra_[0].get() == nullptr) && (adjacent_tetrahedra_[1].get() == nullptr);
}


bool Triangle3D::isClosed() const {
  return (adjacent_tetrahedra_[0].get() != nullptr) && (adjacent_tetrahedra_[1].get() != nullptr);
}


std::shared_ptr<ExactVector> Triangle3D::getExactNormalVector() const {
  return calculateExactNormalVector(getExactPositionVectors());
}


void Triangle3D::updateNormalVector(const std::array<double, 3>& new_normal_vector) {
  this->normal_vector_ = new_normal_vector;
  this->offset_ = Matrix::dot(this->normal_vector_, nodes_[0]->getPosition());
  this->normal_vector_updated_ = true;
}


std::shared_ptr<ExactVector> Triangle3D::calculateCircumCenterExact(
    const std::array<std::shared_ptr<ExactVector>, 3>& points,
    const std::shared_ptr<ExactVector>& normal_vector) {
  auto a = points[0];
  // Start by calculating the normal vectors:
  std::array<std::shared_ptr<ExactVector>, 3> n = { points[1]->subtract(a), points[2]->subtract(a),
      normal_vector };
  std::array<std::shared_ptr<Rational>, 3> rationals = { points[1]->add(a)->dotProduct(n[0])
      ->divideBy(Rational::create(2, 1)), points[2]->add(a)->dotProduct(n[1])->divideBy(
      Rational::create(2, 1)), a->dotProduct(n[2]) };
  return calculate3PlaneXPoint(n, rationals, ExactVector::det(n));
}


double Triangle3D::calculateSDDistance(const std::array<double, 3>& fourth_point) const {
  if (!isInfinite()) {
    // calc that distance within 6 subtractions, 3 additions, 1 division
    // and 9 multiplications. Beat that!
    auto ad = Matrix::subtract(nodes_[0]->getPosition(), fourth_point);
    double denominator = Matrix::dot(ad, this->normal_vector_);
    if ((denominator != 0.0) && (std::abs(denominator) < this->tolerance_)) {
      auto n0_vector = ExactVector::create(nodes_[0]->getPosition());
      auto v1 = n0_vector->subtract(ExactVector::create(nodes_[1]->getPosition()));
      auto v2 = n0_vector->subtract(ExactVector::create(nodes_[2]->getPosition()));
      auto normal_vector = v1->crossProduct(v2);
      auto dot = normal_vector->dotProduct(n0_vector->subtract(ExactVector::create(fourth_point)));
      if (dot->isZero()) {
        denominator = 0.0;
      } else {
        denominator = dot->doubleValue();
        dot = normal_vector->dotProduct(ExactVector::create(this->normal_vector_));
        if (dot->compareTo(Rational::create(0, 1)) < 0) {
          denominator = 0 - denominator;
        }
      }
    }
    if (denominator != 0) {
      double sd_distance = Matrix::dot(
          ad,
          Matrix::subtract(
              Matrix::scalarMult(0.5, Matrix::add(nodes_[0]->getPosition(), fourth_point)),
              circum_center_)) / denominator;
      return sd_distance;
    }
  }
  return std::numeric_limits<double>::max();
}


std::shared_ptr<Rational> Triangle3D::calculateSDDistanceExact(
    const std::array<std::shared_ptr<ExactVector>, 4>& points,
    const std::shared_ptr<ExactVector>& normal_vector) const {
  if (!isInfinite()) {
    // calc that distance within 6 subtractions, 3 additions, 1 division
    // and 9 multiplications. Beat that!
    auto ad = points[0]->subtract(points[3]);
    auto denominator = ad->dotProduct(normal_vector);
    if (!denominator->isZero()) {
      std::array<std::shared_ptr<ExactVector>, 3> points_3;
      for (std::size_t i = 0; i < 3; i++) {
        points_3[i] = points[i];
      }
      auto circum_center = calculateCircumCenterExact(points_3, normal_vector);
      return points[0]->add(points[3])->divideBy(Rational::create(2, 1))->decreaseBy(circum_center)
          ->dotProduct(ad)->divideBy(denominator);
    }
  }
  return Rational::create(std::numeric_limits<int64_t>::max());
}


void Triangle3D::updateCircumCenterIfNecessary() {
  if (!circum_center_updated_ && !isInfinite()) {
    circum_center_updated_ = true;
    auto a = nodes_[0]->getPosition();
    // Start by calculating the normal vectors:
    std::array<std::array<double, 3>, 3> n;
    auto line_1 = Matrix::subtract(nodes_[1]->getPosition(), a);
    auto line_2 = Matrix::subtract(nodes_[2]->getPosition(), a);
    n[0] = Matrix::normalize(line_1);
    n[1] = Matrix::normalize(line_2);
    n[2] = Matrix::crossProduct(n[0], n[1]);
    updateNormalVector(n[2]);
    this->normal_vector_updated_ = true;
    this->tolerance_ = Matrix::dot(this->normal_vector_, this->normal_vector_) * 0.000000001;
    this->normal_vector_updated_ = true;
    // cut the three planes:
    circum_center_ = calculate3PlaneXPoint(
        n,
        { Matrix::dot(Matrix::add(a, nodes_[1]->getPosition()), n[0]) * 0.5, Matrix::dot(
            Matrix::add(a, nodes_[2]->getPosition()), n[1]) * 0.5, Matrix::dot(a, n[2]) });
  }
}


std::array<std::shared_ptr<ExactVector>, 3> Triangle3D::getExactPositionVectors() const {
  std::array<std::shared_ptr<ExactVector>, 3> result;
  for (std::size_t i = 0; i < 3; i++) {
    result[i] = ExactVector::create(nodes_[i]->getPosition());
  }
  return result;
}


std::shared_ptr<ExactVector> Triangle3D::calculateExactNormalVector(
    const std::array<std::shared_ptr<ExactVector>, 3>& points) const {
  return points[1]->subtract(points[0])->crossProduct(points[2]->subtract(points[0]));
}

// define templates that should be compiled

}  // namespace spatial_organization
}  // namespace cx3d
