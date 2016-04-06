#include "spatial_organization/open_triangle_organizer.h"

#include <limits>
#include <vector>
#include <cstdlib>
#include <stdexcept>

#include "matrix.h"
#include "string_util.h"
#include "physics/physical_node.h"
#include "spatial_organization/rational.h"
#include "spatial_organization/exact_vector.h"
#include "spatial_organization/space_node.h"
#include "spatial_organization/triangle_3d.h"
#include "spatial_organization/tetrahedron.h"
#include "spatial_organization/edge_hash_key.h"
#include "spatial_organization/triangle_hash_key.h"

namespace cx3d {
namespace spatial_organization {


OpenTriangleOrganizer::OpenTriangleOrganizer(
    int preferred_capacity,
    const std::shared_ptr<SimpleTriangulationNodeOrganizer >& tno)
    : shortest_distance_ { std::numeric_limits<double>::max() },
      tno_ { tno } {
}


void OpenTriangleOrganizer::recoredNewTetrahedra() {
  //todo obsolete function in current implementation
}


std::vector<std::shared_ptr<Tetrahedron > > OpenTriangleOrganizer::getNewTetrahedra() {
  return new_tetrahedra_;
}


std::shared_ptr<Tetrahedron > OpenTriangleOrganizer::getANewTetrahedron() {
  return a_new_tetrahedron_;
}


void OpenTriangleOrganizer::removeAllTetrahedraInSphere(
    const std::shared_ptr<Tetrahedron >& starting_tetrahedron) {
  if (starting_tetrahedron->isValid()) {
    std::vector<std::shared_ptr<Tetrahedron>>tetrahedrons_to_remove;
    for(auto triangle : starting_tetrahedron->getAdjacentTriangles()) {
      auto opposite_tetrahedron = triangle->getOppositeTetrahedron(starting_tetrahedron);
      if ((opposite_tetrahedron.get() != nullptr)
          && (!starting_tetrahedron->isInfinite() ^ opposite_tetrahedron->isInfinite())
          && (starting_tetrahedron->isInsideSphere(opposite_tetrahedron->getOppositeNode(triangle)->getPosition()))) {
        tetrahedrons_to_remove.push_back(opposite_tetrahedron);
      }
      if (triangle->isClosed()) {
        putTriangle(triangle);
      } else {
        removeTriangle(triangle);
      }
    }
    starting_tetrahedron->remove();
    for(auto tetrahedron_remove : tetrahedrons_to_remove) {
      removeAllTetrahedraInSphere(tetrahedron_remove);
    }
  }
}


void OpenTriangleOrganizer::putTriangle(
    const std::shared_ptr<Triangle3D >& triangle) {
  auto nodes = triangle->getNodes();
  TriangleHashKey key(nodes[0], nodes[1], nodes[2]);
  map_[key] = triangle;
  tno_->addTriangleNodes(triangle);
  open_triangles_.push(triangle);
}


void OpenTriangleOrganizer::removeTriangle(
    const std::shared_ptr<Triangle3D >& triangle) {
  auto nodes = triangle->getNodes();
  TriangleHashKey key(nodes[0], nodes[1], nodes[2]);
  map_.erase(key);
}


std::shared_ptr<Triangle3D > OpenTriangleOrganizer::getTriangle(
    const std::shared_ptr<SpaceNode >& a,
    const std::shared_ptr<SpaceNode >& b,
    const std::shared_ptr<SpaceNode >& c) {
  TriangleHashKey key(a, b, c);
  std::shared_ptr<Triangle3D> ret(nullptr);
  if (map_.find(key) == map_.end()) {
    std::shared_ptr<Tetrahedron> null_tetraherdon(nullptr);
    ret = Triangle3D::create(a, b, c, null_tetraherdon, null_tetraherdon);
    map_[key] = ret;
    open_triangles_.push(ret);
  } else {
    ret = map_[key];
    if (ret->isCompletelyOpen()) {
      open_triangles_.push(ret);
    } else {
      map_.erase(key);
    }
  }
  return ret;
}


std::shared_ptr<Triangle3D > OpenTriangleOrganizer::getTriangleWithoutRemoving(
    const std::shared_ptr<SpaceNode >& a,
    const std::shared_ptr<SpaceNode >& b,
    const std::shared_ptr<SpaceNode >& c) {
  TriangleHashKey key(a, b, c);
  std::shared_ptr<Triangle3D> ret(nullptr);
  if (map_.find(key) == map_.end()) {
    std::shared_ptr<Tetrahedron> null_tetraherdon(nullptr);
    ret = Triangle3D::create(a, b, c, null_tetraherdon, null_tetraherdon);
    map_[key] = ret;
    open_triangles_.push(ret);
  } else {
    ret = map_[key];
  }
  return ret;
}


void OpenTriangleOrganizer::triangulate() {
  if (open_triangles_.empty())
    createInitialTriangle();
  double upper_bound = 0, lower_bound = 0;
  std::deque<std::shared_ptr<SpaceNode>>similar_distance_nodes;
  std::deque<std::shared_ptr<SpaceNode>>onCircle_nodes;
  auto open_triangle = getAnOpenTriangle();
  int security_counter = 0;
  while (open_triangle.get() != nullptr) {
    open_triangle->update();
    std::shared_ptr<SpaceNode> picked_node(nullptr);
    open_triangle->orientToOpenSide();
    auto last_tetrahedron = open_triangle->getOppositeTetrahedron(
        std::shared_ptr<Tetrahedron>(nullptr));
    shortest_distance_ = std::numeric_limits<double>::max();
    upper_bound = shortest_distance_;
    lower_bound = shortest_distance_;
    double tolerance = open_triangle->getTypicalSDDistance() * 0.0000001;
    for (auto node : tno_->getNodes(open_triangle->getNodes()[0])) {
      if (!open_triangle->isAdjacentTo(node)) {
        double current_distance = open_triangle->getSDDistance(
            node->getPosition());
        if ((current_distance < upper_bound)) {
          bool smaller = false;
          if (current_distance > lower_bound) {
            auto last_sd_distance = open_triangle->getSDDistanceExact(
                picked_node->getPosition());
            auto new_sd_distance = open_triangle->getSDDistanceExact(
                node->getPosition());
            int comparison = last_sd_distance->compareTo(new_sd_distance);
            if (comparison == 0) {
              similar_distance_nodes.push_back(node);
            } else if (comparison > 0) {
              smaller = true;
            }
          } else {
            smaller = true;
          }
          if (smaller) {
            similar_distance_nodes.clear();
            shortest_distance_ = current_distance;
            // bounds to determine if another node causes the 'same' signed delaunay distance:
            upper_bound = shortest_distance_ + tolerance;
            lower_bound = shortest_distance_ - tolerance;
            picked_node = node;
          }
        } else if (open_triangle->orientationToUpperSide(node->getPosition())
            == 0
            && open_triangle->circleOrientation(node->getPosition()) == 0) {
          onCircle_nodes.push_back(node);
        }
      }
    }
    if (picked_node.get() == nullptr
        || (similar_distance_nodes.empty() && onCircle_nodes.empty())) {
      createNewTetrahedron(open_triangle, picked_node);
    } else {
      similar_distance_nodes.push_back(picked_node);
      triangulatePointsOnSphere(similar_distance_nodes, onCircle_nodes,
                                open_triangle);
    }
    similar_distance_nodes.clear();
    onCircle_nodes.clear();
    open_triangle = getAnOpenTriangle();
    security_counter++;
    if (security_counter > 2000) {
      //fnoexceptionthrow std::runtime_error("Am I in an infinite loop?");
    }
  }
}


std::string OpenTriangleOrganizer::toString() const {
  std::ostringstream str_stream;
  str_stream << "(";
  str_stream << StringUtil::toStr(static_cast<int>(map_.size()));
  str_stream << ", ";
  str_stream << StringUtil::toStr(new_tetrahedra_);
  str_stream << ", ";
  str_stream << StringUtil::toStr(tno_);
  str_stream << ", ";
  str_stream << StringUtil::toStr(static_cast<int>(open_triangles_.size()));
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(shortest_distance_);
  str_stream << ", ";
  str_stream << StringUtil::toStr(a_new_tetrahedron_);
  str_stream << ")";
  return str_stream.str();
}


bool OpenTriangleOrganizer::equalTo(
    const std::shared_ptr<OpenTriangleOrganizer>& other) const {
  return this == other.get();
}


bool OpenTriangleOrganizer::contains(
    const std::shared_ptr<SpaceNode >& a,
    const std::shared_ptr<SpaceNode >& b,
    const std::shared_ptr<SpaceNode >& c) const {
  TriangleHashKey key(a, b, c);
  return map_.find(key) != map_.end();
}


bool OpenTriangleOrganizer::isEmpty() const {
  map_.empty();
}


std::shared_ptr<Triangle3D > OpenTriangleOrganizer::getAnOpenTriangle() {
  if (open_triangles_.empty()) {
    return std::shared_ptr<Triangle3D>(nullptr);
  }
  auto ret = open_triangles_.top();
  open_triangles_.pop();
  while ((ret->isInfinite() || ret->isClosed() || ret->isCompletelyOpen())) {
    if (open_triangles_.empty()) {
      return std::shared_ptr<Triangle3D>(nullptr);
    }
    ret = open_triangles_.top();
    open_triangles_.pop();
  }
  return ret;
}


std::shared_ptr<EdgeHashKey> OpenTriangleOrganizer::putEdgeOnMap(
    const std::shared_ptr<SpaceNode>& a,
    const std::shared_ptr<SpaceNode>& b,
    const std::shared_ptr<SpaceNode>& opposite_node,
    const std::shared_ptr<EdgeHashKey>& old_open_edge,
    std::unordered_map<EdgeHashKey, std::shared_ptr<EdgeHashKey>,
        EdgeHashKeyHash, EdgeHashKeyEqual >& map) {
  auto hk1 = std::shared_ptr<EdgeHashKey>(
      new EdgeHashKey(a, b, opposite_node));
  if (map.find(*hk1) != map.end()) {
    map.erase(*hk1);
    return old_open_edge;
  } else {
    map[*hk1] = hk1;
    return hk1;
  }
}


std::shared_ptr<SpaceNode > OpenTriangleOrganizer::findCenterNode(
    const std::deque<std::shared_ptr<SpaceNode > >& nodes) {
  std::shared_ptr<SpaceNode> center_node(nullptr);
  int minID = std::numeric_limits<int>::max();
  for (auto node : nodes) {
    if (node->getId() < minID) {
      minID = node->getId();
      center_node = node;
    }
  }
  return center_node;
}


std::shared_ptr<EdgeHashKey> OpenTriangleOrganizer::triangulateSortedCirclePoints(
    const std::vector<std::shared_ptr<SpaceNode > >& sorted_nodes,
    const std::shared_ptr<SpaceNode >& center_node,
    std::unordered_map<EdgeHashKey, std::shared_ptr<EdgeHashKey>,
        EdgeHashKeyHash, EdgeHashKeyEqual >& map,
    std::vector<std::shared_ptr<Triangle3D > >& triangle_list) {
  auto it = sorted_nodes.begin();
  auto last = *it;
  it++;
  auto current = *it;
  std::shared_ptr<EdgeHashKey> ret_value;
  while (it != sorted_nodes.end()) {
    last = current;
    it++;
    current = *it;
    auto triangle = getTriangleWithoutRemoving(last, current, center_node);
    triangle_list.push_back(triangle);
    putEdgeOnMap(center_node, last, current, nullptr, map);
    ret_value = putEdgeOnMap(last, current, center_node, ret_value, map);
    putEdgeOnMap(current, center_node, last, nullptr, map);
  }
  return ret_value;
}


void OpenTriangleOrganizer::removeForbiddenTriangles(
    const std::vector<std::shared_ptr<SpaceNode > >& sorted_nodes) {  //todo list is inefficient here - use vector
  std::shared_ptr<Tetrahedron> null_tetrahedron(nullptr);
  // Special treatment for situation with 4 nodes only:
  if (sorted_nodes.size() == 4) {
    auto it = sorted_nodes.begin();
    auto center = *it++;
    auto a = *it++;
    auto b = *it++;
    auto c = *it++;
    // in case there is only one valid triangle, remove the others:
    if (contains(center, a, b)) {
      if (!contains(center, b, c)) {
        auto tetrahedron = getTriangleWithoutRemoving(center, a, b)
            ->getOppositeTetrahedron(null_tetrahedron);
        removeAllTetrahedraInSphere(tetrahedron);
      }
    } else if (contains(center, b, c)) {
      auto tetrahedron = getTriangleWithoutRemoving(center, b, c)
          ->getOppositeTetrahedron(null_tetrahedron);
      removeAllTetrahedraInSphere(tetrahedron);
    }
    // otherwise, check if there are triangles that are not allowed (and remove them if necessary
    else {
      if (contains(a, b, c))
        removeAllTetrahedraInSphere(
            getTriangleWithoutRemoving(a, b, c)->getOppositeTetrahedron(
                null_tetrahedron));
      if (contains(center, a, c))
        removeAllTetrahedraInSphere(
            getTriangleWithoutRemoving(center, a, c)->getOppositeTetrahedron(
                null_tetrahedron));
    }
  }
  // general case:
  else {
    bool remove_all_circle_triangles = false;
    // first, copy nodes to array for faster access:
    std::vector<std::shared_ptr<SpaceNode>>nodes(sorted_nodes.begin(), sorted_nodes.end());  //todo inefficient
    // check if any valid triangle is missing, if yes, set removeAllCircleTriangles = true
    for (int i = 1; (i < nodes.size() - 1) && (!remove_all_circle_triangles);
        i++) {
      if (!contains(nodes[0], nodes[i], nodes[i + 1])) {
        remove_all_circle_triangles = true;
      }
    }
    // if any valid triangle was missing...
    if (remove_all_circle_triangles) {
      // ...test if any triangle imaginable exists...
      for (int i = 0; i < nodes.size() - 2; i++) {
        for (int j = i + 1; j < nodes.size() - 1; j++) {
          for (int k = j + 1; k < nodes.size(); k++) {
            if (contains(nodes[i], nodes[j], nodes[k])) {
              // and remove it together with incident tetrahedra:
              auto tetrahedron = getTriangleWithoutRemoving(nodes[i], nodes[j],
                                                            nodes[k])
                  ->getOppositeTetrahedron(null_tetrahedron);
              removeAllTetrahedraInSphere(tetrahedron);
            }
          }
        }
      }
    }
  }
}


std::vector<std::shared_ptr<SpaceNode > > OpenTriangleOrganizer::sortCircleNodes(
    std::deque<std::shared_ptr<SpaceNode > >& nodes,
    std::shared_ptr<EdgeHashKey> starting_edge,
    const std::shared_ptr<SpaceNode >& center_node) {
  std::deque<std::shared_ptr<SpaceNode>>sorted_nodes;
  std::shared_ptr<SpaceNode> null_space_node(nullptr);
  auto search_node = null_space_node;
  auto last_search_node = null_space_node;
  auto removed_node_1 = null_space_node;
  auto removed_node_2 = null_space_node;
  if (starting_edge.get() == nullptr) {
    last_search_node = nodes.front();
    nodes.pop_front();
    double min_distance = std::numeric_limits<double>::max();
    for (auto node : nodes) {
      auto dummy = Matrix::subtract(last_search_node->getPosition(), node->getPosition());
      double dot = Matrix::dot(dummy, dummy);
      if (dot < min_distance) {
        search_node = node;
        min_distance = dot;
      }
    }
    STLUtil::dequeRemove(nodes, search_node);
    removed_node_1 = last_search_node;
    removed_node_2 = search_node;
  } else {
    search_node = starting_edge->getEndpointB();
    last_search_node = starting_edge->getEndpointA();
  }
  while (!nodes.empty()) {
    auto last_vector = Matrix::normalize(Matrix::subtract(search_node->getPosition(), last_search_node->getPosition()));
    double biggest_cosinus = -2.0;
    auto picked_node = null_space_node;
    for (auto node : nodes) {
      auto dummy = Matrix::normalize(Matrix::subtract(node->getPosition(), search_node->getPosition()));
      double current_cosinus = Matrix::dot(dummy, last_vector);
      if (current_cosinus > biggest_cosinus) {
        biggest_cosinus = current_cosinus;
        picked_node = node;
      }
    }
    sorted_nodes.push_back(picked_node);
    last_search_node = search_node;
    search_node = picked_node;
    STLUtil::dequeRemove(nodes, picked_node);
  }
  if (starting_edge.get() != nullptr) {
    sorted_nodes.push_front(starting_edge->getEndpointB());
    sorted_nodes.push_front(starting_edge->getEndpointA());
  } else {
    sorted_nodes.push_front(removed_node_2);
    sorted_nodes.push_front(removed_node_1);
  }
  while (!sorted_nodes.empty() && sorted_nodes.front() != center_node) {
    auto node = sorted_nodes.front();
    sorted_nodes.pop_front();
    nodes.push_back(node);
  }
  for(auto node : nodes) {
    sorted_nodes.push_back(node);
  }
  //fixme IW cpy
  std::vector<std::shared_ptr<SpaceNode>>ret;
  for(auto el : sorted_nodes) {
    ret.push_back(el);
  }
  return ret;
}


std::shared_ptr<EdgeHashKey> OpenTriangleOrganizer::triangulatePointsOnCircle(
    std::deque<std::shared_ptr<SpaceNode > >& similar_distance_nodes,
    const std::shared_ptr<EdgeHashKey>& starting_edge,
    std::unordered_map<EdgeHashKey, std::shared_ptr<EdgeHashKey>,
        EdgeHashKeyHash, EdgeHashKeyEqual >& map,
    std::vector<std::shared_ptr<Triangle3D > >& triangle_list) {
  if (starting_edge.get() != nullptr) {
    similar_distance_nodes.push_front(starting_edge->getEndpointA());
    similar_distance_nodes.push_front(starting_edge->getEndpointB());
  }
  auto center_node = findCenterNode(similar_distance_nodes);
  if (starting_edge.get() != nullptr) {
    similar_distance_nodes.pop_front();
    similar_distance_nodes.pop_front();
  }
  auto sorted_nodes = sortCircleNodes(similar_distance_nodes, starting_edge,
                                      center_node);
  removeForbiddenTriangles(sorted_nodes);
  return triangulateSortedCirclePoints(sorted_nodes, center_node, map,
                                       triangle_list);
}


void OpenTriangleOrganizer::triangulatePointsOnSphere(
    std::deque<std::shared_ptr<SpaceNode > >& nodes,
    std::deque<std::shared_ptr<SpaceNode > >& on_circle_nodes,
    const std::shared_ptr<Triangle3D >& startingTriangle) {
  std::vector<std::shared_ptr<Triangle3D>>surface_triangles;
  auto starting_triangle_nodes = startingTriangle->getNodes();
  nodes.push_back(starting_triangle_nodes[0]);
  nodes.push_back(starting_triangle_nodes[1]);
  nodes.push_back(starting_triangle_nodes[2]);
  for(auto node : on_circle_nodes) {
    nodes.push_back(node);
  }
  std::unordered_map<EdgeHashKey, std::shared_ptr<EdgeHashKey>, EdgeHashKeyHash, EdgeHashKeyEqual > map;
  std::shared_ptr<EdgeHashKey> an_open_edge(nullptr);
  if (on_circle_nodes.empty()) {
    surface_triangles.push_back(startingTriangle);
    for (auto i = 0; i < 3; i++) {
      an_open_edge = putEdgeOnMap(
          starting_triangle_nodes[i],
          starting_triangle_nodes[(i + 1) % 3],
          starting_triangle_nodes[(i + 2) % 3],
          an_open_edge, map);
    }
  } else {
    on_circle_nodes.push_back(starting_triangle_nodes[0]);
    on_circle_nodes.push_back(starting_triangle_nodes[1]);
    on_circle_nodes.push_back(starting_triangle_nodes[2]);
    an_open_edge = triangulatePointsOnCircle(on_circle_nodes, std::shared_ptr<EdgeHashKey>(nullptr), map, surface_triangles);
  }
  std::deque<std::shared_ptr<SpaceNode>> similar_distance_nodes;
  double upper_bound, lower_bound;
  while (!map.empty()) {
    auto a = an_open_edge->getEndpointA(), b = an_open_edge->getEndpointB();
    double smallest_cosinus = std::numeric_limits<double>::max();
    upper_bound = smallest_cosinus;
    lower_bound = smallest_cosinus;
    std::shared_ptr<SpaceNode> picked_node(nullptr);
    double tolerance = 0.000000001;
    for (auto current_node : nodes) {
      if (current_node.get() != an_open_edge->getEndpointA().get()
          && current_node.get() != an_open_edge->getEndpointB().get()) {
        double cosinus = an_open_edge->getCosine(current_node->getPosition());
        if (cosinus < upper_bound) {
          if (cosinus > lower_bound) {
            similar_distance_nodes.push_back(current_node);
          } else {
            picked_node = current_node;
            smallest_cosinus = cosinus;
            upper_bound = smallest_cosinus + tolerance;
            lower_bound = smallest_cosinus - tolerance;
            similar_distance_nodes.clear();
          }
        }
      }
    }
    if (similar_distance_nodes.empty()) {
      auto new_triangle = getTriangleWithoutRemoving( a, b, picked_node);
      surface_triangles.push_back(new_triangle);
      // add the new edges to the hashmap:
      map.erase(*an_open_edge);
      an_open_edge = putEdgeOnMap(a, picked_node, b, std::shared_ptr<EdgeHashKey>(nullptr), map);
      an_open_edge = putEdgeOnMap(b, picked_node, a, an_open_edge, map);

    } else {
      similar_distance_nodes.push_back(picked_node);
      an_open_edge = triangulatePointsOnCircle( similar_distance_nodes, an_open_edge, map, surface_triangles);
      similar_distance_nodes.clear();
    }
    if ((an_open_edge.get() == nullptr) && (!map.empty())) {
      auto it = map.begin()++;
      an_open_edge = (*it).second;
    }

  }
  auto center_node = findCenterNode(nodes);
  for (auto triangle : surface_triangles) {
    if (!triangle->isAdjacentTo(center_node))
    createNewTetrahedron(triangle, center_node);
  }
}


std::shared_ptr<Rational> OpenTriangleOrganizer::calc2DSDDistanceExact(
    const std::array<double, 3>& av, const std::array<double, 3>& bv,
    const std::array<double, 3>& third_point) {
  auto av_x = ExactVector::create(av);
  auto bv_y = ExactVector::create(bv);
  auto third_point_x = ExactVector::create(third_point);
  auto av_to_third_point_x = third_point_x->subtract(av_x);
  std::array<std::shared_ptr<ExactVector>, 3> normals_x = { bv_y->subtract(
      av_x), bv_y->subtract(av_x)->crossProduct(av_to_third_point_x),
      av_to_third_point_x };

  std::array<std::shared_ptr<Rational>, 3> offsets_x = { normals_x[0]
      ->dotProduct(av_x->add(bv_y))->multiply(Rational::create(1, 2)),
      normals_x[1]->dotProduct(av_x), normals_x[2]->dotProduct(
          av_x->add(third_point_x))->multiply(Rational::create(1, 2)) };
  auto circum_center = Triangle3D::calculate3PlaneXPoint(
      normals_x, offsets_x, ExactVector::det(normals_x));
  return circum_center->subtract(
      av_x->add(bv_y)->multiply(Rational::create(1, 2)))->squaredLength();
}


void OpenTriangleOrganizer::createInitialTriangle() {
  // find a starting node:
  auto a = tno_->getFirstNode();

  double tolerance = 0.000000001;  //todo bad design
  // find the second node by minimizing the distance to the first node:
  shortest_distance_ = std::numeric_limits<double>::max();
  std::shared_ptr<SpaceNode> b(nullptr);
  for (auto dummy : tno_->getNodes(a)) {
    auto vector = Matrix::subtract(dummy->getPosition(), a->getPosition());
    double distance = Matrix::dot(vector, vector);
    if (distance < shortest_distance_ + tolerance) {
      if (distance > shortest_distance_ - tolerance) {
        auto dist_new = (ExactVector::create(a->getPosition()))->subtract(
            ExactVector::create(dummy->getPosition()))->squaredLength();
        auto dist_last = (ExactVector::create(a->getPosition()))->subtract(
            ExactVector::create(b->getPosition()))->squaredLength();
        if (dist_last->compareTo(dist_new) > 0) {
          b = dummy;
          shortest_distance_ = std::min(shortest_distance_, distance);
        }
      } else {
        b = dummy;
        shortest_distance_ = distance;
        tolerance = 0.000000001 * distance;
      }
    }
  }

  // find the third node by minimizing the distance between the center of
  // the circumcircle and the middle point between a and b:

  shortest_distance_ = std::numeric_limits<double>::max();
  auto av = a->getPosition();
  auto bv = b->getPosition();
  std::array<std::array<double, 3>, 3> normals;
  normals[0] = Matrix::subtract(bv, av);
  std::array<double, 3> offsets;
  offsets[0] = 0.5 * Matrix::dot(normals[0], Matrix::add(av, bv));
  std::shared_ptr<SpaceNode> c(nullptr);
  tolerance = Matrix::dot(normals[0], normals[0]) * 0.000000001;

  for (auto dummy : tno_->getNodes(a)) {
    auto dummy_pos = dummy->getPosition();
    auto av_to_dummy_pos = Matrix::subtract(dummy_pos, av);
    normals[1] = Matrix::crossProduct(normals[0], av_to_dummy_pos);
    offsets[1] = Matrix::dot(normals[1], av);
    normals[2] = av_to_dummy_pos;
    offsets[2] = 0.5 * Matrix::dot(normals[2], Matrix::add(av, dummy_pos));
    // find the circumcenter by cutting 3 planes:
    // the plane describing all points with equal distance to a and b,
    // the plane defined by a, b and dummy and
    // the plane describing all points with equal distance to a and
    // dummy
    auto circum_center = Triangle3D::calculate3PlaneXPoint(normals, offsets);
    auto vector = Matrix::subtract(
        circum_center, Matrix::scalarMult(0.5, Matrix::add(av, bv)));
    double distance = Matrix::dot(vector, vector);
    if (distance < shortest_distance_ + tolerance) {
      if (distance > shortest_distance_ - tolerance) {
        auto dist_1 = calc2DSDDistanceExact(av, bv, dummy_pos);
        auto dist_2 = calc2DSDDistanceExact(av, bv, c->getPosition());
        int comparison = dist_1->compareTo(dist_2);
        if ((comparison < 0)
            || ((comparison == 0) && (dummy->getId() < c->getId()))) {
          c = dummy;
          shortest_distance_ = std::min(shortest_distance_, distance);
        }
      } else {
        c = dummy;
        shortest_distance_ = distance;
      }
    }
  }
  std::shared_ptr<Tetrahedron> null_tetrahedron(nullptr);
  putTriangle(
      Triangle3D::create(a, b, c, null_tetrahedron, null_tetrahedron));
}


void OpenTriangleOrganizer::createNewTetrahedron(
    const std::shared_ptr<Triangle3D >& open_triangle,
    const std::shared_ptr<SpaceNode >& opposite_node) {
  a_new_tetrahedron_ = Tetrahedron::create(open_triangle, opposite_node,
                                              this->shared_from_this());
  new_tetrahedra_.push_back(a_new_tetrahedron_);
}

}  // namespace cx3d
}  // namespace spatial_organization
