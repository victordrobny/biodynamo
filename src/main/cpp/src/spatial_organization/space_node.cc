#include "spatial_organization/space_node.h"

#include <cmath>
#include <limits>
#include <sstream>

#include "physics/physical_node.h"
#include "matrix.h"
#include "stl_util.h"
#include "string_util.h"
#include "java_util.h"
#include "physics/physical_node.h"
#include "sim_state_serialization_util.h"
#include "spatial_organization/edge.h"
#include "spatial_organization/triangle_3d.h"
#include "spatial_organization/tetrahedron.h"
#include "spatial_organization/open_triangle_organizer.h"
#include "spatial_organization/spatial_organization_edge.h"
#include "spatial_organization/spatial_organization_node_movement_listener.h"

namespace cx3d {
namespace spatial_organization {


JavaUtil* SpaceNode::java_ = nullptr;


int SpaceNode::checking_index_ = 0;


int SpaceNode::id_counter_ = 0;


Tetrahedron* SpaceNode::searchInitialInsertionTetrahedron(
    Tetrahedron* start,
    const std::array<double, 3>& coordinate) {
  Tetrahedron* current = start;
  SpaceNode* null_space_node = nullptr;
  Tetrahedron* null_tetrahedron = nullptr;
  if (current->isInfinite()) {
    current = current->getOppositeTriangle(null_space_node)
        ->getOppositeTetrahedron(current);
  }
  auto last = null_tetrahedron;
  while (current != last && !current->isInfinite()) {
    last = current;
    current = current->walkToPoint(coordinate, java_->generateTriangleOrder());
  }
  return current;
}


SpaceNode::SpaceNode()
    : id_(SpaceNode::id_counter_++),
      content_(),
      listeners_(),
      position_({ 0.0, 0.0, 0.0 }),
      adjacent_edges_(),
      adjacent_tetrahedra_(),
      volume_(0.0) {
}


SpaceNode::SpaceNode(const std::array<double, 3>& position,
                        physics::PhysicalNode* content)
    : id_(SpaceNode::id_counter_++),
      content_(content),
      listeners_(),
      position_(position),
      adjacent_edges_(),
      adjacent_tetrahedra_(),
      volume_(0.0) {
}


SpaceNode::SpaceNode(double x, double y, double z,
                        physics::PhysicalNode* content)
    : id_(SpaceNode::id_counter_++),
      content_(content),
      listeners_(),
      position_({ x, y, z }),
      adjacent_edges_(),
      adjacent_tetrahedra_(),
      volume_(0.0) {
}


void SpaceNode::addSpatialOrganizationNodeMovementListener(
    SpatialOrganizationNodeMovementListener* listener) {
  listeners_.push_back(listener);
}


std::vector<Edge* > SpaceNode::getEdges() const {
  return adjacent_edges_;
}


std::vector<physics::PhysicalNode* > SpaceNode::getNeighbors() const {
  std::vector<physics::PhysicalNode*> result;
  for (auto e : adjacent_edges_) {
    result.push_back(e->getOpposite(this)->getUserObject());
  }
  return result;
}


SpaceNode* SpaceNode::getNewInstance(
    const std::array<double, 3>& position,
    physics::PhysicalNode* user_object) {
  // create a new SpaceNode:
  auto insert_point = SpaceNode::create(position, user_object);

  // the new instance should have the same listeners!
  insert_point->setListenerList(listeners_);
  // check if this point is capable of inserting a new point:
  if (adjacent_tetrahedra_.empty()) {
    // enough nodes collected:
    if (adjacent_edges_.size() == 2) {
      // collect the nodes:
      auto a = adjacent_edges_.front()->getOpposite(this);
      auto b = adjacent_edges_.back()->getOpposite(this);
      // clear the edge lists:
      adjacent_edges_.clear();
      a->adjacent_edges_.clear();
      b->adjacent_edges_.clear();
      // now create the first tetrahedron:
      auto oto = java_->oto_createSimpleOpenTriangleOrganizer();
      Tetrahedron::createInitialTetrahedron(this,
                                               insert_point, a, b, oto);
    } else {
      Edge::create(this, insert_point);
      if (adjacent_edges_.size() == 2) {
        Edge::create(
            adjacent_edges_.back()->getOpposite(this),
            insert_point);
      }
    }
  } else {
    // insert point:
    insert_point->insert(adjacent_tetrahedra_.front());
  }
  return insert_point;
}

std::vector<physics::PhysicalNode* > SpaceNode::getPermanentListOfNeighbors() const {
  std::vector<physics::PhysicalNode* > ret;
  for (auto e : adjacent_edges_) {
    auto opp = e->getOpposite(this);
    if (opp != nullptr) {
      ret.push_back(opp->getUserObject());
    }
  }
  return ret;
}


std::array<double, 3> SpaceNode::getPosition() const {
  return position_;
}


physics::PhysicalNode* SpaceNode::getUserObject() const {
  return content_;
}


std::array<physics::PhysicalNode*, 4> SpaceNode::getVerticesOfTheTetrahedronContaining(
    const std::array<double, 3>& position,
    std::array<int, 1>& returned_null) const {
  returned_null[0] = 0;
  if (adjacent_tetrahedra_.empty()) {
    returned_null[0] = 1;
    return std::array<physics::PhysicalNode*, 4>();
  }
  auto insertion_tetrahedron = adjacent_tetrahedra_.front();
  if (insertion_tetrahedron->isInfinite()) {
    auto opposite_triangle = insertion_tetrahedron->getOppositeTriangle(nullptr);
    insertion_tetrahedron = opposite_triangle->getOppositeTetrahedron(
        insertion_tetrahedron);
  }
  Tetrahedron* last(nullptr);
  while (!insertion_tetrahedron->equalTo(last)
      && (!insertion_tetrahedron->isInfinite())) {
    last = insertion_tetrahedron;
    // TODO walktoPoint doesn't //fnoexceptionthrow exception
    //    try {
    std::array<int, 4> triangleOrder = java_->generateTriangleOrder();
    insertion_tetrahedron = insertion_tetrahedron->walkToPoint(position,
                                                               triangleOrder);
//    } catch (PositionNotAllowedException e) {
//      insertionTetrahedron = last;
//    }
  }
  if (insertion_tetrahedron->isInfinite()) {
    returned_null[0] = 1;
    return std::array<physics::PhysicalNode*, 4>();
  }
  std::array<physics::PhysicalNode*, 4> ret;
  auto cnt = 0;
  for (auto node : insertion_tetrahedron->getAdjacentNodes()) {
    if (node != nullptr) {
      ret[cnt++] = node->getUserObject();
    }
  }
  return ret;
}


double SpaceNode::getVolume() const {
  return volume_;
}


void SpaceNode::moveFrom(const std::array<double, 3>& delta) {
  moveTo(Matrix::add(position_, delta));
}


void SpaceNode::remove() {
  removeAndReturnCreatedTetrahedron();
}


StringBuilder& SpaceNode::simStateToJson(StringBuilder& sb) const {
  sb.append("{");

  SimStateSerializationUtil::keyValue(sb, "id", id_);
  SimStateSerializationUtil::keyValue(sb, "position", position_);
  SimStateSerializationUtil::keyValue(sb, "volume", volume_);

  SimStateSerializationUtil::removeLastChar(sb);
  sb.append("}");
  return sb;
}


std::vector<Tetrahedron* > SpaceNode::getAdjacentTetrahedra() const {
  return adjacent_tetrahedra_;
}


void SpaceNode::addAdjacentTetrahedron(
    Tetrahedron* tetrahedron) {
  adjacent_tetrahedra_.push_back(tetrahedron);
}


void SpaceNode::removeTetrahedron(
    Tetrahedron* tetrahedron) {
  STLUtil::vectorRemove(adjacent_tetrahedra_, tetrahedron);
}


void SpaceNode::moveTo(const std::array<double, 3>& new_position) {
  if (checkIfTriangulationIsStillValid(new_position)) {
    auto delta = Matrix::subtract(new_position, position_);
    for (auto listener : listeners_) {
      listener->nodeAboutToMove(this, delta);
    }
    position_ = new_position;
    restoreDelaunay();
    for (auto listener : listeners_) {
      listener->nodeMoved(this);
    }
  } else {
    auto insert_position = SpaceNode::searchInitialInsertionTetrahedron(
        adjacent_tetrahedra_.front(), new_position);
    auto aNewTetrahedron = removeAndReturnCreatedTetrahedron();
    if (!insert_position->isValid()) {
      insert_position = aNewTetrahedron;
    }
    auto oldPosition = position_;
    position_ = new_position;
    // TODO positionnotallowedexception
//    try {
    insert(insert_position);
//    } catch (PositionNotAllowedException e) {
//       position_ = oldPosition;
//       insert(insert_position);
//      //fnoexceptionthrow e;
  }
}


void SpaceNode::changeVolume(double change) {
  volume_ += change;
}


void SpaceNode::addEdge(Edge* edge) {
  adjacent_edges_.push_back(edge); //fixme IW breaks
}


int SpaceNode::getId() const {
  return id_;
}


Edge* SpaceNode::searchEdge(
    SpaceNode* opposite_node) {
  for (auto e : adjacent_edges_) {
    auto opp = e->getOpposite(this);
    if ((opp != nullptr && opp->equalTo(opposite_node))
        || (opposite_node != nullptr && opposite_node->equalTo(opp))) {
      return e;
    }
  }
  return Edge::create(this, opposite_node);
}


void SpaceNode::removeEdge(Edge* edge) {
  STLUtil::vectorRemove(adjacent_edges_, edge);
}


void SpaceNode::setListenerList(
    const std::vector<SpatialOrganizationNodeMovementListener* >& listeners) {
  listeners_ = listeners;
}


Tetrahedron* SpaceNode::searchInitialInsertionTetrahedron(
    Tetrahedron* start) {
  return SpaceNode::searchInitialInsertionTetrahedron(start, position_);
}


Tetrahedron* SpaceNode::insert(
    Tetrahedron* start) {
  auto insertion_start = searchInitialInsertionTetrahedron(start);

  if (listeners_.size() != 0) {
    // tell the listeners that there will be a new node:
    auto vertice_contents = insertion_start->getVerticeContents();
    for (auto listener : listeners_) {
      listener->nodeAboutToBeAdded(this, position_,
                                   vertice_contents);
    }
  }

  auto oto = java_->oto_createSimpleOpenTriangleOrganizer();
  std::vector<Triangle3D* > queue;
  std::vector<Triangle3D* > outer_triangles;

  processTetrahedron(insertion_start, queue, oto);
//  while (!queue.empty()) {
  for(auto i = 0; i < queue.size(); i++) {
    auto currentTriangle = queue[i];
//    queue.pop_front();
    auto oppositeTetrahedron = currentTriangle->getOppositeTetrahedron(
        nullptr);
    if ((oppositeTetrahedron != nullptr)) {
      if (oppositeTetrahedron->isTrulyInsideSphere(position_)) {
        processTetrahedron(oppositeTetrahedron, queue, oto);
      } else {
        outer_triangles.push_back(currentTriangle);
      }
    }
  }
  Tetrahedron* ret(nullptr);
  // create a star-shaped triangulation:
  for (auto current_triangle : outer_triangles) {
    if (!current_triangle->isCompletelyOpen()) {
      ret = Tetrahedron::create(current_triangle, this,
                                   oto);
    }
  }
  // }
  // tell the listeners that the node was added:
  for (auto listener : listeners_) {
    listener->nodeAdded(this);
  }

  return ret;
}


void SpaceNode::restoreDelaunay() {
  std::vector<Tetrahedron*>active_tetrahedra;
  std::vector<Tetrahedron*> tmp;
  for (auto tet : adjacent_tetrahedra_) {
    tmp.push_back(tet);
  }
  int upper_bound = tmp.size();
  auto sptr =  static_cast<SpaceNode*>(this);
#pragma omp simd
  for (auto i = 0; i < upper_bound; i++) {
    auto tetrahedron = tmp[i];
#pragma forceinline
    static_cast<Tetrahedron*>(tetrahedron)->updateCirumSphereAfterNodeMovement(sptr);
    active_tetrahedra.push_back(tetrahedron);
  }
  while (!active_tetrahedra.empty()) {
    int checking_index = createNewCheckingIndex();
    std::vector<Tetrahedron*>problem_tetrahedra;
    std::vector<Tetrahedron*>flat_tetrahedra;
    while (!active_tetrahedra.empty()) {
      auto tetrahedron = active_tetrahedra.front();
//      active_tetrahedra.pop_front();
      STLUtil::vectorRemove(active_tetrahedra, active_tetrahedra[0]);
      if (tetrahedron->isValid()) {
        tetrahedron->isInfinite();
        auto tet_adjacent_triangles = tetrahedron->getAdjacentTriangles();
        auto tet_adjacent_nodes = tetrahedron->getAdjacentNodes();
        for (auto i = (tetrahedron->isInfinite() ? 1 : 0); i < 4; i++) {
          auto triangle_i = tet_adjacent_triangles[i];
          // Check whether or not we already tested this
          // combination:
          if (!triangle_i->wasCheckedAlready(checking_index)) {
            auto tetrahedron_i = triangle_i->getOppositeTetrahedron(tetrahedron);
            auto node_i = tetrahedron_i->getOppositeNode(triangle_i);
            // is there a violation of the Delaunay criterion?
            if ((node_i != nullptr)
                && ((tetrahedron->isTrulyInsideSphere(node_i->getPosition())
                        || (tetrahedron->isFlat() && tetrahedron_i->isFlat())))) {
              std::vector<Tetrahedron*>new_tetrahedra;
              // check if there is a neighboring tetrahedron
              // also violating the Delaunay criterion
              for (auto j = (tetrahedron->isInfinite() ? 1 : 0); j < 4; j++) {
                if (i != j) {
                  auto triangle_j = tetrahedron->getAdjacentTriangles()[j];
                  auto tetrahedron_j = triangle_j->getOppositeTetrahedron(tetrahedron);
                  // is there also a violation of the Delaunay criterion
                  // between tetrahedron_i and tetrahedron_j?
                  if (tetrahedron_j->isNeighbor(tetrahedron_i)) {
                    auto opp_j = tet_adjacent_nodes[j];
                    auto opp_i = tet_adjacent_nodes[i];
                    if (opp_i != nullptr && opp_j != nullptr) {
                      // Either all 3 tetrahedra are flat & they are all neighbors
                      // of each other, or they are not flat but
                      // their spheres include the other tetrahedra's points
                      if ((tetrahedron->isFlat()
                              && tetrahedron_i->isFlat()
                              && tetrahedron_j->isFlat() && tetrahedron_i != tetrahedron_j)
                          || (tetrahedron_j->isTrulyInsideSphere(opp_j->getPosition())
                              && tetrahedron_i->isTrulyInsideSphere(opp_i->getPosition()))) {
                        auto array = Tetrahedron::flip3to2(tetrahedron, tetrahedron_i, tetrahedron_j);
                        STLUtil::arrayToList(array, new_tetrahedra);
                        break;
                      }
                    }
                  }
                }
              }  // for j
              // if no 3->2 flip was found, perform a 2->3  flip, if possible (convex!)
              if (new_tetrahedra.empty()) {
                if (tetrahedron->isFlat() && tetrahedron_i->isFlat() && tetrahedron->isAdjacentTo(node_i)) {
                  new_tetrahedra = Tetrahedron::remove2FlatTetrahedra(tetrahedron, tetrahedron_i);
                } else if (!(tetrahedron->isFlat() || tetrahedron_i->isFlat())) {
                  auto array = Tetrahedron::flip2to3(tetrahedron, tetrahedron_i);
                  STLUtil::arrayToList(array, new_tetrahedra);
                  // TODO resolve hack to return null array
                  bool all_null = true;
                  for (auto t : new_tetrahedra) {
                    if (t != nullptr) {
                      all_null = false;
                    }
                  }
                  if (all_null) {
                    new_tetrahedra.clear();
                  }
                  // TODO end
                }
              }
              if (!new_tetrahedra.empty()) {
                for (auto tet : new_tetrahedra) {
                  active_tetrahedra.push_back(tet);
                  if (tet->isFlat()) {
                    flat_tetrahedra.push_back(tet);
                  }
                }
                break;
              } else {
                problem_tetrahedra.push_back(tetrahedron);
                problem_tetrahedra.push_back(tetrahedron_i);
                active_tetrahedra.push_back(tetrahedron_i);
              }
            }  // if nodeI inside sphere of tetrahedron
          }  // if was checked before
        }  // for i
        // if (messedUp && tetrahedron.isValid())
        // messedUpTetrahedra.offer(tetrahedron);
      }  // if tetrahedron is valid
    }  // while

    // special case: in some situation (like an octahedron), some false
    // tetrahedra might not have been removed. We solve this problem by
    // simply removing all invalid tetrahedra and triangulating the
    // holes
    std::vector<Tetrahedron*>messed_up_tetrahedra;
    // check if there are flat tetrahedra left:
    for (auto flat_tetrahedron : flat_tetrahedra) {
      if (flat_tetrahedron->isValid()
          && !STLUtil::vectorContains(messed_up_tetrahedra, flat_tetrahedron)) {
        for (auto triangle : flat_tetrahedron->getAdjacentTriangles()) {
          auto opposite_tetrahedron = triangle->getOppositeTetrahedron(
              flat_tetrahedron);
          if (opposite_tetrahedron->isValid()
              && !STLUtil::vectorContains(messed_up_tetrahedra, opposite_tetrahedron)) {
            messed_up_tetrahedra.push_back(opposite_tetrahedron);
          }
        }
        messed_up_tetrahedra.push_back(flat_tetrahedron);
      }
    }
    // filter for messedUpTetrahedra, that are still valid and still
    // messed up:
    for (auto tetrahedron : problem_tetrahedra) {
      if (tetrahedron->isValid() && !tetrahedron->isFlat()
          && !STLUtil::vectorContains(messed_up_tetrahedra, tetrahedron)) {
        for (auto adjacent_triangle : tetrahedron->getAdjacentTriangles()) {
          auto opposite_tetrahedron = adjacent_triangle->getOppositeTetrahedron(tetrahedron);
          if (!opposite_tetrahedron->isInfinite()) {
            auto opposite_node = opposite_tetrahedron->getOppositeNode(adjacent_triangle);
            if (tetrahedron->isTrulyInsideSphere(opposite_node->getPosition())) {
              messed_up_tetrahedra.push_back(tetrahedron);
              break;
            }
          }
        }
      }
    }
    if (!messed_up_tetrahedra.empty()) {
      cleanUp(messed_up_tetrahedra);
    }
  }
}


std::array<double, 3> SpaceNode::proposeNewPosition() {
  auto min_distance = std::numeric_limits<double>::max();
  std::array<double, 3> farthest_away_diff;
  auto max_distance = std::numeric_limits<double>::min();
  for (auto edge : adjacent_edges_) {
    auto other_node = edge->getOpposite(this);
    if (other_node != nullptr) {
      auto diff = Matrix::subtract(other_node->getPosition(), position_);
      double distance = Matrix::dot(diff, diff);
      if (distance < min_distance) {
        min_distance = distance;
      }
      if (distance > max_distance) {
        max_distance = distance;
        farthest_away_diff = diff;
      }
    } else if (max_distance < std::numeric_limits<double>::max()) {
      max_distance = std::numeric_limits<double>::max();
      auto some_adjacent_tetrahedron = edge->getAdjacentTetrahedra().front();
      auto triangle = some_adjacent_tetrahedron->getAdjacentTriangles()[0];
      triangle->updatePlaneEquationIfNecessary();
      auto opposite_tetrahedron = triangle->getOppositeTetrahedron(
          some_adjacent_tetrahedron);
      farthest_away_diff = triangle->getNormalVector();
      if (!opposite_tetrahedron->isInfinite()) {
        auto outer_position = Matrix::add(position_, farthest_away_diff);
        auto position = opposite_tetrahedron->getOppositeNode(triangle)
            ->getPosition();
        if (triangle->onSameSide(outer_position, position)) {
          farthest_away_diff = Matrix::scalarMult(-1, farthest_away_diff);
        }
      }
    }
  }
  return Matrix::add(
      position_,
      Matrix::scalarMult(std::sqrt(min_distance) * 0.5,
                         Matrix::normalize(farthest_away_diff)));
}


std::vector<Edge* > SpaceNode::getAdjacentEdges() const {
  return adjacent_edges_;
}


std::string SpaceNode::toString() const {
  std::ostringstream str_stream;
  str_stream << "(";
  str_stream << StringUtil::toStr(id_);
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(volume_);
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(content_);
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(position_);
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(adjacent_tetrahedra_);
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(adjacent_edges_);
//  str_stream << ", ";
//  str_stream << StringUtil::toStr(listeners_);
  str_stream << ")";
  return str_stream.str();
}


bool SpaceNode::equalTo(const SpaceNode* other) const {
  return other == this;
}


int SpaceNode::createNewCheckingIndex() {
  checking_index_ = (checking_index_ + 1) % 2000000000;
  return checking_index_;
}


Tetrahedron* SpaceNode::removeAndReturnCreatedTetrahedron() {
  for (auto listener : listeners_) {
    listener->nodeAboutToBeRemoved(this);
  }
  auto oto = java_->oto_createSimpleOpenTriangleOrganizer();
  std::vector<Tetrahedron*>messed_up_tetrahedra;
  // Collect the triangles that are opened by removing the point and
  // remove the corresponding tetrahedrons:
  std::vector<Tetrahedron*>adjacent_tetrahedra_cpy;
  for (auto tetrahedron : adjacent_tetrahedra_) {
    adjacent_tetrahedra_cpy.push_back(tetrahedron);
  }
  for (auto tetrahedron : adjacent_tetrahedra_cpy) {
    if (tetrahedron->isValid()) {
      auto opposite_triangle = tetrahedron->getOppositeTriangle(
          this);
      oto->putTriangle(opposite_triangle);
      auto opposite_tetrahedron = opposite_triangle->getOppositeTetrahedron(
          tetrahedron);
      tetrahedron->remove();
      if ((opposite_tetrahedron != nullptr)
          && !opposite_tetrahedron->isInfinite()
          && (opposite_tetrahedron->isInsideSphere(getPosition()))) {
        messed_up_tetrahedra.push_back(opposite_tetrahedron);
      }
    }
  }
  for (auto tetrahedron : messed_up_tetrahedra) {
    if (tetrahedron->isValid()) {
      oto->removeAllTetrahedraInSphere(tetrahedron);
    }
  }
  oto->triangulate();
  for (auto listener : listeners_) {
    listener->nodeRemoved(this);
  }
  return oto->getANewTetrahedron();
}


void SpaceNode::processTetrahedron(
    Tetrahedron* tetrahedron,
    std::vector<Triangle3D* >& queue,
    OpenTriangleOrganizer* oto) {
  tetrahedron->remove();
  for (auto current_triangle : tetrahedron->getAdjacentTriangles()) {
    if (current_triangle->isCompletelyOpen()) {
      oto->removeTriangle(current_triangle);
    } else {
      queue.push_back(current_triangle);
      oto->putTriangle(current_triangle);
    }
  }
}


bool SpaceNode::checkIfTriangulationIsStillValid(
    const std::array<double, 3>& new_position) {
  for (auto tetrahedron : adjacent_tetrahedra_) {
    if (tetrahedron->isFlat()) {
      return false;
    }
    if (tetrahedron->isInfinite()) {
      auto inner_tet = tetrahedron->getAdjacentTetrahedron(0);
      if (inner_tet->getAdjacentTetrahedron(0)->isInfinite()
          && inner_tet->getAdjacentTetrahedron(1)->isInfinite()
          && inner_tet->getAdjacentTetrahedron(2)->isInfinite()
          && inner_tet->getAdjacentTetrahedron(3)->isInfinite()) {
        return true;
      } else {
        return false;
      }
    } else {
      // check for each adjacent triangle if this node would remain on the same side:
      auto triangle = tetrahedron->getOppositeTriangle(
          this);
      triangle->updatePlaneEquationIfNecessary();
      if (!triangle->trulyOnSameSide(position_, new_position)) {
        tetrahedron->testPosition(new_position);
        return false;
      }
    }
  }
  return true;
}


bool SpaceNode::removeTetrahedronDuringCleanUp(
    Tetrahedron* tetrahedron_to_remove,
    std::vector<Tetrahedron* >& list,
    std::vector<SpaceNode* >& node_list,
    OpenTriangleOrganizer* oto) {
  bool ret = false;
  for (auto node : tetrahedron_to_remove->getAdjacentNodes()) {
    if ((node != nullptr) && (!STLUtil::vectorContains(node_list, node))) {
      ret = true;
      node_list.push_back(node);
    }
  }
  for (auto adjacent_triangle : tetrahedron_to_remove->getAdjacentTriangles()) {
    auto opposite = adjacent_triangle->getOppositeTetrahedron(
        tetrahedron_to_remove);
    if ((opposite != nullptr) && (!STLUtil::vectorContains(list, opposite))) {
      list.push_back(opposite);
      ret = true;
    }
  }

  tetrahedron_to_remove->remove();
  for (auto current_triangle : tetrahedron_to_remove->getAdjacentTriangles()) {
    if (current_triangle->isCompletelyOpen()) {
      oto->removeTriangle(current_triangle);
    } else {
      oto->putTriangle(current_triangle);
    }
  }
  return ret;
}


void SpaceNode::cleanUp(
    const std::vector<Tetrahedron* >& messed_up_tetrahedra) {
  std::vector<Tetrahedron*>outer_tetrahedra;
  std::vector<SpaceNode*> problem_nodes;
  auto oto = java_->oto_createSimpleOpenTriangleOrganizer();
  for (auto tetrahedron : messed_up_tetrahedra) {
    if (tetrahedron->isValid()) {
      removeTetrahedronDuringCleanUp(tetrahedron, outer_tetrahedra, problem_nodes, oto);
      if (STLUtil::vectorContains(outer_tetrahedra, tetrahedron)) {
        STLUtil::vectorRemove(outer_tetrahedra, tetrahedron);
      }
    }
  }
  bool done = false;
  while (!done) {
    Tetrahedron* problem_tetrahedron(nullptr);
    for (auto outer_tetrahedron : outer_tetrahedra) {
      if (outer_tetrahedron->isValid())
      for (auto node : problem_nodes) {
        if (!outer_tetrahedron->isAdjacentTo(node)) {
          if (outer_tetrahedron->isFlat() || outer_tetrahedron->isInsideSphere(node->getPosition())) {
            removeTetrahedronDuringCleanUp( outer_tetrahedron, outer_tetrahedra, problem_nodes, oto);
            problem_tetrahedron = outer_tetrahedron;
            break;
          }
        }
      }
      if (problem_tetrahedron != nullptr) {
        break;
      }
    }
    if (problem_tetrahedron != nullptr) {
      STLUtil::vectorRemove(outer_tetrahedra, problem_tetrahedron);
    } else {
      done = true;
    }
  }
  oto->triangulate();
}


}  // namespace spatial_organization
}  // namespace cx3d
