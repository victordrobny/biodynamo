#ifndef SPATIAL_ORGANIZATION_EDGE_H_
#define SPATIAL_ORGANIZATION_EDGE_H_

#include <memory>
#include <string>
#include <vector> //former list

#include "spatial_organization/spatial_organization_edge.h"

#ifdef EDGE_DEBUG
#include "spatial_organization/debug/edge_debug.h"
#endif


namespace cx3d {
namespace spatial_organization {

 class SpaceNode;
 class Tetrahedron;

/**
 * This class is used to represent an edge in a triangulation. Each edge has two endpoint and
 * additionally stores links to all tetrahedra adjacent to it.
 *
 * @param  The type of the user objects stored in the endpoints of an edge.
 */

class Edge : public SpatialOrganizationEdge {
 public:
#ifndef EDGE_NATIVE
  Edge()
      : a_(),
        b_(),
        cross_section_area_(0.0),
        adjacent_tetrahedra_() {
  }
#endif
  /**
   * Creates a new Edge object and returns it within a <code>std::shared_ptr</code>
   * @see Edge(...)
   *
   * If functions return a std::shared_ptr of <code>*this</code> using
   * <code>return this;</code>, the following precondition must be met:
   * There must be at least one std::shared_ptr p that owns *this!
   * Calling <code>shared_from_this</code> on a non-shared object results in undefined behaviour.
   * http://mortoray.com/2013/08/02/safely-using-enable_shared_from_this/
   *
   * Therefore, all constructors are private and accessed through static factory methods that return
   * a std::shared_ptr.
   *
   * TODO(lukas) SWIG doesn't seem to support variadic templates and perfect forwarding system.
   * Once mapping to Java is not needed anymore, replace following create functions with:
   * <code>
   * template<typename ... T>
   * static Edge* create(T&& ... all) {
   *   return Edge*(new Edge(std::forward(all)...));
   * }
   * </code>
   */
  static Edge* create(SpaceNode* a,
                                         SpaceNode* b) {
#ifdef EDGE_DEBUG
    Edge* edge(new EdgeDebug(a, b));
#else
    Edge* edge = new Edge(a, b);
#endif
    edge->initializationHelper();
    return edge;
  }


  virtual ~Edge() {
  }

  /**
   * {@inheritDoc}
   */
  SpaceNode* getOpposite(const SpaceNode* node) const
      override;

  /**
   * {@inheritDoc}
   */
  physics::PhysicalNode* getOppositeElement(const physics::PhysicalNode* first) const override;

  /**
   * {@inheritDoc}
   */
  physics::PhysicalNode* getFirstElement() const override;

  /**
   * {@inheritDoc}
   */
  physics::PhysicalNode* getSecondElement() const override;

  /**
   * {@inheritDoc}
   */
  double getCrossSection() const override;

  /**
   *  @return A string representation of this edge
   */
  virtual const std::string toString() const;

  /**
   * Determines if two instances of this object are equal
   */
  virtual bool equalTo(const Edge* other);

  /**
   * Tests whether this edge is connecting a pair of points.
   * @param a The first node.
   * @param b The second node.
   * @return <code>true</code>, if this edge connects <code>a</code> and <code>b</code>.
   */
  virtual bool equals(const SpaceNode* a,
                      const SpaceNode* b) const;
  /**
   * Removes a tetrahedron from this edge's list of tetrahedra. If this edge is not incident to
   * any tetrahedra after the removal of the specified tetrahedron, the edge removes itself from
   * the triangulation by calling {@link #remove()}.
   * @param tetrahedron A tetrahedron incident to this edge which should be removed.
   */
  virtual void removeTetrahedron(Tetrahedron* tetrahedron);

  /**
   * Adds a tetrahedron to this edge's list of tetrahedra.
   * @param tetrahedron A tetrahedron incident to this edge which should be added.
   */
  virtual void addTetrahedron(Tetrahedron* tetrahedron);

  /**
   * Removes this edge from the triangulation. To do so, the two endpoints are informed
   * that the edge was removed.
   */
  virtual void remove();

  /**
   * Returns the list of incident tetrahedra.
   * @return The list of incident tetrahedra.
   */
  virtual std::vector<Tetrahedron* > getAdjacentTetrahedra() const;

  /**
   * Changes the cross section area of this edge.
   * @param change The value by which the cross section area of this edge has changed.
   * At initialization, this area is set to zero and all tetrahedra that are registered as
   * incident tetrahedra increase the cross section area.
   */
  virtual void changeCrossSectionArea(double change);

 protected:
  /**
   * Initializes a new edge with two specified endpoints.
   * @param a The first endpoint of the new edge.
   * @param b The second endpoint of the new edge.
   */
  Edge(SpaceNode* a, SpaceNode* b);

 private:
#ifdef EDGE_NATIVE
  Edge() = delete;
#endif
  Edge(const Edge&) = delete;
  Edge& operator=(const Edge&) = delete;

  /**
   * The two endpoints of this edge.
   */
  SpaceNode* a_;
  SpaceNode* b_;

  /**
   * A list of all tetrahedra that are adjacent to this edge.
   */
  std::vector<Tetrahedron* > adjacent_tetrahedra_;

  /**
   * Stores the cross section area associated with this edge.
   */
  double cross_section_area_;

  /**
   * Initialization code that cannot be called inside the constructor, because it passes
   * a std::shared_ptr of itself to another function.
   * It does that using this. @see create function documentation for a more
   * detailed explanation of the requirements before calling this
   */
  void initializationHelper();
};

}  // namespace spatial_organization
}  // namespace cx3d

#endif  // SPATIAL_ORGANIZATION_EDGE_H_
