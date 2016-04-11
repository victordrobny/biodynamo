#ifndef SPATIAL_ORGANIZATION_BINARY_TREE_ELEMENT_H_
#define SPATIAL_ORGANIZATION_BINARY_TREE_ELEMENT_H_

#include <vector> //former list
#include <memory>
#include <string>

namespace cx3d {
namespace spatial_organization {

 class SpaceNode;
 class SimpleTriangulationNodeOrganizer;


class BinaryTreeElement {
  friend class SimpleTriangulationNodeOrganizer;
 public:
  static BinaryTreeElement* generateTreeHead();

  BinaryTreeElement(SpaceNode* content);
  virtual ~BinaryTreeElement();

  virtual bool contains(const SpaceNode* content) const;

  virtual void insert(SpaceNode* content);

  virtual void remove(const SpaceNode* content,
                      BinaryTreeElement* parent);

  // todo replace with STL iterator
  virtual std::vector<SpaceNode*>inOrderTraversal() const;

  std::string toString() const;

 protected:
  SpaceNode* content_;
  BinaryTreeElement* smaller_;
  BinaryTreeElement* bigger_;
  int content_id_ = 0;

 private:
  BinaryTreeElement() = delete;
  BinaryTreeElement(const BinaryTreeElement&) = delete;
  BinaryTreeElement& operator=(const BinaryTreeElement&) = delete;

  int getHash(const SpaceNode* content) const;

  bool contains(int id, const SpaceNode* content) const;

  void insert(BinaryTreeElement* element);

  void remove(int id, const SpaceNode* content, BinaryTreeElement* parent);

  void changeLink(BinaryTreeElement* old_el, BinaryTreeElement* new_el);
};


class TreeHead : public BinaryTreeElement {
 public:
  TreeHead();

  bool contains(const SpaceNode* content) const override;

  void insert(SpaceNode* content) override;

  void remove(const SpaceNode* content,
              BinaryTreeElement* parent) override;

  std::vector<SpaceNode*>inOrderTraversal() const override;

private:
  TreeHead(const TreeHead&) = delete;
  TreeHead& operator=(const TreeHead&) = delete;
};

} // namespace spatial_organization
} // namespace cx3d

#endif // SPATIAL_ORGANIZATION_BINARY_TREE_ELEMENT_H_
