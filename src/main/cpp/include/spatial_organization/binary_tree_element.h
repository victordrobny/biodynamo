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

  BinaryTreeElement(std::shared_ptr<SpaceNode> content);
  virtual ~BinaryTreeElement();

  virtual bool contains(const std::shared_ptr<SpaceNode>& content) const;

  virtual void insert(const std::shared_ptr<SpaceNode>& content);

  virtual void remove(const std::shared_ptr<SpaceNode>& content,
                      BinaryTreeElement* parent);

  // todo replace with STL iterator
  virtual std::vector<std::shared_ptr<SpaceNode>>inOrderTraversal() const;

  std::string toString() const;

 protected:
  std::shared_ptr<SpaceNode> content_;
  BinaryTreeElement* smaller_;
  BinaryTreeElement* bigger_;
  int content_id_ = 0;

 private:
  BinaryTreeElement() = delete;
  BinaryTreeElement(const BinaryTreeElement&) = delete;
  BinaryTreeElement& operator=(const BinaryTreeElement&) = delete;

  int getHash(std::shared_ptr<SpaceNode> content) const;

  bool contains(int id, const std::shared_ptr<SpaceNode>& content) const;

  void insert(BinaryTreeElement* element);

  void remove(int id, const std::shared_ptr<SpaceNode>& content, BinaryTreeElement* parent);

  void changeLink(BinaryTreeElement* old_el, BinaryTreeElement* new_el);
};


class TreeHead : public BinaryTreeElement {
 public:
  TreeHead();

  bool contains(const std::shared_ptr<SpaceNode>& content) const override;

  void insert(const std::shared_ptr<SpaceNode>& content) override;

  void remove(const std::shared_ptr<SpaceNode>& content,
              BinaryTreeElement* parent) override;

  std::vector<std::shared_ptr<SpaceNode>>inOrderTraversal() const override;

private:
  TreeHead(const TreeHead&) = delete;
  TreeHead& operator=(const TreeHead&) = delete;
};

} // namespace spatial_organization
} // namespace cx3d

#endif // SPATIAL_ORGANIZATION_BINARY_TREE_ELEMENT_H_
