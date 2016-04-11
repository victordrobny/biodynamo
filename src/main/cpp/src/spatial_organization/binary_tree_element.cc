#include "spatial_organization/binary_tree_element.h"

#include <stdint.h>

#include <stack>
#include <sstream>

#include "string_util.h"
#include "physics/physical_node.h"
#include "spatial_organization/space_node.h"

namespace cx3d {
namespace spatial_organization {


BinaryTreeElement* BinaryTreeElement::generateTreeHead() {
  return new TreeHead();
}


BinaryTreeElement::BinaryTreeElement(SpaceNode* content)
    : content_ { content },
      bigger_ { nullptr },
      smaller_ { nullptr } {
  if (content_ != nullptr) {
    content_id_ = getHash(content);
  } else {
    content_id_ = -1;
  }
}


BinaryTreeElement::~BinaryTreeElement() {
  delete smaller_;
  delete bigger_;
}


bool BinaryTreeElement::contains(
    const SpaceNode* content) const {
  return contains(getHash(content), content);
}


void BinaryTreeElement::insert(
    SpaceNode* content) {
  if (content != nullptr) {
    insert(new BinaryTreeElement(content));
  }
}


void BinaryTreeElement::remove(const SpaceNode* content,
                                  BinaryTreeElement* parent) {
  remove(getHash(content), content, parent);
}


std::string BinaryTreeElement::toString() const {
  stringstream str_stream;
  str_stream << StringUtil::toStr(smaller_);
  str_stream << ", ";
  str_stream << StringUtil::toStr(content_);
  str_stream << ", ";
  str_stream << StringUtil::toStr(bigger_);
  return str_stream.str();
}


int BinaryTreeElement::getHash(const SpaceNode* content) const {
  uint64_t id = content_->getId();
  uint64_t c = 7481;
  return (id * c) % 74317;
}


bool BinaryTreeElement::contains(
    int id, const SpaceNode* content) const {
  return contains(getHash(content), content);
}


void BinaryTreeElement::insert(BinaryTreeElement* element) {
  if (content_id_ == element->content_id_
      && content_ == element->content_) {
    return;
  } else if ((content_id_ >= element->content_id_)) {
    if ((smaller_ != nullptr)) {
      smaller_->insert(element);
    } else {
      smaller_ = element;
    }
  } else if (content_id_ < element->content_id_) {
    if (bigger_ != nullptr) {
      bigger_->insert(element);
    } else {
      bigger_ = element;
    }
  }
}


void BinaryTreeElement::remove(int id,
                                  const SpaceNode* content,
                                  BinaryTreeElement* parent) {
  if ((content_id_ == id) && (content_ == content)) {
    if ((smaller_ == nullptr) && (bigger_ == nullptr)) {
      parent->changeLink(this, nullptr);
      //use of randomization in the next if showed no influence on the simulation outcome
    } else if ((smaller_ != nullptr) || (bigger_ == nullptr)) {
      parent->changeLink(this, smaller_);
      if (bigger_ != nullptr) {
        smaller_->insert(bigger_);
      }
    } else {
      parent->changeLink(this, bigger_);
      if (smaller_ != nullptr) {
        bigger_->insert(smaller_);
      }
    }
  } else {
    if ((content_id_ >= id) && (smaller_ != nullptr)) {
      smaller_->remove(id, content, this);

    } else if ((content_id_ < id) && (bigger_ != nullptr)) {
      bigger_->remove(id, content, this);
    }
  }
}


void BinaryTreeElement::changeLink(BinaryTreeElement* old_el,
                                      BinaryTreeElement* new_el) {
  if (smaller_ == old_el) {
    smaller_ = new_el;
  } else if (bigger_ == old_el) {
    bigger_ = new_el;
  }
}


std::vector<SpaceNode*>BinaryTreeElement::inOrderTraversal() const {
  std::vector<SpaceNode*> traversal;
  std::stack<const BinaryTreeElement*> stack;
  const BinaryTreeElement* dummy = this;
  while(dummy != nullptr) {
    stack.push(dummy);
    dummy = dummy->smaller_;
  }

  while(!stack.empty()) {
    dummy = stack.top();
    stack.pop();
    auto it = dummy->bigger_;
    while(it != nullptr) {
      stack.push(it);
      it = it->smaller_;
    }
    traversal.push_back(dummy->content_);
  }
  return traversal;
}

//-------------------------------------------------------------------------------------------------
// TreeHead


TreeHead::TreeHead()
    : BinaryTreeElement(nullptr) {

}


bool TreeHead::contains(const SpaceNode* content) const {
  return
      BinaryTreeElement::bigger_ != nullptr ?
          BinaryTreeElement::bigger_->contains(content) : false;
}


void TreeHead::insert(SpaceNode* content) {
  if (BinaryTreeElement::bigger_ != nullptr) {
    BinaryTreeElement::bigger_->insert(content);
  } else {
    BinaryTreeElement::bigger_ = new BinaryTreeElement(content);
  }
}


void TreeHead::remove(const SpaceNode* content,
                         BinaryTreeElement* parent) {
  if (BinaryTreeElement::bigger_ != nullptr) {
    BinaryTreeElement::bigger_->remove(content, this);
  }
}


std::vector<SpaceNode*>TreeHead::inOrderTraversal() const {
  if(BinaryTreeElement::bigger_ != nullptr) {
    return BinaryTreeElement::bigger_->inOrderTraversal();
  }
  return std::vector<SpaceNode*>();
}


}  // namespace spatial_organization
}  // namespace cx3d
