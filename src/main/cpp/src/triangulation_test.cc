#include <iostream>
#include <array>
#include <memory>
#include <vector>
#include <cmath>
#include <chrono>
#include <ctime>

#include "java_util.h"
#include "physics/physical_node.h"
#include "spatial_organization/space_node.h"
#include "spatial_organization/edge.h"

using namespace cx3d::spatial_organization;
using cx3d::physics::PhysicalNode;
using cx3d::JavaUtil;

double random(double last) {
  return std::fmod(last * 313, 1009.0);
}

std::array<double, 3> nextPos(const std::array<double, 3>& last) {
  return {random(last[0]), random(last[1]), random(last[2])};
}

int main() {
  SpaceNode::setJavaUtil(std::shared_ptr<JavaUtil>{new JavaUtil()});
  std::vector<std::shared_ptr<SpaceNode > > space_nodes;
  std::array<double, 3> pos { 463.7047970232077, 439.8653887819098, 447.1949176631939 };
  auto initial_sn = SpaceNode::create(pos, std::shared_ptr<PhysicalNode> { nullptr });
  auto last_pos = pos;
  auto start_ts = std::chrono::high_resolution_clock::now();
  for (auto i = 0; i < 2000; i++) {
    pos = nextPos(last_pos);
    space_nodes.push_back(initial_sn->getNewInstance(pos, std::shared_ptr<PhysicalNode> { nullptr }));
    last_pos = pos;
  }

  auto create_ts = std::chrono::high_resolution_clock::now();
  std::array<double, 3> direction {0.08741642977919392,-0.020565131563058878, -0.03049639415937795};
  for (auto i = 0; i < 100; i++) {
    for (auto sn : space_nodes) {
      sn->moveFrom(direction);
    }
  }
  auto move_ts = std::chrono::high_resolution_clock::now();
  std::cout << "dt create " << std::chrono::duration_cast<std::chrono::milliseconds>(create_ts - start_ts).count() << std::endl;
  std::cout << "dt move " << std::chrono::duration_cast<std::chrono::milliseconds>(move_ts - create_ts).count() << std::endl;

  // validate
  int hash = 1;
  for (auto sn : space_nodes) {
    for (auto se : sn->getEdges()) {
      hash = hash * 17 + se->getOpposite(sn)->getId();
    }
  }

  // assert
  if (hash == -1669558583) {
    std::cout << "Test successful" << std::endl;
    return 0;
  } else {
    std::cout << "wrong test result:" << hash << std::endl;
    return -1;
  }
}



