#ifndef JAVA_UTIL_H_
#define JAVA_UTIL_H_

#include <array>
#include <memory>

#include <spatial_organization/open_triangle_organizer.h>

namespace cx3d {

/**
 * Contains functions to access former static methods that are still implemented in Java
 */
template<class T>
class JavaUtil {
 public:
  virtual ~JavaUtil() {
  }

  /**
   * returns a
   */
  virtual std::array<int, 4> generateTriangleOrder() {
    return {0, 1, 2, 3};
  }

  /**
   * redirects call, because static methods cannot be handled by SWIG direcotr
   */
  virtual std::shared_ptr<spatial_organization::OpenTriangleOrganizer<T>> oto_createSimpleOpenTriangleOrganizer() {
    return spatial_organization::OpenTriangleOrganizer<T>::createSimpleOpenTriangleOrganizer();
  }
};

/**
 *  java provided functions without templating
 */
class JavaUtil2 {
 public:
  virtual ~JavaUtil2(){
  }

  virtual std::array<double, 3> matrixRandomNoise3(double k) {
//    //fnoexceptionthrow std::logic_error(
//        "JavaUtil::matrixRandomNoise must never be called - Java must provide implementation at this point");
  }
};

}  // namespace cx3d

#endif  // JAVA_UTIL_H_
