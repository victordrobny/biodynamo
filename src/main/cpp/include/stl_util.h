#ifndef STL_UTIL_H_
#define STL_UTIL_H_

#include <vector> //former list
#include <array>
#include <list>
#include <deque>
#include <unordered_map>
#include <algorithm>

class STLUtil {
 public:
  /**
   * returns whether an element is contained in a std::vector
   * Caution: linear runtime
   */
  template<typename C>
  static bool listContains(const std::list<C>& l, C element) {
    return std::find(l.begin(), l.end(), element) != l.end();
  }

  template<typename C>
  static bool vectorContains(const std::vector<C>& l, C element) {
    return std::find(l.begin(), l.end(), element) != l.end();
  }

  /**
   * copies the elements of the std::array into std::vector
   */
  template<typename C, std::size_t N>
  static void arrayToList(const std::array<C, N>& arr, std::vector<C>& list) {
    list.clear();
    for (auto el : arr) {
      list.push_back(el);
    }
  }

  /**
   * returns whether an element with Key key is stored in this map
   * @param map
   * @param key
   */
  template<typename K, typename V>
  static bool mapContains(const std::unordered_map<K, V>& map, const K& key) {
    return map.find(key) != map.end();
  }

  /**
   * returns -1 if val < 0
   *          0 if val == 0
   *          1 if val > 0
   */
  template<typename T>
  static int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

  template<typename T>
  static void vectorRemove(std::vector<T>& vector, const T& element) {
    auto it = std::find(vector.begin(), vector.end(), element);
    vector.erase(it);
  }

  template<typename T>
  static void dequeRemove(std::deque<T>& vector, const T& element) {
    auto it = std::find(vector.begin(), vector.end(), element);
    vector.erase(it);
  }

  template<typename T>
  static std::vector<T> dequeToVector(const std::deque<T>& deque) {
    std::vector<T> ret;
    for (auto el : deque) {
      ret.push_back(ret);
    }
    return ret;
  }

  template<typename T>
  static std::deque<T> vectorToDeque(const std::vector<T>& vector) {
    std::deque<T> ret;
    for (auto el : vector) {
      ret.push_back(el);
    }
    return ret;
  }

 private:
  STLUtil() = delete;
};

#endif // STL_UTIL_H_
