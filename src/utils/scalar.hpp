#ifndef UTILS_SCALAR
#define UTILS_SCALAR

namespace Utils {

/** calculates the scalar product of two vectors a nd b */
template <typename T1, typename T2> double scalar(const T1 &a, const T2 &b) {
  double d2 = 0.0;
  for (int i = 0; i < 3; i++)
    d2 += a[i] * b[i];
  return d2;
}


}

#endif
