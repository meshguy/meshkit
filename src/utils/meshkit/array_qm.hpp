#ifdef __APPLE__
#include <tr1/array>
namespace meshkit{
  using std::tr1::array;
}
#endif

#ifdef __linux__
#include <array>
namespace meshkit{
  using std::array;
}
#endif
