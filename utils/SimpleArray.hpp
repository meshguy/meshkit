#ifndef SIMP_ARR_H
#define SIMP_ARR_H

#include <algorithm>
#include <cstdlib>

template <typename T> 
class SimpleArray
{
private:
  T* arr;
  int arrSize;
  int arrAllocated;
     
public:
  SimpleArray() : arr(NULL) , arrSize(0), arrAllocated(0) {}

  SimpleArray( unsigned s ) : arr(NULL), arrSize(s), arrAllocated(s) 
  {
    resize(s);
  }
    
  ~SimpleArray() {
    clear();
  }

  T**  ptr()            { return &arr; }
  int& size()           { return arrSize; }
  int  size()     const { return arrSize; }
  int& capacity()       { return arrAllocated; }
  int  capacity() const { return arrAllocated; }

  void resize( unsigned s )
  {
    clear();
    if (s == 0) return;

    arr = static_cast<T*>(malloc(s*sizeof(T)));
    for (unsigned i = 0; i < s; ++i)
      new (arr+i) T();
    arrSize = arrAllocated = s;
  }

  void clear()
  {
    if (arr == NULL) return;
    for (int i = 0; i < size(); ++i)
      arr[i].~T();
    free(arr); 
    arr = NULL;
    arrSize = arrAllocated = 0;
  }

  void swap(SimpleArray &other)
  {
    using std::swap;
    swap(arr, other.arr);
    swap(arrSize, other.arrSize);
    swap(arrAllocated, other.arrAllocated);
  }

  typedef T* iterator;
  typedef const T* const_iterator;
  iterator       begin()       { return arr; }
  const_iterator begin() const { return arr; }
  iterator         end()       { return arr + arrSize; }
  const_iterator   end() const { return arr + arrSize; }

  T& operator[]( unsigned idx )       { return arr[idx]; }
  T  operator[]( unsigned idx ) const { return arr[idx]; }
};

#define ARRAY_INOUT( A ) A.ptr(), &A.capacity(), &A.size()
#define ARRAY_IN( A ) &A[0], A.size()

namespace std
{
  // Alas, this is non-standard, but it works everywhere, and there's no other
  // way to do it.
  template<typename T>
  void swap(SimpleArray<T> &lhs, SimpleArray<T> &rhs)
  {
    lhs.swap(rhs);
  }
}

#endif
