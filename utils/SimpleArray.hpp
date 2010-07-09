#ifndef SIMPLEARRAY_H
#define SIMPLEARRAY_H

#include <algorithm>
#include <cstdlib>

template <typename T> 
class SimpleArray
{
private:
  SimpleArray(const SimpleArray &);
  SimpleArray & operator = (const SimpleArray &);

  T *arr;
  int arrSize;
  int arrAllocated;
public:
  SimpleArray() : arr(NULL), arrSize(0), arrAllocated(0) {}

  SimpleArray( unsigned s ) : arr(NULL), arrSize(s), arrAllocated(s) 
  {
    resize(s);
  }
    
  ~SimpleArray() {
    clear();
  }

  int size()     const { return arrSize; }
  int capacity() const { return arrAllocated; }

  void resize( unsigned s )
  {
    clear();
    if (s == 0) return;
    arr = (T*)malloc(s*sizeof(T));
    for (unsigned i = 0; i < s; ++i)
      new (arr+i) T();
    arrSize = arrAllocated = s;
  }

  void clear() 
  {
    if (arr == NULL) return;
    for (int i = 0; i < arrSize; ++i)
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

  T&       operator[]( unsigned idx )       { return arr[idx]; }
  const T& operator[]( unsigned idx ) const { return arr[idx]; }

  // Don't use these directly!
  T**  ptr_()      { return &arr; }
  int* size_()     { return &arrSize; }
  int* capacity_() { return &arrAllocated; }
};

#define ARRAY_INOUT( A ) (A).ptr_(), (A).capacity_(), (A).size_()
#define ARRAY_IN( A ) &(A)[0], (A).size()

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
