
// SIMPLE ARRAY TEMPLATE
/* Frees allocated arrays for us */
template <typename T> class SimpleArray
{
private:
  T* arr;
  int arrSize;
  int arrAllocated;
  int i;
public:
  SimpleArray() : arr(0) , arrSize(0), arrAllocated(0) {}
  SimpleArray( unsigned s ) :arrSize(s), arrAllocated(s) {
    arr = (T*)malloc(s*sizeof(T));
    for (unsigned i = 0; i < s; ++i)
      new (arr+i) T();
  }
  void setSize( int &s ){
    arrSize = s; 
    arrAllocated = s;
    arr = (T*)malloc(s*sizeof(T));
    for (i = 0; i < s; ++i)
      new (arr+i) T();
  }    
  ~SimpleArray() {
    for (i = 0; i < size(); ++i)
      arr[i].~T();
    free(arr);
  }

  T**  ptr()            { return &arr; }
  int& size()           { return arrSize; }
  int  size()     const { return arrSize; }
  int& capacity()       { return arrAllocated; }
  int  capacity() const { return arrAllocated; }
    
  typedef T* iterator;
  typedef const T* const_iterator;
  iterator       begin()       { return arr; }
  const_iterator begin() const { return arr; }
  iterator         end()       { return arr + arrSize; }
  const_iterator   end() const { return arr + arrSize; }
    
    
  T& operator[]( unsigned idx )       { return arr[idx]; }
  T  operator[]( unsigned idx ) const { return arr[idx]; }
};
