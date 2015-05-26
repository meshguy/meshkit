//////////////////////////////////////////
// myany.hpp
#include <memory>
#include <stdexcept>
#include <utility>

#ifdef __APPLE__
#  define unique_ptr auto_ptr
#endif

struct myany
{
    myany() = default;
    template <typename T> myany(T const& v) : _storage(new storage<T>(v)) { }
    myany(myany const& other) : _storage()
    {
        if (other._storage)
          #ifdef __APPLE__
                  if (other._storage.get())
                      _storage = other._storage->clone();
          #else // __APPLE__
                  if (other._storage.get())
                      _storage = std::move(other._storage->clone());
          #endif // __APPLE__
    }

    //void swap(myany& other)               { _storage.swap(other._storage); }
        void swap(myany& other)
          {
    #ifdef __APPLE__
          std::unique_ptr<storage_base> tmp = _storage;
          _storage = other._storage;
          other._storage = tmp;
    #else // __APPLE__
          _storage.swap(other._storage);
    #endif // __APPLE__
          }
    friend void swap(myany& a, myany& b) { a.swap(b); }
    myany& operator=(myany other)        { swap(other); return *this; }

    // todo move semantics
private:
    struct storage_base {
        virtual std::unique_ptr<storage_base> clone() = 0;
        virtual ~storage_base() { }
    };
    template <typename T>
    struct storage : storage_base {
        T value;
        explicit storage(T const& v) : value(v) {}
        std::unique_ptr<storage_base> clone() { return std::unique_ptr<storage_base>(new storage<T>(value)); }
    };
    std::unique_ptr<storage_base> _storage;
    template<typename T> friend T      & any_cast(myany      &);
    template<typename T> friend T const& any_cast(myany const&);
};

template <typename T> T& any_cast(myany& a) {
    if (auto p = dynamic_cast<myany::storage<T>*>(a._storage.get()))
        return p->value;
    else
        throw std::bad_cast();
}

template <typename T> T const& any_cast(myany const& a) {
    if (auto p = dynamic_cast<myany::storage<T> const*>(a._storage.get()))
        return p->value;
    else
        throw std::bad_cast();
}
