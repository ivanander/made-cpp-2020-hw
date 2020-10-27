#ifndef CHUCK_ALLOCATOR__CHUCK_ALLOCATOR_H_
#define CHUCK_ALLOCATOR__CHUCK_ALLOCATOR_H_
#include <cstdint>
#include <utility>
#include <limits>
#include <exception>

constexpr std::size_t chuck_size = 10'000'000;

template <typename T>
struct chuck {
  chuck() = delete;
  chuck(chuck* prev) : data_(new char[chuck_size]), prev_(prev), used_(0) {}
  ~chuck() { delete [] data_; }
  char* data_;
  chuck* prev_;
  std::size_t used_;
};

template <typename T>
struct chuck_allocator {
 public:
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;
  using difference_type	= std::ptrdiff_t;
  using value_type = T;
  template<class U> struct rebind { typedef chuck_allocator<U> other; };

  chuck_allocator()
      : chuck_list_(new chuck<T>(nullptr)),
        counter_ptr_(new std::size_t(1)) {}

  chuck_allocator(const chuck_allocator& other) {
    chuck_list_ = other.chuck_list_;
    counter_ptr_ = other.counter_ptr_;
    ++(*counter_ptr_);
  }

  void swap(chuck_allocator& other) {
    std::swap(chuck_list_, other.chuck_list_);
    std::swap(counter_ptr_, other.counter_ptr_);
  }

  ~chuck_allocator() {
    --(*counter_ptr_);

    if (!(*counter_ptr_)) {
      while (chuck_list_) {
        chuck<T>* tmp = chuck_list_;
        chuck_list_ = chuck_list_->prev_;
        delete tmp;
      }
      delete chuck_list_;
      delete counter_ptr_;
    }
  }

  chuck_allocator& operator=(const chuck_allocator& other) {
    if (this == &other)
      return *this;

    chuck_allocator<T> tmp(other);
    swap(tmp);
    return *this;
  }

  pointer allocate(std::size_t n) {
    if (std::numeric_limits<std::size_t>::max() / sizeof(T) < n)
      throw std::bad_array_new_length();

    std::size_t need_to_allocate = n * sizeof(T);
    for (auto ptr = chuck_list_; ptr; ptr = ptr->prev_) {
      if (chuck_size > ptr->used_ + need_to_allocate) {
        auto res_ptr = ptr->data_ + ptr->used_;
        ptr->used_ += need_to_allocate;
        return (T*) res_ptr;
      }
    }

    chuck_list_ = new chuck<T>(chuck_list_);
    chuck_list_->used_ += need_to_allocate;
    return (T*) (chuck_list_->data_);
  }

  void deallocate(T* p, std::size_t n) {
    //delete [] (char*)p;
  }


  template <class... Args>
  void construct(T* p, Args&&... args) {
    new (p) T(std::forward<Args>(args)...);
  }

  void destroy(T* p) {
    p->~T();
  }

  std::size_t get_counter() const {
    if (counter_ptr_) {
      return *counter_ptr_;
    }
  }

 private:
  chuck<T>* chuck_list_;
  std::size_t* counter_ptr_;
};

#endif //CHUCK_ALLOCATOR__CHUCK_ALLOCATOR_H_
