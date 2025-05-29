#include <iostream>
#include <memory>
#include <type_traits>

struct BaseControlBlock {
  size_t shared_count = 0;
  size_t weak_count = 0;

  virtual void destroy_obj() = 0;
  virtual void destroy_block() = 0;

  virtual ~BaseControlBlock() = default;
};

template <typename T, typename Deleter = std::default_delete<T>,
          typename Alloc = std::allocator<T>>
struct ControlBlockRegular : BaseControlBlock {

  using AllocTraits = std::allocator_traits<Alloc>;
  using BlockAlloc = typename AllocTraits::template rebind_alloc<
      ControlBlockRegular<T, Deleter, Alloc>>;
  using BlockAllocTraits = std::allocator_traits<BlockAlloc>;

  T* ptr = nullptr;
  Deleter del;
  Alloc alloc;

  ControlBlockRegular(T* ptr, Deleter del = Deleter(), Alloc alloc = Alloc())
      : ptr(ptr), del(del), alloc(alloc) {}

  void destroy_obj() override {
    del(ptr);
    ptr = nullptr;
  }

  void destroy_block() override {
    BlockAlloc temp(alloc);
    BlockAllocTraits::deallocate(temp, this, 1);
  }

  ~ControlBlockRegular() override = default;
};

template <typename T, typename Alloc = std::allocator<T>>
struct ControlBlockMakeShared : BaseControlBlock {

  using AllocTraits = std::allocator_traits<Alloc>;
  using BlockAlloc = typename AllocTraits::template rebind_alloc<
      ControlBlockMakeShared<T, Alloc>>;
  using BlockAllocTraits = std::allocator_traits<BlockAlloc>;

  T value;
  Alloc alloc = alloc();

  template <typename... Args>
  ControlBlockMakeShared(const Alloc& alloc, Args&&...args)
      : value(std::forward<Args>(args)...), alloc(alloc) {}

  void destroy_obj() override { AllocTraits::destroy(alloc, &value); }

  void destroy_block() override {
    BlockAlloc temp(alloc);
    BlockAllocTraits::deallocate(temp, this, 1);
  }

  ~ControlBlockMakeShared() override = default;
};


template<typename T>
class WeakPtr;

template <typename T> class SharedPtr {

  template <typename>
  friend class SharedPtr;

  template <typename U>
  friend class WeakPtr;

  template <typename U, typename Alloc, typename... Args>
  friend SharedPtr<U> allocateShared(const Alloc& alloc, Args&&...args);

  template <typename U, typename... Args>
  friend SharedPtr<U> makeShared(Args&&...args);

  T* ptr = nullptr;
  BaseControlBlock* pcb = nullptr;

  template <typename U, typename Alloc = std::allocator<U>>
  SharedPtr(ControlBlockMakeShared<U, Alloc>* cb) : ptr(&(cb->value)), pcb(cb) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

  template <typename U>
  SharedPtr(const WeakPtr<U>& other) : ptr(other.ptr), pcb(other.pcb) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

public:
  SharedPtr() = default;

  template <typename U>
  explicit SharedPtr(U* ptr) : ptr(ptr), pcb(new ControlBlockRegular<T>(ptr)) {
    pcb->shared_count++;
  }

  SharedPtr(const SharedPtr& other) : ptr(other.ptr), pcb(other.pcb) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

  template <typename U>
  SharedPtr(const SharedPtr<U>& other) : ptr(other.ptr), pcb(other.pcb) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

  template <typename U>
  SharedPtr(SharedPtr<U>&& other) : ptr(other.ptr), pcb(other.pcb)  {
    other.pcb = nullptr;
    other.ptr = nullptr;
  }

  SharedPtr(SharedPtr&& other) : ptr(other.ptr), pcb(other.pcb) {
    other.pcb = nullptr;
    other.ptr = nullptr;
  }

  template <typename U, typename Deleter>
  SharedPtr(U* ptr, Deleter del) 
      : ptr(ptr), pcb(new ControlBlockRegular<U, Deleter>(ptr, del)) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

  template <typename U, typename Deleter, typename Alloc>
  SharedPtr(U* ptr, Deleter del, Alloc alloc) : ptr(ptr) {

    using AllocTraits = std::allocator_traits<Alloc>;
    using BlockAlloc = typename AllocTraits::template rebind_alloc<
        ControlBlockRegular<U, Deleter, Alloc>>;
    using BlockAllocTraits = std::allocator_traits<BlockAlloc>;

    BlockAlloc balloc(alloc);
    auto* cb = BlockAllocTraits::allocate(balloc, 1);
    new (cb) ControlBlockRegular<U, Deleter, Alloc>(ptr, del, alloc);
    pcb = cb;
    if (pcb) {
      pcb->shared_count++;
    }
  }

  template <typename U>
  SharedPtr(SharedPtr<U>& other, U* ptr) : ptr(ptr), pcb(other.pcb) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

  template <typename U>
  SharedPtr(SharedPtr& other, U* ptr) : ptr(ptr), pcb(other.pcb) {
    if (pcb) {
      pcb->shared_count++;
    }
  }

  void swap(SharedPtr& other) {
    std::swap(ptr, other.ptr);
    std::swap(pcb, other.pcb);
  }

  template <typename U> 
  SharedPtr& operator=(const SharedPtr<U>& other) {
    SharedPtr temp(other);
    swap(temp);
    return *this;
  }

  SharedPtr& operator=(const SharedPtr& other) {
    if (this == &other) {
      return *this;
    }
    SharedPtr temp(other);
    swap(temp);
    return *this;
  }

  template <typename U>
  SharedPtr& operator=(SharedPtr<U>&& other) {
    SharedPtr temp(std::move(other));
    swap(temp);
    return *this;
  }

  SharedPtr& operator=(SharedPtr&& other) {
    SharedPtr temp(std::move(other));
    swap(temp);
    return *this;
  }

  ~SharedPtr() {
    if (pcb == nullptr) {
      return;
    }
    pcb->shared_count--;
    if (pcb->shared_count == 0) {
      pcb->destroy_obj();
      if (pcb->weak_count == 0) {
        pcb->destroy_block();
      }
    }
  }

  size_t use_count() const { return pcb->shared_count; }

  void reset() {
    *this = SharedPtr<T>();
  }

  template <typename U>
  void reset(U* ptr = nullptr) {
    SharedPtr<U> temp(ptr);
    swap(temp);
  }

  T* get() const {
    return ptr;
  }

  T& operator*() const {
    return *ptr;
  }

  T* operator->() const {
    return ptr;
  }
};

template <typename U, typename Alloc, typename... Args>
SharedPtr<U> allocateShared(const Alloc& alloc, Args&&...args) {

  using AllocTraits = std::allocator_traits<Alloc>;
  using BlockAlloc = typename AllocTraits::template rebind_alloc<ControlBlockMakeShared<U, Alloc>>;
  using BlockAllocTraits = std::allocator_traits<BlockAlloc>;

  BlockAlloc balloc = alloc;
  auto* pcb = BlockAllocTraits::allocate(balloc, 1);
  BlockAllocTraits::construct(balloc, pcb, alloc, std::forward<Args>(args)...);

  return SharedPtr<U>(pcb);
}

template <typename U, typename... Args>
SharedPtr<U> makeShared(Args&&...args) {
  return allocateShared<U, std::allocator<U>>(std::allocator<U>(),
                                              std::forward<Args>(args)...);
}

template <typename T>
class WeakPtr {
private:

  template <typename>
  friend class WeakPtr;

  template <typename U>
  friend class SharedPtr;

  T* ptr = nullptr;
  BaseControlBlock* pcb = nullptr;

 public:

  WeakPtr() = default;

  WeakPtr(const SharedPtr<T>& other) : ptr(other.ptr), pcb(other.pcb) {
    if(pcb) {
      pcb->weak_count++;
    }
  }

  template <typename U>
  WeakPtr(const SharedPtr<U>& other) : ptr(other.ptr), pcb(other.pcb) {
    if(pcb) {
      pcb->weak_count++;
    }
  }

  WeakPtr(const WeakPtr& other) : ptr(other.ptr), pcb(other.pcb) {
    if(pcb) {
      pcb->weak_count++;
    }
  }

  template <typename U>
  WeakPtr(const WeakPtr<U>& other) : ptr(other.ptr), pcb(other.pcb) {
    if(pcb) {
      pcb->weak_count++;
    }
  }  

  WeakPtr(WeakPtr&& other) : ptr(other.ptr), pcb(other.pcb) {
    other.pcb = nullptr;
    other.ptr = nullptr;
  }

  template <typename U>
  WeakPtr(SharedPtr<U>&& other) : ptr(other.ptr), pcb(other.pcb) {
    other.pcb = nullptr;
    other.ptr = nullptr;
  }

  void swap(WeakPtr& other) {
    std::swap(ptr, other.ptr);
    std::swap(pcb, other.pcb);
  }

  template <typename U> 
  WeakPtr& operator=(const WeakPtr<U>& other) {
    if (this == &other) {
      return this;
    }
    WeakPtr temp(other);
    swap(temp);
    return *this;
  }

  WeakPtr& operator=(const WeakPtr& other) {
    if (this == &other) {
      return this;
    }
    WeakPtr temp(other);
    swap(temp);
    return *this;
  }

  template <typename U>
  WeakPtr& operator=(WeakPtr<U>&& other) {
    WeakPtr temp = std::move(other);
    swap(temp);
    return *this;
  }

  WeakPtr& operator=(WeakPtr&& other) {
    WeakPtr temp = std::move(other);
    swap(temp);
    return *this;
  }

  size_t use_count() const {
    return pcb ? pcb->shared_count : 0;
  }

  bool expired() const {
    return use_count() == 0;
  }

  SharedPtr<T> lock() const {
    return expired() ? SharedPtr<T>() : SharedPtr<T>(*this);
  }

  ~WeakPtr() {
    if (pcb == nullptr) {
      return;
    }
    pcb->weak_count--;
    if (pcb->weak_count == 0 && pcb->shared_count == 0) {
      pcb->destroy_block();
    }
  }
};