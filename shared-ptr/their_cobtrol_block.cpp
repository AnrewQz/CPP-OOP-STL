#include<iostream>
#include<memory>
#include<type_traits>

struct BaseControlBlock {
  size_t shared_count;
  size_t weak_count;

  virtual void destroy_obj() = 0;
  virtual void destroy_block() = 0;

  size_t use_count() const {
    return shared_count;
  }

  virtual ~BaseControlBlock() = default;
};

template <typename U, typename Alloc = std::allocator<U>, typename Deleter = std::default_delete<U>>
struct ControlBlockRegular: BaseControlBlock {

    using BlockAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockRegular<U, Alloc, Deleter>>;

    U* ptr = nullptr;
    Deleter d;
    BlockAlloc allocThisCopy;

    void destroy_obj() override {
        d(ptr);
        ptr = nullptr;
    }

    void destroy_block() override {
        auto allocThis = allocThisCopy;
        this->~ControlBlockRegular();
        std::allocator_traits<BlockAlloc>::deallocate(allocThis, this, 1);
    }

    ControlBlockRegular(U* ptr, BlockAlloc alloc = BlockAlloc(), Deleter d = Deleter()): ptr(ptr), d(d), allocThisCopy(alloc) {}

    ~ControlBlockRegular() override = default;
};