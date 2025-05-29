#pragma once

#include <deque>

#include <iostream>
#include <memory>
#include <type_traits>

template <typename U, typename V>
using is_base_or_derived =
    std::enable_if_t<std::is_base_of_v<U, V> || std::is_same_v<U, V>>;

template <typename T>
class WeakPtr;

struct BaseControlBlock {
    size_t shared_count = 0;
    size_t weak_count = 0;

    virtual void destroy() = 0;
    virtual void deallocate() = 0;
    virtual ~BaseControlBlock() = default;
};

template <typename T, typename Deleter, typename Alloc>
struct ControlBlockRegular : BaseControlBlock {
    T* ptr = nullptr;
    Deleter deleter;
    Alloc alloc;

    ControlBlockRegular(T* pt, const Deleter& deleter, const Alloc& alloc)
        : ptr(pt), deleter(deleter), alloc(alloc) {}

    void destroy() override {
        deleter(ptr);
        ptr = nullptr;
    }

    void deallocate() override {
        using allocTraits =
            typename std::allocator_traits<Alloc>::template rebind_traits<
                ControlBlockRegular<T, Deleter, Alloc>>;
        using allocType =
            typename std::allocator_traits<Alloc>::template rebind_alloc<
                ControlBlockRegular<T, Deleter, Alloc>>;
        allocType alloc_copy = alloc;
        allocTraits::deallocate(alloc_copy, this, 1);
    }
};

template <typename T, typename Alloc = std::allocator<T>>
struct ControlBlockMakeShared : public BaseControlBlock {
    T ptr;
    Alloc alloc;

    template <typename... Args>
    ControlBlockMakeShared(Alloc alloc, Args&&... args)
        : ptr(std::forward<Args>(args)...), alloc(alloc) {}

    void destroy() override {
        using allocTraits = typename std::allocator_traits<
            Alloc>::template rebind_traits<ControlBlockMakeShared<Alloc>>;
        using allocType = typename std::allocator_traits<
            Alloc>::template rebind_alloc<ControlBlockMakeShared<Alloc>>;
        allocType alloc_copy = alloc;
        allocTraits::destroy(alloc_copy, &ptr);
    }

    void deallocate() override {
        using allocTraits = typename std::allocator_traits<
            Alloc>::template rebind_traits<ControlBlockMakeShared<T, Alloc>>;
        using allocType = typename std::allocator_traits<
            Alloc>::template rebind_alloc<ControlBlockMakeShared<T, Alloc>>;
        allocType alloc_copy = alloc;
        allocTraits::deallocate(alloc_copy, this, 1);
    }
};

template <typename T>
class SharedPtr {
  private:
    template <typename Alloc>
    SharedPtr(ControlBlockMakeShared<T, Alloc>* cb) : ptr(nullptr), cb(cb) {
        if (cb != nullptr) {
            ++(cb->shared_count);
        }
    }

    SharedPtr(const WeakPtr<T>& p) : ptr(p.ptr), cb(p.cb) {
        if (cb != nullptr) {
            ++(cb->shared_count);
        }
    }

    T* ptr = nullptr;
    BaseControlBlock* cb = nullptr;

  public:
    template <typename U, typename = is_base_or_derived<T, U>>
    void swap(SharedPtr<U>& other) {
        std::swap(ptr, other.ptr);
        std::swap(cb, other.cb);
    }

    template <typename U, typename... Args>
    friend SharedPtr<U> makeShared(Args&&... args);

    template <typename U, typename Alloc, typename... Args>
    friend SharedPtr<U> allocateShared(const Alloc& alloc, Args&&... args);

    template <typename U>
    friend class WeakPtr;

    template <typename U>
    friend class SharedPtr;

    SharedPtr(){};

    template <typename U, typename Deleter = std::default_delete<T>,
              typename Alloc = std::allocator<T>>
    SharedPtr(U* ptr, Deleter deleter = Deleter(), Alloc alloc = Alloc())
        : ptr(ptr) {
        using ControlBlockAllocator =
            typename std::allocator_traits<Alloc>::template rebind_alloc<
                ControlBlockRegular<T, Deleter, Alloc>>;
        using ControlBlockAllocatorTraits =
            typename std::allocator_traits<Alloc>::template rebind_traits<
                ControlBlockRegular<T, Deleter, Alloc>>;

        ControlBlockAllocator controlBlockAlloc = alloc;
        auto pt = ControlBlockAllocatorTraits::allocate(controlBlockAlloc, 1);
        new (pt) ControlBlockRegular<T, Deleter, Alloc>(ptr, deleter, alloc);

        cb = pt;
        ++(cb->shared_count);
    }

    SharedPtr(const SharedPtr& other) : ptr(other.ptr), cb(other.cb) {
        if (other.cb != nullptr) {
            ++(cb->shared_count);
        }
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    SharedPtr(const SharedPtr<U>& other) : ptr(other.ptr), cb(other.cb) {
        if (other.cb != nullptr) {
            ++(cb->shared_count);
        }
    }

    SharedPtr(SharedPtr&& other) : ptr(other.ptr), cb(other.cb) {
        other.ptr = nullptr;
        other.cb = nullptr;
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    SharedPtr(SharedPtr<U>&& other) : ptr(other.ptr), cb(other.cb) {
        other.ptr = nullptr;
        other.cb = nullptr;
    }

    SharedPtr& operator=(const SharedPtr& other) {
        if (this == &other) {
            return *this;
        }
        SharedPtr<T> copy = SharedPtr(other);
        swap(copy);
        return *this;
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    SharedPtr& operator=(const SharedPtr<U>& other) {
        SharedPtr<T> copy(other);
        swap(copy);
        return *this;
    }

    SharedPtr& operator=(SharedPtr&& other) {
        SharedPtr copy(std::move(other));
        swap(copy);
        return *this;
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    SharedPtr& operator=(SharedPtr<U>&& other) {
        SharedPtr<T> copy(std::move(other));
        swap(copy);
        return *this;
    }

    size_t use_count() const {
        if (cb == nullptr) {
            return 0;
        }
        return cb->shared_count;
    }

    void reset() {
        SharedPtr<T>().swap(*this);
    }

    template <typename U, typename Deleter = std::default_delete<T>,
              typename Alloc = std::allocator<T>>
    void reset(U* ptr, Deleter del = Deleter(), Alloc alloc = Alloc()) {
        SharedPtr<T>(ptr, del, alloc).swap(*this);
    }

    T& operator*() const {
        if (ptr != nullptr) {
            return *ptr;
        }
        return static_cast<ControlBlockMakeShared<T>*>(cb)->ptr;
    }

    T* operator->() const {
        if (ptr != nullptr) {
            return ptr;
        }
        return &(static_cast<ControlBlockMakeShared<T>*>(cb)->ptr);
    }

    T* get() const {
        return ptr;
    }

    ~SharedPtr() {
        if (cb == nullptr) {
            return;
        }
        --(cb->shared_count);
        if (cb->shared_count == 0) {
            cb->destroy();
            if (cb->weak_count == 0) {
                cb->deallocate();
            }
        }
    }
};

template <typename T>
class WeakPtr {
  private:
    T* ptr = nullptr;
    BaseControlBlock* cb = nullptr;

  public:
    void swap(WeakPtr<T>& other) {
        std::swap(ptr, other.ptr);
        std::swap(cb, other.cb);
    }

    template <typename U>
    friend class WeakPtr;

    template <typename U>
    friend class SharedPtr;

    WeakPtr(){};

    WeakPtr(const SharedPtr<T>& shared) : ptr(shared.ptr), cb(shared.cb) {
        if (cb != nullptr) {
            ++cb->weak_count;
        }
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    WeakPtr(const SharedPtr<U>& shared) : ptr(shared.ptr), cb(shared.cb) {
        if (cb != nullptr) {
            ++cb->weak_count;
        }
    }

    WeakPtr(const WeakPtr& other) : ptr(other.ptr), cb(other.cb) {
        if (cb != nullptr) {
            ++cb->weak_count;
        }
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    WeakPtr(const WeakPtr<U>& other) : ptr(other.ptr), cb(other.cb) {
        if (cb != nullptr) {
            ++cb->weak_count;
        }
    }

    WeakPtr(WeakPtr&& other) : ptr(other.ptr), cb(other.cb) {
        other.ptr = nullptr;
        other.cb = nullptr;
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    WeakPtr(WeakPtr<U>&& other) : ptr(other.ptr), cb(other.cb) {
        other.ptr = nullptr;
        other.cb = nullptr;
    }

    WeakPtr& operator=(const WeakPtr& other) {
        WeakPtr<T> copy(other);
        swap(copy);
        return *this;
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    WeakPtr& operator=(const WeakPtr<U>& other) {
        WeakPtr<T> copy(other);
        swap(copy);
        return *this;
    }

    WeakPtr& operator=(WeakPtr&& other) {
        WeakPtr<T> copy(std::move(other));
        swap(copy);
        return *this;
    }

    template <typename U, typename = is_base_or_derived<T, U>>
    WeakPtr& operator=(WeakPtr<U>&& other) {
        WeakPtr<T> copy(std::move(other));
        swap(copy);
        return *this;
    }

    size_t use_count() const {
        if (cb == nullptr) {
            return 0;
        }
        return cb->shared_count;
    }

    bool expired() const {
        return use_count() == 0;
    }

    SharedPtr<T> lock() const {
        if (expired()) {
            return SharedPtr<T>();
        }
        return SharedPtr<T>(*this);
    }

    ~WeakPtr() {
        if (cb == nullptr) {
            return;
        }
        --(cb->weak_count);
        if (cb->weak_count == 0 && cb->shared_count == 0) {
            cb->deallocate();
        }
    }
};

template <typename U, typename Alloc, typename... Args>
SharedPtr<U> allocateShared(const Alloc& alloc, Args&&... args) {
    using SharedAlloc = typename std::template allocator_traits<
        Alloc>::template rebind_alloc<ControlBlockMakeShared<U, Alloc>>;
    using SharedAllocTraits =
        typename std::template allocator_traits<SharedAlloc>;
    SharedAlloc sharedAlloc = alloc;
    auto* pt = SharedAllocTraits::allocate(sharedAlloc, 1);
    SharedAllocTraits::construct(sharedAlloc, pt, std::move(alloc),
                                 std::forward<Args>(args)...);
    return SharedPtr<U>(pt);
}

template <typename U, typename... Args>
SharedPtr<U> makeShared(Args&&... args) {
    return allocateShared<U>(std::allocator<U>(), std::forward<Args>(args)...);
}