#pragma once

#include <memory>
#include <functional>
#include <stdexcept>
#include <atomic>

class ControlBlock {
protected:
    size_t shared_counter_ = 0;
    size_t weak_counter_ = 0;

    virtual void destroy_field() = 0;
    virtual void destroy_this() = 0;

public:
    void inc_shared() {
        ++shared_counter_;
    }

    void inc_weak() {
        ++weak_counter_;
    }

    void dec_shared() {
        if (shared_counter_ == 0) throw std::logic_error("shared_counter_ < 0");
        --shared_counter_;
        if (shared_counter_ == 0) {
            destroy_field();
        }
        if (shared_counter_ == 0 && weak_counter_ == 0) {
            destroy_this();
        }
    }

    void dec_weak() {
        if (weak_counter_ == 0) throw std::logic_error("weak_counter_ < 0");
        --weak_counter_;
        if (shared_counter_ == 0 && weak_counter_ == 0) {
            destroy_this();
        }
    }

    size_t use_count() const {
        return shared_counter_;
    }

    virtual ~ControlBlock() = default;
};

template <typename T, typename Alloc = std::allocator<T>>
class ControlBlockAllocateShared: public ControlBlock {
private:
    using AllocThis = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockAllocateShared<T, Alloc>>;

    alignas(T) char value[sizeof(T)];
    Alloc alloc;
    AllocThis allocThisCopy;

    void destroy_field() override {
        std::allocator_traits<Alloc>::destroy(alloc, ptr);
    }

    void destroy_this() override {
        auto allocThis = allocThisCopy;
        this->~ControlBlockAllocateShared();
        std::allocator_traits<AllocThis>::deallocate(allocThis, this, 1);
    }

public:
    T* ptr;

    template <typename... Args>
    ControlBlockAllocateShared(Alloc alloc, AllocThis allocThis, Args&&... args): alloc(alloc), allocThisCopy(allocThis) {
        ptr = reinterpret_cast<T*>(value);
        std::allocator_traits<Alloc>::construct(alloc, ptr, std::forward<Args>(args)...);
    }

    ~ControlBlockAllocateShared() override = default;
};

template <typename T, typename Alloc = std::allocator<T>, typename Deleter = std::default_delete<T>>
class ControlBlockPointer: public ControlBlock {
private:
    using AllocThis = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockPointer<T, Alloc, Deleter>>;

    T* ptr = nullptr;
    Deleter d;
    AllocThis allocThisCopy;

    void destroy_field() override {
        d(ptr);
        ptr = nullptr;
    }

    void destroy_this() override {
        auto allocThis = allocThisCopy;
        this->~ControlBlockPointer();
        std::allocator_traits<AllocThis>::deallocate(allocThis, this, 1);
    }

public:
    ControlBlockPointer(T* ptr, AllocThis alloc = AllocThis(), Deleter d = Deleter()): ptr(ptr), d(d), allocThisCopy(alloc) {}

    ~ControlBlockPointer() override = default;
};

template <typename T>
class WeakPtr;
template <typename T>
class SharedPtr;
template <typename T>
class EnableSharedFromThis;

template <typename T, typename Alloc, typename... Args>
SharedPtr<T> allocateShared(Alloc alloc, Args&&... args);

template <typename T>
class SharedPtr {
private:
    ControlBlock* control_block_ = nullptr;
    T* ptr_ = nullptr;

    SharedPtr(ControlBlock* data, T* ptr): control_block_(data), ptr_(ptr) {
        data->inc_shared();
    }

    template <typename U, typename Alloc, typename... Args>
    friend SharedPtr<U> allocateShared(Alloc alloc, Args&&... args);

    template <typename U>
    friend class WeakPtr;

    template <typename U>
    friend class SharedPtr;

public:
    SharedPtr() = default;

    template <typename U>
    explicit SharedPtr(U* ptr): control_block_(new ControlBlockPointer<U>(ptr)), ptr_(static_cast<T*>(ptr)) {
        if (control_block_) control_block_->inc_shared();
    }

    SharedPtr(const SharedPtr<T>& that): control_block_(that.control_block_), ptr_(that.ptr_) {
        if (control_block_) control_block_->inc_shared();
    }

    SharedPtr(SharedPtr<T>&& that): control_block_(that.control_block_), ptr_(that.ptr_) {
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
    }

    template <typename U>
    SharedPtr(const SharedPtr<U>& other) : ptr_(other.ptr_), control_block_(other.control_block_) {
        if (control_block_) {
            control_block_->inc_shared();
        }
    }

    template <typename U>
    SharedPtr(SharedPtr<U>&& that): control_block_(that.control_block_), ptr_(static_cast<T*>(that.ptr_)) {
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
    }

    template <typename Deleter>
    SharedPtr(T* ptr, Deleter d): SharedPtr(ptr, d, std::allocator<T>()) {}

    template <typename Deleter, typename Alloc>
    SharedPtr(T* ptr, Deleter d, Alloc alloc): ptr_(ptr) {
        using AllocThis = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockPointer<T, Alloc, Deleter>>;
        AllocThis allocThis = alloc;
        auto holder = std::allocator_traits<AllocThis>::allocate(allocThis, 1);
        new (holder) ControlBlockPointer<T, Alloc, Deleter>(ptr, allocThis, d);
        control_block_ = holder;
        control_block_->inc_shared();
    }

    SharedPtr<T>& operator=(const SharedPtr<T>& that) {
        if (this == &that) return *this;
        if (control_block_) control_block_->dec_shared();
        control_block_ = that.control_block_;
        if (control_block_) control_block_->inc_shared();
        ptr_ = that.ptr_;
        return *this;
    }

    SharedPtr<T>& operator=(SharedPtr<T>&& that) {
        if (control_block_) control_block_->dec_shared();
        control_block_ = that.control_block_;
        ptr_ = that.ptr_;
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
        return *this;
    }

    template <typename U> 
    SharedPtr& operator=(const SharedPtr<U>& other) {
        if (this == &other) {
        return *this;
        }
        SharedPtr temp(other);
        swap(temp);
        return *this;
    }

    template <typename U>
    SharedPtr<T>& operator=(SharedPtr<U>&& that) {
        if (control_block_) control_block_->dec_shared();
        control_block_ = that.control_block_;
        ptr_ = static_cast<T*>(that.ptr_);
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
        return *this;
    }

    void swap(SharedPtr<T>& that) {
        std::swap(control_block_, that.control_block_);
        std::swap(ptr_, that.ptr_);
    }

    size_t use_count() const {
        if (control_block_) {
            return control_block_->use_count();
        }
        return 0;
    }

    T& operator*() const {
        return *get();
    }

    T* operator->() const {
        return get();
    }

    void reset(T* nptr = nullptr) {
        if (control_block_) control_block_->dec_shared();
        control_block_ = nullptr;
        ptr_ = nullptr;
        if (!nptr) return;
        *this = SharedPtr(nptr);
    }

    T* get() const {
        return ptr_;
    }

    ~SharedPtr() {
        if (control_block_) control_block_->dec_shared();
    }
};

template <typename T, typename Alloc, typename... Args>
SharedPtr<T> allocateShared(Alloc alloc, Args&&... args) {
    using AllocThis = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockAllocateShared<T, Alloc>>;
    AllocThis allocThis = alloc;
    auto ptr = std::allocator_traits<AllocThis>::allocate(allocThis, 1);
    new (ptr) ControlBlockAllocateShared<T, Alloc>(alloc, allocThis, std::forward<Args>(args)...);
    SharedPtr<T> ans(ptr, ptr->ptr);
    if constexpr (std::is_base_of_v<EnableSharedFromThis<T>, T>) {
        ptr->ptr->ptr_ = ans;
    }
    return ans;
}

template <typename T, typename... Args>
SharedPtr<T> makeShared(Args&&... args) {
    return allocateShared<T, std::allocator<T>>(std::allocator<T>(), std::forward<Args>(args)...);
}

template <typename T>
class WeakPtr {
    ControlBlock* control_block_ = nullptr;
    T* ptr_ = nullptr;

    template <typename U>
    friend class WeakPtr;

public:
    WeakPtr() = default;

    WeakPtr(const WeakPtr<T>& that): control_block_(that.control_block_), ptr_(that.ptr_) {
        if (control_block_) control_block_->inc_weak();
    }

    WeakPtr(WeakPtr<T>&& that): control_block_(that.control_block_), ptr_(that.ptr_) {
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
    }

    template <typename U>
    WeakPtr(const WeakPtr<U>& that): control_block_(that.control_block_), ptr_(static_cast<T*>(that.ptr_)) {
        if (control_block_) control_block_->inc_weak();
    }

    template <typename U>
    WeakPtr(WeakPtr<U>&& that): control_block_(that.control_block_), ptr_(static_cast<T*>(that.ptr_)) {
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
    }

    WeakPtr<T>& operator=(const WeakPtr<T>& that) {
        if (control_block_) control_block_->dec_weak();
        control_block_ = that.control_block_;
        if (control_block_) control_block_->inc_weak();
        ptr_ = that.ptr_;
        return *this;
    }

    WeakPtr<T>& operator=(WeakPtr<T>&& that) {
        if (control_block_) control_block_->dec_weak();
        control_block_ = that.control_block_;
        ptr_ = that.ptr_;
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
        return *this;
    }

    template <typename U>
    WeakPtr<T>& operator=(const WeakPtr<U>& that) {
        if (control_block_) control_block_->dec_weak();
        control_block_ = that.control_block_;
        if (control_block_) control_block_->inc_weak();
        ptr_ = static_cast<T*>(that.ptr_);
        return *this;
    }

    template <typename U>
    WeakPtr<T>& operator=(WeakPtr<U>&& that) {
        if (control_block_) control_block_->dec_weak();
        control_block_ = that.control_block_;
        ptr_ = static_cast<T*>(that.ptr_);
        that.control_block_ = nullptr;
        that.ptr_ = nullptr;
    }

    template <typename U>
    WeakPtr(const SharedPtr<U>& that): control_block_(that.control_block_), ptr_(static_cast<T*>(that.ptr_)) {
        if (control_block_) control_block_->inc_weak();
    }

    bool expired() const {
        if (!control_block_) return true;
        return control_block_->use_count() == 0;
    }

    SharedPtr<T> lock() const {
        if (!control_block_) {
            throw std::runtime_error("this object is not allocated");
        }
        SharedPtr<T> ans(this->control_block_, this->ptr_);
        return ans;
    }

    size_t use_count() const {
        if (!control_block_) return 0;
        return control_block_->use_count();
    }

    ~WeakPtr() {
        if (control_block_) {
            control_block_->dec_weak();
        }
    }
};

template <typename T>
class EnableSharedFromThis {
private:
    template <typename U, typename Alloc, typename... Args>
    friend SharedPtr<U> allocateShared(Alloc alloc, Args&&... args);

    WeakPtr<T> ptr_;

public:
    SharedPtr<T> shared_from_this() const {
        return ptr_.lock();
    }
};