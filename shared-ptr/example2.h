
#include <iostream>
#include <memory>

template <typename U>
class WeakPtr;

template <typename T>
class SharedPtr {

    template <typename U, typename Alloc, typename... Args>
    friend SharedPtr<U> allocateShared(const Alloc&, Args&&...);

    template <typename U, typename... Args>
    friend SharedPtr<U> makeShared(Args&&...);

    template <typename U>
    friend class SharedPtr;

    template <typename U>
    friend class WeakPtr;

private:
    struct BaseControlBlock {
        size_t shared_count = 0;
        size_t weak_count = 0;

        virtual T* get_ptr() = 0;

        virtual void delete_ptr() = 0;

        virtual void destroy_block() = 0;

        virtual ~BaseControlBlock() = default;
    };

    template <typename U, typename Deleter = std::default_delete<U>, typename Alloc = std::allocator<U>>
    struct ControlBlockSimple : BaseControlBlock {
        using AllocTraits = std::allocator_traits<Alloc>;
        using BlockAllocType = typename AllocTraits::template rebind_alloc<ControlBlockSimple<U, Deleter, Alloc>>;
        using BlockAllocTraits = std::allocator_traits<BlockAllocType>;

        T* ptr;
        Deleter deleter;
        BlockAllocType block_alloc;

        ControlBlockSimple(T* ptr, Deleter deleter, BlockAllocType block_alloc)
                : ptr(ptr), deleter(deleter), block_alloc(block_alloc) {}

        T* get_ptr() override {
            return ptr;
        }

        void delete_ptr() override {
            deleter(ptr);
        }

        void destroy_block() override {
            BlockAllocType temp = block_alloc;
            deleter.~Deleter();
            block_alloc.~BlockAllocType();
            BlockAllocTraits::deallocate(temp, this, 1);
        }

        ~ControlBlockSimple() = default;
    };

    template <typename U, typename Alloc = std::allocator<U>>
    struct ControlBlockMakeShared : BaseControlBlock {
        using AllocTraits = std::allocator_traits<Alloc>;
        using BlockAllocType = typename AllocTraits::template rebind_alloc<ControlBlockMakeShared<U, Alloc>>;
        using BlockAllocTraits = std::allocator_traits<BlockAllocType>;

        T object;
        Alloc alloc;
        BlockAllocType block_alloc;

        template <typename... Args>
        ControlBlockMakeShared(Alloc alloc, BlockAllocType block_alloc, Args&&... args)
                : object(std::forward<Args>(args)...), alloc(alloc), block_alloc(block_alloc) {}

        T* get_ptr() override {
            return &object;
        }

        void delete_ptr() override {
            AllocTraits::destroy(alloc, &object);
        }

        void destroy_block() override {
            BlockAllocType temp = block_alloc;
            alloc.~Alloc();
            block_alloc.~BlockAllocType();
            BlockAllocTraits::deallocate(temp, this, 1);
        }

        ~ControlBlockMakeShared() = default;
    };

    BaseControlBlock* cb = nullptr;

    struct MakeSharedConstructorTag {};

    void plus_counter() {
        if (cb) {
            cb->shared_count++;
        }
    }

    SharedPtr(MakeSharedConstructorTag, BaseControlBlock* cb) : cb(cb) {
        plus_counter();
    }

public:
    template <typename U, typename Deleter = std::default_delete<U>, typename Alloc = std::allocator<U>>
    SharedPtr(U* ptr, Deleter deleter = Deleter(), const Alloc& alloc = Alloc());

    SharedPtr(const SharedPtr& shar);

    SharedPtr& operator=(const SharedPtr& shar);

    SharedPtr(SharedPtr&& shar);

    SharedPtr& operator=(SharedPtr&& shar);

    template <typename U>
    SharedPtr(const SharedPtr<U>& shar);

    void swap(SharedPtr& shar);

    template <typename U>
    SharedPtr& operator=(const SharedPtr<U>& shar);

    template <typename U>
    SharedPtr(SharedPtr<U>&& shar);

    template <typename U>
    SharedPtr& operator=(SharedPtr<U>&& shar);

    SharedPtr() = default;

    template <typename U>
    SharedPtr(const WeakPtr<U>& weak);

    ~SharedPtr();

    void reset();

    template <typename U>
    void reset(U* pointer);

    size_t use_count() const;

    T* get() const;

    T* operator->() const;

    T& operator*() const;
};

template <typename T>
template <typename U, typename Deleter, typename Alloc>
SharedPtr<T>::SharedPtr(U* ptr, Deleter deleter, const Alloc& alloc) {
    using AllocTraits = std::allocator_traits<Alloc>;
    using BlockAllocType = typename AllocTraits::template rebind_alloc<ControlBlockSimple<U, Deleter, Alloc>>;
    using BlockAllocTraits = std::allocator_traits<BlockAllocType>;

    BlockAllocType block_alloc = alloc;
    ControlBlockSimple<U, Deleter, Alloc>* cur = BlockAllocTraits::allocate(block_alloc, 1);
    new (cur) ControlBlockSimple<U, Deleter, Alloc>(ptr, deleter, block_alloc);
    cb = cur;
    plus_counter();
}

template <typename T>
SharedPtr<T>::SharedPtr(const SharedPtr& shar) : cb(reinterpret_cast<BaseControlBlock*>(shar.cb)) {
    plus_counter();
}

template <typename T>
SharedPtr<T>& SharedPtr<T>::operator=(const SharedPtr<T>& shar) {
    SharedPtr<T> copy = shar;
    copy.swap(*this);
    return *this;
}

template <typename T>
SharedPtr<T>::SharedPtr(SharedPtr&& shar) : cb(reinterpret_cast<BaseControlBlock*>(shar.cb)) {
    shar.cb = nullptr;
}

template <typename T>
SharedPtr<T>& SharedPtr<T>::operator=(SharedPtr<T>&& shar) {
    SharedPtr copy = std::move(shar);
    copy.swap(*this);
    return *this;
}

template <typename T>
template <typename U>
SharedPtr<T>::SharedPtr(const SharedPtr<U>& shar)
    : cb(reinterpret_cast<BaseControlBlock*>(shar.cb)) {
    plus_counter();
}

template <typename T>
void SharedPtr<T>::swap(SharedPtr& shar) {
    std::swap(cb, shar.cb);
}

template <typename T>
template <typename U>
SharedPtr<T>& SharedPtr<T>::operator=(const SharedPtr<U>& shar) {
    SharedPtr<T> copy = shar;
    copy.swap(*this);
    return *this;
}

template <typename T>
template <typename U>
SharedPtr<T>::SharedPtr(SharedPtr<U>&& shar) : cb(reinterpret_cast<BaseControlBlock*>(shar.cb)) {
    shar.cb = nullptr;
}

template <typename T>
template <typename U>
SharedPtr<T>& SharedPtr<T>::operator=(SharedPtr<U>&& shar) {
    SharedPtr copy = std::move(shar);
    copy.swap(*this);
    return *this;
}

template <typename T>
template <typename U>
SharedPtr<T>::SharedPtr(const WeakPtr<U>& weak) : cb(reinterpret_cast<BaseControlBlock*>(weak.cb)) {
    plus_counter();
}

template <typename T>
SharedPtr<T>::~SharedPtr() {
    if (!cb || !cb->shared_count) return;
    cb->shared_count--;
    if (!cb->shared_count) {
        cb->delete_ptr();
    }
    if (!cb->shared_count && !cb->weak_count) {
        cb->destroy_block();
    }
}

template <typename T>
void SharedPtr<T>::reset() {
    *this = SharedPtr();
}

template <typename T>
template <typename U>
void SharedPtr<T>::reset(U* pointer) {
    *this = SharedPtr(pointer);
}

template <typename T>
size_t SharedPtr<T>::use_count() const {
    return (cb ? cb->shared_count : 0);
}

template <typename T>
T* SharedPtr<T>::get() const {
    return cb ? cb->get_ptr() : nullptr;
}

template <typename T>
T* SharedPtr<T>::operator->() const {
    return get();
}

template <typename T>
T& SharedPtr<T>::operator*() const {
    return *get();
}

template <typename T>
class WeakPtr {

private:
    template <typename U>
    friend class WeakPtr;

    template <typename U>
    friend class SharedPtr;

    typename SharedPtr<T>::BaseControlBlock* cb = nullptr;

    void plus_counter() {
        if (cb) {
            cb->weak_count++;
        }
    }

public:
    template <typename U>
    WeakPtr(const SharedPtr<U>& shar);

    WeakPtr() = default;

    template <typename U>
    WeakPtr(const WeakPtr<U>& weak);

    void swap(WeakPtr<T>& weak);

    template <typename U>
    WeakPtr& operator=(const WeakPtr<U>& weak);

    template <typename U>
    WeakPtr(WeakPtr<U>&& weak);

    template <typename U>
    WeakPtr& operator=(WeakPtr<U>&& weak);

    WeakPtr(const WeakPtr& weak);

    WeakPtr& operator=(const WeakPtr& weak);

    WeakPtr(WeakPtr&& weak);

    WeakPtr& operator=(WeakPtr&& weak);

    void reset();

    size_t use_count() const;

    bool expired() const;

    SharedPtr<T> lock() const;

    ~WeakPtr();
};

template <typename T>
template <typename U>
WeakPtr<T>::WeakPtr(const SharedPtr<U>& shar) : cb(reinterpret_cast<typename SharedPtr<T>::BaseControlBlock*>(shar.cb)) {
    plus_counter();
}

template <typename T>
template <typename U>
WeakPtr<T>::WeakPtr(const WeakPtr<U>& weak) : cb(reinterpret_cast<typename SharedPtr<T>::BaseControlBlock*>(weak.cb)) {
    plus_counter();
}

template <typename T>
void WeakPtr<T>::swap(WeakPtr<T>& weak) {
    std::swap(cb, weak.cb);
}

template <typename T>
template <typename U>
WeakPtr<T>& WeakPtr<T>::operator=(const WeakPtr<U>& weak) {
    WeakPtr copy = weak;
    copy.swap(*this);
    return *this;
}

template <typename T>
template <typename U>
WeakPtr<T>::WeakPtr(WeakPtr<U>&& weak) : cb(weak.cb) {
    weak.cb = nullptr;
}

template <typename T>
template <typename U>
WeakPtr<T>& WeakPtr<T>::operator=(WeakPtr<U>&& weak) {
    WeakPtr copy = std::move(weak);
    copy.swap(*this);
    return *this;
}

template <typename T>
WeakPtr<T>::WeakPtr(const WeakPtr& weak) : cb(reinterpret_cast<typename SharedPtr<T>::BaseControlBlock*>(weak.cb)) {
    plus_counter();
}

template <typename T>
WeakPtr<T>& WeakPtr<T>::operator=(const WeakPtr& weak) {
    WeakPtr copy = weak;
    copy.swap(*this);
    return *this;
}

template <typename T>
WeakPtr<T>::WeakPtr(WeakPtr&& weak) : cb(weak.cb) {
    weak.cb = nullptr;
}

template <typename T>
WeakPtr<T>& WeakPtr<T>::operator=(WeakPtr&& weak) {
    WeakPtr copy = std::move(weak);
    copy.swap(*this);
    return *this;
}

template <typename T>
void WeakPtr<T>::reset() {
    *this = WeakPtr<T>();
}

template <typename T>
size_t WeakPtr<T>::use_count() const {
    return (cb ? cb->shared_count : 0);
}

template <typename T>
bool WeakPtr<T>::expired() const {
    return use_count() == 0;
}

template <typename T>
SharedPtr<T> WeakPtr<T>::lock() const {
    return expired() ? SharedPtr<T>() : SharedPtr<T>(*this);
}

template <typename T>
WeakPtr<T>::~WeakPtr() {
    if (!cb || !cb->weak_count) return;
    cb->weak_count--;
    if (cb->shared_count + cb->weak_count == 0) cb->destroy_block();
}