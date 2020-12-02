/**
 * @file CallbackResultBuilder.hpp
 */


#ifndef PROJET_IPS1_CALLBACKRESULTBUILDER_HPP
#define PROJET_IPS1_CALLBACKRESULTBUILDER_HPP

#include <list>

/**
 * Implemented operators
 */
enum class operation_type {
    Add,
    Mul,
};

/**
 * @class CallbackResultBuilder
 * Thread safe accumulator to make commutative operations on a single buffer/result
 *
 * If the accumulator is not accessible (if we are making expensive operations (matrix multiplication ?)), the argument is buffered
 * and will be computed later.
 * @tparam T
 */
template<typename T>
class CallbackResultBuilder {
public:
    /**
     * Constructor to use
     * @param acc the initial value, like an empty matrix, or whatever
     * @param type commutative operation that is going to be performed
     */
    explicit CallbackResultBuilder(const T &acc, operation_type type);

    /**
     * Computes buffered ops if needed
     * @return the result of the build
     */
    T GetResult();

    /**
     * Pushes a value to the accumulator.
     * @param arg
     */
    inline void push(const T &arg);

private:
    inline void acquire_lock();

    inline void release_lock();

    inline void acquire_buffer_lock();

    inline void release_buffer_lock();

    inline void exec_operation(const T &arg);

    T result;
    operation_type op;
    std::list<T> buffer{};
    volatile bool result_locked = false;
    volatile bool buffer_locked = false;
};


template<typename T>
T CallbackResultBuilder<T>::GetResult() {
    for (auto r : buffer) {
        exec_operation(r);
    }
    return result;
}

template<typename T>
CallbackResultBuilder<T>::CallbackResultBuilder(const T &acc, operation_type type) : result(acc), op(type) {}

template<typename T>
inline void CallbackResultBuilder<T>::push(const T &arg) {
    if (result_locked) {
        acquire_buffer_lock();
        buffer.push_back(arg);
        release_buffer_lock();
    } else {
        acquire_lock();
        exec_operation(arg);
        release_lock();
    }
}

template<typename T>
inline void CallbackResultBuilder<T>::acquire_lock() {
    while (result_locked) {// we wait
    }
    result_locked = true;
}

template<typename T>
inline void CallbackResultBuilder<T>::acquire_buffer_lock() {
    while (buffer_locked) {// we wait
    }
    buffer_locked = true;
}

template<typename T>
inline void CallbackResultBuilder<T>::release_lock() {
    result_locked = false;
}

template<typename T>
inline void CallbackResultBuilder<T>::release_buffer_lock() {
    buffer_locked = false;
}

/**
 * We must ensure the locks are properly placed etc
 * @tparam T
 * @param arg
 */
template<typename T>
inline void CallbackResultBuilder<T>::exec_operation(const T &arg) {
    if (this->op == operation_type::Add) {
        result += arg;
    } else if (this->op == operation_type::Mul) {
        result *= arg;
    } else {
        //TODO CRASH BECAUSE NOT IMPLEMENTED
    }
}


#endif //PROJET_IPS1_CALLBACKRESULTBUILDER_HPP
