/**
 * @file ThreadSafeAccumulator.hpp
 */


#ifndef PROJET_IPS1_THREADSAFEACCUMULATOR_HPP
#define PROJET_IPS1_THREADSAFEACCUMULATOR_HPP

#include <list>

/**
 * Implemented operators
 */
enum class operation_type {
    Add,
    Mul,
    Nop,
};

/**
 * @class ThreadSafeAccumulator
 * Thread safe accumulator to make commutative operations on a single buffer/result
 *
 * If the accumulator is not accessible (if we are making expensive operations (matrix multiplication ?)), the argument is buffered
 * and will be computed later.
 * @tparam T
 */
template<typename T>
class ThreadSafeAccumulator {
public:
    typedef T (* acc_function)(T accumulator, T arg);
    /**
     * Constructor to use
     * @param acc the initial value, like an empty matrix, or whatever
     * @param type commutative operation that is going to be performed
     */
    ThreadSafeAccumulator(const T& acc, operation_type type);

    ThreadSafeAccumulator(const T& acc, acc_function function);

    /**
     * Computes buffered ops if needed
     * @return the result of the build
     */
    T GetResult();

    /**
     * Pushes a value to the accumulator.
     * @param arg
     */
    inline void push(const T& arg);

private:
    inline void acquire_lock();

    inline void release_lock();

    inline void acquire_buffer_lock();

    inline void release_buffer_lock();

    inline void exec_operation(const T& arg);

    const acc_function fun = nullptr;
    T result;
    const operation_type op = operation_type::Nop;
    std::list<T> buffer{};
    volatile bool result_locked = false;
    volatile bool buffer_locked = false;
};

template<typename T>
T ThreadSafeAccumulator<T>::GetResult()
{
    for (auto r : buffer) {
        exec_operation(r);
    }
    return result;
}

template<typename T>
ThreadSafeAccumulator<T>::ThreadSafeAccumulator(const T& acc, operation_type type)
        : result(acc), op(type) { }

template<typename T>
ThreadSafeAccumulator<T>::ThreadSafeAccumulator(const T& acc, acc_function function)
        : result(acc), fun(function) { }

template<typename T>
inline void ThreadSafeAccumulator<T>::push(const T& arg)
{
    if (result_locked) {
        acquire_buffer_lock();
     //   std::cout << "buffering " << std::endl;
        buffer.push_back(arg);
        release_buffer_lock();
    }
    else {
        acquire_lock();
        exec_operation(arg);
        release_lock();
    }
}

template<typename T>
inline void ThreadSafeAccumulator<T>::acquire_lock()
{
    while (result_locked) {// we wait
    }
    result_locked = true;
}

template<typename T>
inline void ThreadSafeAccumulator<T>::acquire_buffer_lock()
{
    while (buffer_locked) {// we wait
    }
    buffer_locked = true;
}

template<typename T>
inline void ThreadSafeAccumulator<T>::release_lock()
{
    result_locked = false;
}

template<typename T>
inline void ThreadSafeAccumulator<T>::release_buffer_lock()
{
    buffer_locked = false;
}

/**
 * We must ensure the locks are properly placed etc
 * @tparam T
 * @param arg
 */
template<typename T>
inline void ThreadSafeAccumulator<T>::exec_operation(const T& arg)
{
    if (this->op==operation_type::Add) {
        result += arg;
    }
    else if (this->op==operation_type::Mul) {
        result *= arg;
    }
    else if (fun) {
        result = fun(result, arg);
    }
    else if (this->op==operation_type::Nop) {
        std::cerr << "No op defined in accumulator" << std::endl;
    }
}

#endif //PROJET_IPS1_THREADSAFEACCUMULATOR_HPP
