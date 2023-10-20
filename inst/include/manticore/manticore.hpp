#ifndef MANTICORE_MANTICORE_HPP
#define MANTICORE_MANTICORE_HPP

#include <mutex>
#include <condition_variable>
#include <functional>
#include <string>
#include <stdexcept>

/**
 * @file manticore.hpp
 * @brief Defines the **manticore** namespace. 
 */

/**
 * @namespace manticore
 * @brief Defines **manticore** classes and functions.
 */
namespace manticore {

/**
 * @brief Execute arbitrary functions on the main thread.
 *
 * An instance of this class should be created on the main thread and initialized with `initialize()`.
 * A corresponding number of worker threads can then be created using `std::thread`.
 * Each worker can request functions for main thread execution via `run()`.
 * The main thread should call `listen()` to wait for such requests, blocking until all workers have finished by calling `finish_thread()`.
 *
 * Outside of parallel contexts, the same class can be used to execute functions directly on the main thread by just calling `run()`.
 * It is not necessary to call `initialize()` beforehand, nor is it required to `listen()`, or to call `finish_thread()`.
 * This is helpful to enable the same code to be used in parallel and serial contexts.
 */
class Executor {
    std::mutex run_lock;
    std::condition_variable cv;

    size_t nthreads;
    size_t ncomplete;
    std::string fallback_error;
    std::string error_message;

    enum class Status : char { FREE, PRIMED, FINISHED };
    Status status;
    std::function<void()> fun;

    bool initialized = false;
    bool done() const { 
        return ncomplete == nthreads;
    }

public:
    /**
     * Initialize a parallel execution session where some tasks must be run on the main thread.
     *
     * @param n Number of worker threads in this session.
     * @param e Default error message if a non-standard exception is thrown.
     */
    void initialize(size_t n, std::string e) {
        nthreads = n;
        ncomplete = 0;
        fallback_error = std::move(e);
        error_message.clear();
        status = Status::FREE;
        initialized = true;
    }

    /**
     * Initialize a parallel execution session where some tasks must be run on the main thread.
     *
     * @param n Number of worker threads in this session.
     */
    void initialize(size_t n) {
        initialize(n, "failed main thread execution");
    }

    /**
     * Indicate that one of the worker threads has finished its execution.
     * This is thread-safe and should be called by each worker thread upon completion.
     * 
     * @param notify Whether to notify the main thread that execution has finished on a worker thread.
     * Set to `false` for more efficiency if this method is called before `listen()`.
     */
    void finish_thread(bool notify = true) {
        // Lock everything before bumping ncomplete to avoid problems with
        // simultaneous writes. We use a lock rather than an atomic variable to
        // ensure that the change propagates correctly to the main thread when
        // it calls done() after notification, otherwise the order of events
        // across multiple threads is not guaranteed.
        {
            std::lock_guard lck(run_lock);
            ++ncomplete;
        }

        if (notify) {
            cv.notify_all(); // possibly trigger loop exit in listen(). 
        }
    }

public:
    /**
     * Make a request to the main thread to run the specified function.
     * 
     * If `initialize()` was previously called on the main thread, this method should be called from a worker thread.
     * It is expected that the main thread is (or will be) running `listen()`.
     *
     * If `initialize()` was not previously called, it is assumed that this method is being called from the main thread.
     * In such cases, `f` is executed immediately.
     *
     * @tparam Function_ Class of the function object, typically a lambda.
     * This should accept no arguments and return no value.
     *
     * @param f Function object or functor to be run on the main thread.
     */
    template<class Function_>
    void run(Function_ f) {
        if (!initialized) {
            f();
            return;
        }

        // Waiting until the main thread executor is free,
        // and then assigning it a task.
        std::unique_lock lk(run_lock);
        cv.wait(lk, [&]{ return status == Status::FREE; });

        fun = std::move(f);
        status = Status::PRIMED;

        // Notifying the main thread that there is a task. Only the
        // main thread waits on PRIMED, so other works should not proceed.
        lk.unlock();
        cv.notify_all();

        lk.lock();
        cv.wait(lk, [&]{ return status == Status::FINISHED; });

        // Making a copy of any error message so we can use it for throwing
        // after the unlock. Also clearing the error message for the next thread.
        auto errcopy = error_message;
        error_message.clear();

        // Unblocking other worker threads, if any are waiting.
        status = Status::FREE;
        lk.unlock();
        cv.notify_all();

        if (!errcopy.empty()) {
            throw std::runtime_error(errcopy);
        }
    }

public:
    /**
     * Listen for `run()` requests from worker threads to execute a function on the main thread.
     * This should be run only on the main thread, and will exit once all worker threads have called `finish_thread()`.
     */
    void listen() {
        while (1) {
            std::unique_lock lk(run_lock);

            cv.wait(lk, [&]{ return status == Status::PRIMED || done(); });
            if (done()) {
                break;
            }

            try {
                fun();
            } catch (std::exception& x) {
                // No throw, we need to make sure we notify the worker of (failed) completion.
                error_message = x.what();
            } catch (...) {
                // For everything else.
                error_message = fallback_error;
            }

            status = Status::FINISHED;

            // Unlock before notifying, see example in https://en.cppreference.com/w/cpp/thread/condition_variable
            lk.unlock();
            cv.notify_all();
        }

        initialized = false;
    }
};

}

#endif
