#ifndef RATICATE_UNKNOWNMATRIX_HPP
#define RATICATE_UNKNOWNMATRIX_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "SimpleMatrix.hpp"
#include "SparseArraySeed.hpp"
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

#ifdef RATICATE_PARALLELIZE_UNKNOWN
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#endif

namespace raticate {

template<typename Data, typename Index>
struct UnknownMatrixCore {
    UnknownMatrixCore(Rcpp::RObject seed) :
        original_seed(seed),
        delayed_env(Rcpp::Environment::namespace_env("DelayedArray")),
        dense_extractor(delayed_env["extract_array"]),
        sparse_extractor(delayed_env["OLD_extract_sparse_array"])
    {
        // We assume the constructor only occurs on the main thread, so we
        // won't bother locking things up. I'm also not sure that the
        // operations in the initialization list are thread-safe.

        {
            auto base = Rcpp::Environment::base_env();
            Rcpp::Function fun = base["dim"];
            Rcpp::RObject output = fun(seed);
            if (output.sexp_type() != INTSXP) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should return an integer vector");
            }
            Rcpp::IntegerVector dims(output);
            if (dims.size() != 2 || dims[0] < 0 || dims[1] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should contain two non-negative integers");
            }
            nrow = dims[0];
            ncol = dims[1];
        }

        {
            Rcpp::Function fun = delayed_env["is_sparse"];
            Rcpp::LogicalVector is_sparse = fun(seed);
            if (is_sparse.size() != 1) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'is_sparse(<" + ctype + ">)' should return a logical vector of length 1");
            }
            sparse = (is_sparse[0] != 0);
        }

        {
            Rcpp::Function fun = delayed_env["colAutoGrid"];
            Rcpp::RObject output = fun(seed);
            Rcpp::IntegerVector spacing = output.slot("spacings");
            if (spacing.size() != 2 || spacing[1] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'spacings' slot of 'colAutoGrid(<" + ctype + ">)' should contain two non-negative integers");
            }
            block_ncol = spacing[1];
        }

        {
            Rcpp::Function fun = delayed_env["rowAutoGrid"];
            Rcpp::RObject output = fun(seed);
            Rcpp::IntegerVector spacing = output.slot("spacings");
            if (spacing.size() != 2 || spacing[0] < 0) {
                auto ctype = get_class_name(original_seed);
                throw std::runtime_error("'spacings' slot of 'rowAutoGrid(<" + ctype + ">)' should contain two non-negative integers");
            }
            block_nrow = spacing[0];
        }
    }

public:
    Index nrow, ncol;
    bool sparse;

private:
    Index block_nrow, block_ncol;

    Rcpp::RObject original_seed;
    Rcpp::Environment delayed_env;
    Rcpp::Function dense_extractor, sparse_extractor;

public:
    template<bool sparse_>
    struct Workspace {
        Workspace(Index full) {
            secondary_indices = R_NilValue;
            secondary_len = full;
        }

        Workspace(Index start, Index length) {
            secondary_indices = create_consecutive_indices(start, length);
            secondary_len = length;
        }

        Workspace(const std::vector<Index>& indices) {
            Rcpp::IntegerVector temp(indices.begin(), indices.end());
            for (auto& x : temp) { ++x; } // 1-based
            secondary_indices = temp;
            secondary_len = indices.size();
        }

    public:
        Index primary_block_start, primary_block_len, secondary_len;
        Rcpp::RObject secondary_indices;
        std::shared_ptr<tatami::Matrix<Data, Index> > buffer;
        std::shared_ptr<tatami::Extractor<tatami::DimensionSelectionType::FULL, sparse_, Data, Index> > bufextractor;
        Rcpp::RObject contents;
    };

private:
    static std::pair<Index, Index> round_indices(Index i, Index interval, Index max) {
        Index new_first = (i / interval) * interval;
        Index new_last = std::min(max, new_first + interval);
        return std::make_pair(new_first, new_last - new_first);
    }

    static Rcpp::IntegerVector create_consecutive_indices(Index start, Index length) {
        Rcpp::IntegerVector out(length);
        std::iota(out.begin(), out.end(), start + 1); // 1-based
        return out;
    }

    template<bool byrow, bool sparse>
    Rcpp::List create_rounded_indices(size_t i, Workspace<sparse>* work) const {
        Rcpp::List indices(2);
        if constexpr(byrow) {
            auto row_rounded = round_indices(i, block_nrow, nrow);
            work->primary_block_start = row_rounded.first;
            work->primary_block_len = row_rounded.second;
            indices[0] = create_consecutive_indices(row_rounded.first, row_rounded.second);
            indices[1] = work->secondary_indices;
        } else {
            auto col_rounded = round_indices(i, block_ncol, ncol);
            work->primary_block_start = col_rounded.first;
            work->primary_block_len = col_rounded.second;
            indices[0] = work->secondary_indices;
            indices[1] = create_consecutive_indices(col_rounded.first, col_rounded.second);
        }
        return indices;
    }

    template<bool byrow, bool sparse, bool sparse_err = sparse>
    void check_buffered_dims(const tatami::Matrix<Data, Index>* parsed, const Workspace<sparse>* work) const {
        size_t parsed_primary = (byrow ? parsed->nrow() : parsed->ncol());
        size_t parsed_secondary = (byrow ? parsed->ncol() : parsed->nrow());

        if (parsed_primary != work->primary_block_len || parsed_secondary != work->secondary_len) {
            auto ctype = get_class_name(original_seed);
            throw std::runtime_error("'" + 
                (sparse_err ? std::string("extract_sparse_array") : std::string("extract_array")) + 
                "(<" + ctype + ">)' returns incorrect dimensions");
        }
    }

public:
    template<bool sparse>
    static bool needs_reset(Index i, const Workspace<sparse>* work) {
        if (work->buffer != nullptr) {
            if (i >= work->primary_block_start && i < work->primary_block_start + work->primary_block_len) {
                return false;
            }
        }
        return true;
    }

    template<bool byrow>
    void run_dense_extractor(Index i, const tatami::Options& options, Workspace<false>* work) const {
        auto indices = create_rounded_indices<byrow>(i, work);
        Rcpp::RObject val0 = dense_extractor(original_seed, indices);

        auto parsed = parse_simple_matrix<Data, Index>(val0);
        check_buffered_dims<byrow, false>(parsed.matrix.get(), work);

        work->buffer = parsed.matrix;
        work->contents = parsed.contents;
        work->bufextractor = tatami::new_extractor<byrow, false>(work->buffer.get(), options);
    }

    template<bool byrow>
    void run_sparse_extractor(Index i, const tatami::Options& options, Workspace<true>* work) const {
        auto indices = create_rounded_indices<byrow>(i, work);

        if (sparse) {
            auto val0 = sparse_extractor(original_seed, indices);
            auto parsed = parse_SparseArraySeed<Data, Index>(val0);
            check_buffered_dims<byrow, true, true>(parsed.matrix.get(), work);

            work->buffer = parsed.matrix;
            work->contents = parsed.contents;
        } else {
            auto val0 = dense_extractor(original_seed, indices);
            auto parsed = parse_simple_matrix<Data, Index>(val0);
            check_buffered_dims<byrow, true, false>(parsed.matrix.get(), work);

            work->buffer = parsed.matrix;
            work->contents = parsed.contents;
        }

        work->bufextractor = tatami::new_extractor<byrow, true>(work->buffer.get(), options);
    }
};

#ifdef RATICATE_PARALLELIZE_UNKNOWN

/*****************************
 *** Main thread evaluator ***
 *****************************/

template<typename Data, typename Index>
struct UnknownEvaluator {
    bool sparse;
    bool byrow;

    Index index;
    const tatami::Options* options;
    typename UnknownMatrixCore<Data, Index>::template Workspace<false>* dwork;
    typename UnknownMatrixCore<Data, Index>::template Workspace<true>* swork;

    const UnknownMatrixCore<Data, Index>* parent;

    bool parallel = false;
    bool ready = false;
    bool finished = false;
    std::string error;

    bool create_work = false;
    tatami::DimensionSelectionType selection;
    Index full_length, block_start, block_length;
    const std::vector<Index>* indices;
    typename UnknownMatrixCore<Data, Index>::template Workspace<false>** new_dwork;
    typename UnknownMatrixCore<Data, Index>::template Workspace<true>** new_swork;

public:
    template<bool B>
    void set(Index i, const tatami::Options& opt, const UnknownMatrixCore<Data, Index>* core) {
        byrow = B;
        index = i;
        options = &opt;
        parent = core;
        ready = true;
        finished = false;
        create_work = false;
    }

    template<bool B>
    void set(Index i, const tatami::Options& opt, typename UnknownMatrixCore<Data, Index>::template Workspace<false>* w, const UnknownMatrixCore<Data, Index>* core) {
        set<B>(i, opt, core);
        sparse = false;
        dwork = w;
    }

    template<bool B>
    void set(Index i, const tatami::Options& opt, typename UnknownMatrixCore<Data, Index>::template Workspace<true>* w, const UnknownMatrixCore<Data, Index>* core) {
        set<B>(i, opt, core);
        sparse = true;
        swork = w;
    }

public:
    template<bool sparse>
    void set(typename UnknownMatrixCore<Data, Index>::template Workspace<sparse>** nw) {
        if constexpr(sparse) {
            new_swork = nw;
        } else {
            new_dwork = nw;
        }
        ready = true;
        finished = false;
        create_work = true;
    }

    template<bool sparse>
    void set(typename UnknownMatrixCore<Data, Index>::template Workspace<sparse>** nw, Index fl) {
        set(nw);
        selection = tatami::DimensionSelectionType::FULL;
        full_length = fl;
    }

    template<bool sparse>
    void set(typename UnknownMatrixCore<Data, Index>::template Workspace<sparse>** nw, Index bs, Index bl) {
        set(nw);
        selection = tatami::DimensionSelectionType::BLOCK;
        block_start = bs;
        block_length = bl;
    }

    template<bool sparse>
    void set(typename UnknownMatrixCore<Data, Index>::template Workspace<sparse>** nw, const std::vector<Index>& idx) {
        set(nw);
        selection = tatami::DimensionSelectionType::INDEX;
        indices = &idx;
    }

public:
    void harvest() {
        if (create_work) {
            if (!sparse) {
                if (selection == tatami::DimensionSelectionType::FULL) {
                    *new_dwork = new typename UnknownMatrixCore<Data, Index>::template Workspace<false>(full_length);
                } else if (selection == tatami::DimensionSelectionType::BLOCK) {
                    *new_dwork = new typename UnknownMatrixCore<Data, Index>::template Workspace<false>(block_start, block_length);
                } else {
                    *new_dwork = new typename UnknownMatrixCore<Data, Index>::template Workspace<false>(*indices);
                }

            } else {
                if (selection == tatami::DimensionSelectionType::FULL) {
                    *new_swork = new typename UnknownMatrixCore<Data, Index>::template Workspace<true>(full_length);
                } else if (selection == tatami::DimensionSelectionType::BLOCK) {
                    *new_swork = new typename UnknownMatrixCore<Data, Index>::template Workspace<true>(block_start, block_length);
                } else {
                    *new_swork = new typename UnknownMatrixCore<Data, Index>::template Workspace<true>(*indices);
                }
            }
            
        } else {
            if (!sparse) {
                if (byrow) {
                    parent->template run_dense_extractor<true>(index, *options, dwork);
                } else {
                    parent->template run_dense_extractor<false>(index, *options, dwork);
                }

            } else {
                if (byrow) {
                    parent->template run_sparse_extractor<true>(index, *options, swork);
                } else {
                    parent->template run_sparse_extractor<false>(index, *options, swork);
                }
            }
        }
        finished = true;
    }

    void reset() {
        finished = false;
        ready = false;
    }
};

template<typename Data, typename Index>
UnknownEvaluator<Data, Index>& unknown_evaluator() {
    static UnknownEvaluator<Data, Index> e;
    return e;
}

/****************************
 *** Parallel coordinator ***
 ****************************/

struct ParallelCoordinator {
    std::mutex coord_lock;
    std::mutex rcpp_lock;
    std::condition_variable cv;

    template<typename Data, typename Index>
    struct OnMainExit {
        UnknownEvaluator<Data, Index> copy;
        OnMainExit() : copy(unknown_evaluator<Data, Index>()) {}
        ~OnMainExit() {
            auto& ex = unknown_evaluator<Data, Index>();
            ex = copy;
        }
    };

    template<typename Data, typename Index, class Function>
    void run(Function f, size_t n, size_t nthreads) {
        // Acquire the evaluator lock to indicate that we're currently in a single
        // parallel context. This avoids wacky messages from other calls to run().
        std::lock_guard<std::mutex> clk(coord_lock);

        // This restores the state of the unknown evaluator to what it was
        // before entry into this function, to enable nested calls to run().
        auto& ex = unknown_evaluator<Data, Index>();
        OnMainExit<Data, Index> copier;

        // Only parallelize if it's strictly necessary.
        ex.parallel = (n > 1 && nthreads > 1);
        ex.error = "";

        if (ex.parallel) {
            size_t jobs_per_worker = std::ceil(static_cast<double>(n) / nthreads);
            size_t start = 0;
            std::vector<std::thread> jobs;
            std::atomic_size_t ncomplete = 0;
            std::vector<std::string> errors(nthreads);

            for (size_t w = 0; w < nthreads; ++w) {
                size_t end = std::min(n, start + jobs_per_worker);
                if (start >= end) {
                    ncomplete++;
                    continue;
                }

                jobs.emplace_back([&](size_t id, size_t s, size_t l) -> void {
                    try {
                        f(id, s, l);
                    } catch (std::exception& x) {
                        // No throw here, we need to make sure we mark the
                        // thread as being completed so that the main loop can quit.
                        errors[id] = x.what();
                    }
                    ncomplete++;
                    cv.notify_all();
                }, w, start, end - start);

                start += jobs_per_worker;
            }

            // Handling all requests from the workers.
            while (1) {
                std::unique_lock lk(rcpp_lock);
                cv.wait(lk, [&]{ return (ex.ready && !ex.finished) || ncomplete.load() == nthreads; });
                if (ncomplete.load() == nthreads) {
                    break;
                }

                try {
                    ex.harvest();
                } catch (std::exception& x) {
                    // No throw, we need to make sure we notify the worker of (failed) completion.
                    ex.finished = true;
                    ex.error = x.what();
                } catch (...) {
                    // Sometimes R throws these weird Rcpp::LongjumpException errors.
                    ex.finished = true;
                    ex.error = "failed extraction from the unknown matrix";
                }

                // Unlock before notifying, see https://en.cppreference.com/w/cpp/thread/condition_variable
                lk.unlock();
                cv.notify_all();
            }

            for (auto& job : jobs) {
                job.join();
            }

            for (auto err : errors) {
                if (!err.empty()) {
                    throw std::runtime_error(err);
                }
            }
        } else {
            f(0, 0, n);
        }
    }

    template<typename Data, typename Index, class ParallelFunction, class SerialFunction>
    void lock(ParallelFunction pfun, SerialFunction sfun) {
        auto& ex = unknown_evaluator<Data, Index>();
        if (!ex.parallel) {
            sfun();
            return;
        }

        // Waiting until the main thread executor is free,
        // and then assigning it a task.
        {
            std::unique_lock lk(rcpp_lock);
            cv.wait(lk, [&]{ return !ex.ready; });

            // We can throw here because we're not obliged to initiate
            // discussion with the main thread; so it's not like we're
            // going to be leaving the main thread in a hanging state.
            if (!ex.error.empty()) {
                throw std::runtime_error(ex.error);
            }

            try {
                pfun();
            } catch (std::exception& e) {
                // Task assignment failed, so we make sure to reset to
                // avoid partial assignment (and unblock all other workers).
                ex.reset();
                ex.error = e.what();
                throw;
            }
        }

        // Notifying everyone that there is a task. At this point,
        // ready = true and finished = false, so the waiting workers
        // should not proceed; only the main thread should respond.
        cv.notify_all();

        // Checking that we get the finished result, and then we set
        // ready = false to give another worker thread the chance to acquire the lock.
        {
            std::unique_lock lk(rcpp_lock);
            cv.wait(lk, [&]{ return ex.finished; });

            ex.reset();
            if (!ex.error.empty()) {
                // Throwing after the reset so that other workers are not blocked.
                throw std::runtime_error(ex.error);
            }
        }
    }
};

inline ParallelCoordinator& parallel_coordinator() {
    static ParallelCoordinator c;
    return c;
}

#endif

template<typename Data, typename Index>
class UnknownMatrix : public tatami::Matrix<Data, Index> {
public:
    UnknownMatrix(Rcpp::RObject seed) : core(seed) {}

    Index nrow() const {
        return core.nrow;
    }

    Index ncol() const {
        return core.ncol;
    }

    bool sparse() const {
        return core.sparse;
    }

    bool prefer_rows() const {
        // All of the individual extract_array outputs are effectively column-major.
        return false;
    }

    bool uses_oracle(bool) const {
        // TODO: add support for oracular extraction.
        return false;
    }

private:
    UnknownMatrixCore<Data, Index> core;

private:
    template<bool byrow, tatami::DimensionSelectionType selection_, bool sparse_>
    struct UnknownExtractor : public tatami::Extractor<selection_, sparse_, Data, Index> {
        template<typename ... Args_>
        static auto setup_workspace(Args_&&... args) {
            typename UnknownMatrixCore<Data, Index>::template Workspace<sparse_>* tmp;

#ifdef RATICATE_PARALLELIZE_UNKNOWN 
            // This involves some Rcpp initializations, so we lock it just in case.
            auto& par = parallel_coordinator();
            auto& ex = unknown_evaluator<Data, Index>();
            par.template lock<Data, Index>(
                [&]() -> void {
                    ex.set(&tmp, std::forward<Args_>(args)...);
                },
                [&]() -> void {
                    tmp = new typename UnknownMatrixCore<Data, Index>::template Workspace<sparse_>(std::forward<Args_>(args)...);
                }
            );
#else
            tmp = new typename UnknownMatrixCore<Data, Index>::template Workspace<sparse_>(std::forward<Args_>(args)...);
#endif

            return tmp;
        }

        UnknownExtractor(const UnknownMatrixCore<Data, Index>* c) : core(c) { 
            if constexpr(selection_ == tatami::DimensionSelectionType::FULL) {
                this->full_length = byrow ? core->ncol : core->nrow;
                work.reset(setup_workspace(this->full_length));
            }
        }

        UnknownExtractor(const UnknownMatrixCore<Data, Index>* c, Index start, Index length) : core(c) {
            if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                this->block_start = start;
                this->block_length = length;
                work.reset(setup_workspace(start, length));
            }
        }

        UnknownExtractor(const UnknownMatrixCore<Data, Index>* c, std::vector<Index> idx) : core(c) {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                indices = std::move(idx);
                this->index_length = indices.size();
                work.reset(setup_workspace(indices));
            }
        }

        const Index* index_start() const {
            if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                return indices.data();
            } else {
                return NULL;
            }
        }

        void set_oracle(std::unique_ptr<tatami::Oracle<Index> >) {
            // TODO: support this.
            return;
        }

    protected:
        const UnknownMatrixCore<Data, Index>* core;
        std::unique_ptr<typename UnknownMatrixCore<Data, Index>::template Workspace<sparse_> > work;
        typename std::conditional<selection_ == tatami::DimensionSelectionType::INDEX, std::vector<Index>, bool>::type indices;
    };

private:
    template<bool byrow, tatami::DimensionSelectionType selection_>
    struct DenseUnknownExtractor : public UnknownExtractor<byrow, selection_, false> {
        template<typename ... Args_>
        DenseUnknownExtractor(const UnknownMatrixCore<Data, Index>* c, tatami::Options opt, Args_... args) : 
            UnknownExtractor<byrow, selection_, false>(c, std::forward<Args_>(args)...), options(std::move(opt)) {}

        const Data* fetch(Index i, Data* buffer) {
            if (this->core->needs_reset(i, this->work.get())) {
#ifndef RATICATE_PARALLELIZE_UNKNOWN 
                this->core->template run_dense_extractor<byrow>(i, options, this->work.get());
#else
                auto& ex = unknown_evaluator<Data, Index>();
                auto& par = parallel_coordinator();
                par.template lock<Data, Index>(
                    [&]() -> void {
                        ex.template set<byrow>(i, options, this->work.get(), this->core);
                    },
                    [&]() -> void {
                        this->core->template run_dense_extractor<byrow>(i, options, this->work.get());
                    }
                );
#endif
            }

            i -= this->work->primary_block_start;
            return this->work->bufextractor->fetch_copy(i, buffer); // making a copy to avoid a possible reference to a transient buffer.
        }

    private:
        tatami::Options options;
    };

private:
    template<bool byrow, tatami::DimensionSelectionType selection_>
    struct SparseUnknownExtractor : public UnknownExtractor<byrow, selection_, true> {
        template<typename ... Args_>
        SparseUnknownExtractor(const UnknownMatrixCore<Data, Index>* c, tatami::Options opt, Args_... args) : 
            UnknownExtractor<byrow, selection_, true>(c, std::forward<Args_>(args)...), options(std::move(opt)) {}

        tatami::SparseRange<Data, Index> fetch(Index i, Data* vbuffer, Index* ibuffer) {
            if (this->core->needs_reset(i, this->work.get())) {
#ifndef RATICATE_PARALLELIZE_UNKNOWN
                this->core->template run_sparse_extractor<byrow>(i, options, this->work.get());
#else
                auto& ex = unknown_evaluator<Data, Index>();
                auto& par = parallel_coordinator();
                par.template lock<Data, Index>(
                    [&]() -> void {
                        ex.template set<byrow>(i, options, this->work.get(), this->core);
                    },
                    [&]() -> void {
                        this->core->template run_sparse_extractor<byrow>(i, options, this->work.get());
                    }
                );
#endif
            }

            i -= this->work->primary_block_start;
            tatami::SparseRange<Data, Index> output = this->work->bufextractor->fetch_copy(i, vbuffer, ibuffer);

            // Need to adjust the indices.
            if (output.index) {
                if constexpr(selection_ == tatami::DimensionSelectionType::BLOCK) {
                    for (size_t i = 0; i < output.number; ++i) {
                        ibuffer[i] = output.index[i] + this->block_start;
                    }
                    output.index = ibuffer;
                } else if constexpr(selection_ == tatami::DimensionSelectionType::INDEX) {
                    for (size_t i = 0; i < output.number; ++i) {
                        ibuffer[i] = this->indices[output.index[i]];
                    }
                    output.index = ibuffer;
                }
            }

            return output;
        }

    private:
        tatami::Options options;
    };

public:
    std::unique_ptr<tatami::FullDenseExtractor<Data, Index> > dense_row(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullDenseExtractor<Data, Index> >(new DenseUnknownExtractor<true, tatami::DimensionSelectionType::FULL>(&core, opt));
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Data, Index> > dense_row(Index block_start, Index block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockDenseExtractor<Data, Index> >(new DenseUnknownExtractor<true, tatami::DimensionSelectionType::BLOCK>(&core, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Data, Index> > dense_row(std::vector<Index> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexDenseExtractor<Data, Index> >(new DenseUnknownExtractor<true, tatami::DimensionSelectionType::INDEX>(&core, opt, std::move(indices)));
    }

    std::unique_ptr<tatami::FullDenseExtractor<Data, Index> > dense_column(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullDenseExtractor<Data, Index> >(new DenseUnknownExtractor<false, tatami::DimensionSelectionType::FULL>(&core, opt));
    }

    std::unique_ptr<tatami::BlockDenseExtractor<Data, Index> > dense_column(Index block_start, Index block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockDenseExtractor<Data, Index> >(new DenseUnknownExtractor<false, tatami::DimensionSelectionType::BLOCK>(&core, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexDenseExtractor<Data, Index> > dense_column(std::vector<Index> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexDenseExtractor<Data, Index> >(new DenseUnknownExtractor<false, tatami::DimensionSelectionType::INDEX>(&core, opt, std::move(indices)));
    }

public:
    std::unique_ptr<tatami::FullSparseExtractor<Data, Index> > sparse_row(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullSparseExtractor<Data, Index> >(new SparseUnknownExtractor<true, tatami::DimensionSelectionType::FULL>(&core, opt));
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Data, Index> > sparse_row(Index block_start, Index block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockSparseExtractor<Data, Index> >(new SparseUnknownExtractor<true, tatami::DimensionSelectionType::BLOCK>(&core, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Data, Index> > sparse_row(std::vector<Index> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexSparseExtractor<Data, Index> >(new SparseUnknownExtractor<true, tatami::DimensionSelectionType::INDEX>(&core, opt, std::move(indices)));
    }

    std::unique_ptr<tatami::FullSparseExtractor<Data, Index> > sparse_column(const tatami::Options& opt) const {
        return std::unique_ptr<tatami::FullSparseExtractor<Data, Index> >(new SparseUnknownExtractor<false, tatami::DimensionSelectionType::FULL>(&core, opt));
    }

    std::unique_ptr<tatami::BlockSparseExtractor<Data, Index> > sparse_column(Index block_start, Index block_length, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::BlockSparseExtractor<Data, Index> >(new SparseUnknownExtractor<false, tatami::DimensionSelectionType::BLOCK>(&core, opt, block_start, block_length));
    }

    std::unique_ptr<tatami::IndexSparseExtractor<Data, Index> > sparse_column(std::vector<Index> indices, const tatami::Options& opt) const {
        return std::unique_ptr<tatami::IndexSparseExtractor<Data, Index> >(new SparseUnknownExtractor<false, tatami::DimensionSelectionType::INDEX>(&core, opt, std::move(indices)));
    }
};

}

#endif
