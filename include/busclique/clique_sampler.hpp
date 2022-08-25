#pragma once
#include<iostream>
#include<ostream>
#include<iomanip>

#include<cstdint>
#include<vector>
#include "util.hpp"
#include "cell_cache.hpp"
#include "bundle_cache.hpp"
#include "clique_cache.hpp"

using std::vector;
using std::min;

namespace busclique{


//! This is an arbitrary-precision unsigned integer class with very very little
//! functionality.  Externally, we support three operations: 
//!     * accumulation
//!     * trial subtraction
//!     * generation of random numbers modulo a bigbadint
//! this functionality is the minimum necessary to provide uniform sampling over
//! the optimal cliques within a clique cache.
//! We use Lemire's "nearly-divisionless" algorithm for random number generation
//! with a slight twist to completely avoid division in the multiple-limb case!
class bigbadint {
    friend std::ostream &operator<<(std::ostream &, const bigbadint &);
    friend class uniform_bigint_distribution;
  public:
    typedef uint64_t word;
    typedef uint32_t half;
    static constexpr auto halfbits = 8*sizeof(half);
    static constexpr auto wordbits = 8*sizeof(word);
    static constexpr auto halfmask = ~half(0);
    static constexpr auto highmask = word(halfmask) << halfbits;
  private:
    vector<word> data;

    static size_t clz(uint64_t x) {
        static_assert(wordbits == 64, "bigbadint needs 64-bit words or a new clz method");
        if (x == 0) return(64);
        size_t n = 0;
        if (x <= 0x00000000FFFFFFFF) {n = n +32; x = x <<32;}
        if (x <= 0x0000FFFFFFFFFFFF) {n = n +16; x = x <<16;}
        if (x <= 0x00FFFFFFFFFFFFFF) {n = n + 8; x = x << 8;}
        if (x <= 0x0FFFFFFFFFFFFFFF) {n = n + 4; x = x << 4;}
        if (x <= 0x3FFFFFFFFFFFFFFF) {n = n + 2; x = x << 2;}
        if (x <= 0x7FFFFFFFFFFFFFFF) {n = n + 1;}
        return n;
    }
    
  public:
    size_t clz() { return clz(size()?data.back():0); }

  private:
    static bool add(word &a, word b) {
        bool carry = (a+b) < a;
        a += b;
        return carry;        
    }

    static bool add(word &a, word b, bool c) {
        bool carry = (a+b+c) < a;
        a += b+c;
        return carry;
    }

  public:
    class zero_tag {};
    bigbadint(zero_tag, size_t s) : data(s, 0) {}
    bigbadint() {}
    bigbadint(word x) : data((x>0), x) {}
    bigbadint(word x, size_t reserve) {
        data.reserve(reserve);
        if (x) data.push_back(x);
    }
    bigbadint(bigbadint x, size_t reserve) {
        data.reserve(reserve);
        std::copy(x.data.begin(), x.data.end(), std::back_inserter(data));
    }
    
    template<typename R>
    bigbadint(R &rng, size_t size) {
        data.reserve(size);
        for(size_t i=size; i--;) {
            data.push_back(rng());
        }
        while (size && !data.back()) {
            data.back() = rng();
        }
    }   

    bigbadint &operator+=(bigbadint other) {
        bool carry = 0;
        size_t size = min(data.size(), other.size());
        for (size_t i = 0; i < size; i++)
            carry = add(data[i], other.data[i], carry);
        for (size_t i = size; i < data.size() && carry; i++)
            carry = add(data[i], carry);
        for (size_t i = size; i < other.size(); i++) {
            data.push_back(other.data[i]);
            carry = add(data[i], carry);
        }
        if (carry)
            data.push_back(carry);
        return *this;
    }

  private:
    bool less_bound(const bigbadint &other, size_t i) const {
        for (;i--;) {
            if (data[i] != other.data[i]) {
                return data[i] < other.data[i];
            }
        }
        return false;
    }

  public:
    bool operator<(const bigbadint &other) const {
        if (size() != other.size()) {
            return size() < other.size();
        }
        return less_bound(other, size());
    }

    size_t size() const {
        return data.size();
    }
    
    void clear() {
        data.clear();
    }

    void reserve(size_t size) {
        data.reserve(size);
    }

  private:
    static word mac(word &w, word u, word v, word c) {
        //maximum: (2^n-1)*(2^n-1) + (2^n-1) = (2^n-1) * 2^n fits in 2 words!
        word x = (u&halfmask)*(v&halfmask);
        word y = (u>>halfbits)*(v&halfmask) + (x>>halfbits);
        word z = (u&halfmask)*(v>>halfbits) + (y&halfmask);
        w = (z<<halfbits) + (x&halfmask);
        return (u>>halfbits)*(v>>halfbits) + (y>>halfbits) + (z>>halfbits) + add(w, c);
    }

    void do_mul(const vector<word> left, const vector<word> right, vector<word> &prod) const {
        for (size_t i = 0; i < left.size(); i++) {
            word carry = 0;
            for (size_t j = 0; j < right.size(); j++)
                carry = mac(prod[i+j], left[i], right[i], carry);
            prod[i + right.size()] = carry;
        }
        while(prod.size() && !prod.back())
            prod.pop_back();
    }

  public:
    void mul(const bigbadint &other, bigbadint &result) const {
        result.data.clear();
        result.data.resize(size() + other.size());
        do_mul(data, other.data, result.data);
    }

    bigbadint operator *(const bigbadint &other) const {
        bigbadint result(zero_tag{}, size() + other.size());
        do_mul(data, other.data, result.data);
        return result;
    }

  private:
    static bool sub(word &a, bool c) {
        bool borrow = (a < c);
        a -= c;
        return borrow;
    }
    static bool sub(word &a, word b, bool c) {
        bool borrow = sub(a, c);
        borrow |= a < b;
        a -= b;
        return borrow;
    }

  public:
    bool trysubtract(const bigbadint &x) {
        if (*this < x) return false;
        bool borrow = 0;
        for (size_t i = 0; i < x.size(); i++)
            borrow = sub(data[i], x.data[i], borrow);
        for (size_t i = x.size(); (i < size()) && borrow; i++)
            borrow = sub(data[i], borrow);
        while (data.size() && !data.back())
            data.pop_back();
        return true;
    }

    template<typename R>
    bigbadint random_mod(R &rng) const {
        if (!size())
            return bigbadint(0);
        //keep the sampling bias under 2^63
        bigbadint x(rng, size()+1);
        bigbadint m = operator*(x);
        m.data.erase(m.data.begin(), m.data.begin()+size()+1);
        return m;
    }
};

std::ostream &operator<<(std::ostream &o, const bigbadint &bbi) {
    auto flags = o.flags();
    auto width = o.width();
    o << std::hex;
    for(size_t i = bbi.size(); i--;) {
        o << bbi.data[i];
        o << std::setw(sizeof(bigbadint::word)*2) << std::setfill('0');
    }
    o.flags(flags);
    o.width(width);
    return o;
}

template<typename topo_spec>
class clique_sampler {
    struct optima {
        optima() : sum(0), c(corner::none), score(0) {}
        bigbadint sum;
        corner c;
        size_t score;
        vector<std::tuple<size_y, size_x, corner, bundle_mask>> prev;
    };

    const topo_spec &topo;
  public:
    const size_t width;
  private:
    vector<vector<vector<optima>>> mem;
  public:
    const size_t yield;
    const bigbadint total;
    bool empty;
  private:
    
    vector<vector<vector<optima>>> init_mem(
        const bundle_cache<topo_spec> &bundles,
        size_t length
    ) {
        vector<vector<vector<optima>>> m;
        for (size_t i = 0; i < width; i++) {
            vector<optima> row(clique_cache<topo_spec>::memcols(topo, width, i));
            vector<vector<optima>> block(clique_cache<topo_spec>::memrows(topo, width, i), row);
            m.push_back(block);
        }
        
        auto spy = [&m](size_t i, size_y y, size_x x, size_y py, size_x px, 
                          size_t score, corner c, bundle_mask b) {
            auto &op = m[i][coordinate_index(y)][coordinate_index(x)];
            if (op.score > score)
                return;
            if (op.score < score) {
                op.score = score;
                op.prev.clear();
                op.sum.clear();
            }
            op.prev.emplace_back(py, px, c, b);
            op.sum += (i==0)?bigbadint(1):m[i-1][coordinate_index(py)][coordinate_index(px)].sum;
        };
        if (length) {
            auto check_length = [&bundles, length](size_y yc, size_x xc,
                                                   size_y y0, size_y y1,
                                                   size_x x0, size_x x1){
                return bundles.length(yc,xc,y0,y1,x0,x1) <= length;
            };
            clique_cache<topo_spec> tmp(bundles.cells, bundles, width, check_length, spy);
        } else {
            clique_cache<topo_spec> tmp(bundles.cells, bundles, width, clique_cache<topo_spec>::nocheck, spy);
        }
        return m;
    }
    size_t init_yield() {
        size_t y = 0;
        for (auto &row: mem.back())
            for (auto &entry: row)
                y = max(entry.score, y);
        return y;
    }
    
    bigbadint init_total() {
        bigbadint t(0);
        for (auto &row: mem.back())
            for (auto &entry: row)
                if (entry.score == yield)
                    t += entry.sum;
        return t;
    }
    
  public:
    clique_sampler(const bundle_cache<pegasus_spec> &bundles, size_t width, size_t length = 0) : 
        topo(bundles.cells.topo),
        width(width),
        mem(init_mem(bundles, length)),
        yield(init_yield()),
        total(init_total()),
        empty(total.size() == 0) {}


    clique_sampler(const bundle_cache<topo_spec> &bundles, size_t width) : 
        topo(bundles.cells.topo),
        width(width),
        mem(init_mem(bundles, 0)),
        yield(init_yield()),
        total(init_total()),
        empty(total.size() == 0) {}

    bool sample(bigbadint &index, vector<vector<size_t>> &emb) const {
        emb.clear();
        size_y y = 0;
        size_x x = 0;
        for (auto &row: mem.back()) {
            x = 0;
            for (auto &entry: row) {
                if (entry.score == yield && !index.trysubtract(entry.sum))
                    goto stop;
                x++;
            }
            y++;
        }
        stop:
        for (size_t i = width; i--;) {
            auto &op = mem[i][y][x];
            for (auto pyxc: op.prev) {
                size_t py = std::get<0>(pyxc);
                size_t px = std::get<1>(pyxc);
                if (!index.trysubtract(i?mem[i-1][py][px].sum:bigbadint(1))) {
                    corner c  = std::get<2>(pyxc);
                    bundle_mask b  = std::get<3>(pyxc);
                    size_y yc, y0, y1;
                    size_x xc, x0, x1;
                    clique_cache<topo_spec>::get_ell_loc(y, x, i, width-i-1, c, yc, y0, y1, xc, x0, x1);
                    if (!(c&corner::skipmask))
                        bundle_cache<topo_spec>::inflate_bundle_mask(topo, yc, xc, y0, y1, x0, x1, b, emb);
                    break;
                }
            }
        }
        return true;
    }

    template<typename R>
    bool sample(R &rng, vector<vector<size_t>> &emb) const {
        if (empty)
            return false;
        bigbadint index = total.random_mod(rng);
        return sample(index, emb);
    }
};

template<typename T>
class const_list {
    typedef struct entry {
        const T item;
        entry *next;
        template<typename ...Args>
        entry(entry *next, Args &&...args) : item(std::forward<Args>(args)...), next(next) {}
    } entry;
    class iterator {
        entry *cur;
      public:
        iterator(entry *cur) : cur(cur) {}
        iterator operator++() { cur = cur->next; return *this; }
        bool operator!=(const iterator &other) { return cur != other.cur; }
        const T &operator*() { return cur->item; }
    };
    entry *head;
  public:
    const_list() : head(nullptr) {}
    ~const_list() {
        if (head != nullptr) {
            decapitate();
            pop();
        }
    }
    template<typename ...Args>
    void emplace(Args &...args) {
        head = new entry(head, args...);
    }
    void decapitate() {
        entry *cur = head->next;
        while (cur != nullptr) {
            entry *next = cur->next;
            delete cur;
            cur = next;
        }
        head->next = nullptr;
    }
    void pop() {
        entry *next = head->next;
        delete head;
        head = next;
    }
    const T &back() { return head->item; }
    iterator begin() { return iterator(head); }
    iterator end() { return iterator(nullptr); }
};

template<typename topo_spec>
class clique_sampler_collection {
    const_list<clique_sampler<topo_spec>> samplers;
    bigbadint total;
    size_t yield;
  public:
    clique_sampler_collection() : samplers(), total(0), yield(0) {}
  
    template<typename ...Args>
    void emplace(Args &&...args) {
        samplers.emplace(std::forward<Args>(args)...);
        auto &back = samplers.back();
        if (back.yield < yield)
            samplers.pop();
        else if (back.yield > yield) {
            samplers.decapitate();
            yield = back.yield;
            total = back.total;
        } else {
            total += back.total;
        }
    }
    template<typename R>
    bool sample(R &rng, vector<vector<size_t>> &emb) {
        bigbadint index = total.random_mod(rng);
        for (auto &sampler: samplers) {
            if (!index.trysubtract(sampler.total))
                return sampler.sample(index, emb);
        }
        return false;
    }
};

}
