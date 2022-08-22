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
    typedef uint64_t word;
    typedef uint32_t half;
    static constexpr auto halfbits = 8*sizeof(half);
    static constexpr auto wordbits = 8*sizeof(word);
    static constexpr auto halfmask = ~half(0);
    static constexpr auto highmask = word(halfmask) << halfbits;
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

    //! compute pow(2, 64*size()) mod *this to aid the "nearly-divisionless"
    //! random number algorithm
    bigbadint magic_constant() const {
        if (size() == 1)
            return bigbadint((-data[0])%data[0]);
            

        bigbadint result(zero_tag{}, size());
        auto n = word(1) << (wordbits-clz(data.back())-1);
        result.data.back() = n;
        while (n > 0) {
            n += n;
            result += result;
            result.trysubtract(*this);
        }
        
        //there you have it, a remainder computed in linear time
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
    bigbadint random_mod(R &rng, const bigbadint &magic_const) const {
        //and here's the good stuff!
    
        if (!size())
            return bigbadint(0);
        bigbadint x(rng, size());
        bigbadint m = operator*(x);
        bool loop = m.size() < size() || m.less_bound(*this, size());
        if (loop)
            while (m.size() < magic_const.size() || m.less_bound(magic_const, size())) {
                x.data.clear();
                for (size_t i = size(); i--;)
                    x.data.push_back(rng());
                while(!x.data.back()) x.data.back() = rng();
                mul(x, m);
            }
        for (size_t hi = m.size(), lo = m.size()-size(); hi--, lo--;)
            m.data[lo] = m.data[hi];
        m.data.resize(m.size()-size());
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
        vector<std::tuple<size_y, size_x, corner>> prev;
    };
    vector<vector<vector<optima>>> mem;

    const clique_cache<topo_spec> &cliques;
  public:
    const size_t yield;
  private:
    bigbadint total;
    bool empty;
    bigbadint magic_constant;
    
  public:
    clique_sampler(const clique_cache<topo_spec> &cliques) : 
        mem([&cliques](){
            vector<vector<vector<optima>>> m;
            for (size_t i = 0; i < cliques.width; i++) {
                vector<optima> row(cliques.memcols(i));
                vector<vector<optima>> block(cliques.memrows(i), row);
                m.push_back(block);
            }
            return m;
        }()),
        cliques(cliques),
        yield([&cliques](){
            const auto &mc = cliques.get(cliques.width-1);
            size_t Y = 0;
            for (size_y y = 0; y < mc.rows; y++)
                for (size_x x = 0; x < mc.cols; x++)
                    Y = max(Y, mc.score(y, x));
            return Y;
        }())
    {
        auto spy = [this](size_t i, size_y y, size_x x, size_y py, size_x px, size_t score, corner c) {
            auto &op = mem[i][coordinate_index(y)][coordinate_index(x)];
            if (op.score > score)
                return;
            if (op.score < score) {
                op.score = score;
                op.prev.clear();
                op.sum.clear();
            }
            op.prev.emplace_back(py, px, c);
            op.sum += (i==0)?bigbadint(1):mem[i-1][coordinate_index(py)][coordinate_index(px)].sum;
        };
        clique_cache<topo_spec> tmp(cliques.cells, cliques.bundles, cliques.width, clique_cache<topo_spec>::nocheck, spy);
        for (auto &row: mem.back())
            for (auto &entry: row)
                if (entry.score == yield)
                    total += entry.sum;
        empty = total.size() == 0;
        if (!empty)
            magic_constant = total.magic_constant();
        std::cout << total << std::endl;
    }


    template<typename R>
    void sample(R &rng, vector<vector<size_t>> &emb) {
        emb.clear();
        if (empty)
            return;
        bigbadint index = total.random_mod(rng, magic_constant);
        size_y y = 0;
        size_x x = 0;
        {
            const auto &mc = cliques.get(cliques.width-1);
            for (auto &row: mem.back()) {
                x = 0;
                for (auto &entry: row) {
                    if (mc.score(y, x) == yield && !index.trysubtract(entry.sum))
                        goto stop;
                    x++;
                }
                y++;
            }
        }
        stop:
        for (size_t i = cliques.width; i--;) {
            auto &op = mem[i][y][x];
            for (auto pyxc: op.prev) {
                size_t py = std::get<0>(pyxc);
                size_t px = std::get<1>(pyxc);
                corner c  = std::get<2>(pyxc);
                if (!index.trysubtract(i?mem[i-1][py][px].sum:bigbadint(1))) {
                    cliques.inflate_first_ell(emb, y, x, i, cliques.width-1-i, c);
                    break;
                }
            }
        }
    }
};
}
