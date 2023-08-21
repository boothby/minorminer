// Copyright 2020 D-Wave Systems Inc.
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#pragma once
#include "util.hpp"

namespace busclique {

class maxcache {
  public:
    const size_y rows;
    const size_x cols;
  private:
    size_t *mem;
    size_t index(size_y y, size_x x) const {
        return coordinate_converter::grid_index(y, x, rows, cols);
    }
  public:
    maxcache(size_t r, size_t c, size_t *m) : rows(r), cols(c), mem(m) {};
    void setmax(size_y y, size_x x, size_t s, corner c) { 
        size_t old_s = score(y, x);
        if(s == old_s) {
            mem[index(y, x)] |= c;
        } else if(s > old_s) {
            mem[index(y, x)] = (s << corner::shift) | c;
        }
    }
    size_t score(size_y y, size_x x) const {
        return mem[index(y, x)] >> corner::shift;
    }
    corner corners(size_y y, size_x x) const {
        return static_cast<corner>(mem[index(y, x)] & corner::mask);
    }
};

class zerocache {
  public:
    inline constexpr size_t score(size_y, size_x) const { return 0; }
};


template<typename topo_spec> class clique_iterator;

template<typename topo_spec>
class clique_cache {
    friend class clique_iterator<topo_spec>;
  public:
    //prevent double-frees by forbidding moving & copying
    clique_cache(const clique_cache&) = delete; 
    clique_cache(clique_cache &&) = delete;
    
  private:
    const cell_cache<topo_spec> &cells;
    const bundle_cache<topo_spec> &bundles;
    const size_t width;
    const size_t bip_v, bip_h;
    size_t *mem;

    size_t memrows(size_t i) const {
        if (i < width) 
            if (coordinate_index(cells.topo.dim_y) >= i) 
                return coordinate_index(cells.topo.dim_y)-i;
            else
                return 0;
        else if (i == width) return 1;
        else throw "memrows";
    }

    size_t memcols(size_t i) const {
        if (i + 1 < width)
            if (coordinate_index(cells.topo.dim_x)+i+2 > width)
                return coordinate_index(cells.topo.dim_x)-width+i+2;
            else
                return 0;
        else if (i + 1 == width) return coordinate_index(cells.topo.dim_x);
        else throw "memcols";
    }

    size_t memsize(size_t i) const {
        return memrows(i)*memcols(i);
    }

    size_t memsize() const {
        size_t size = 0;
        for(size_t i = 0; i < width; i++)
            size += memsize(i) + 1;
        return size;
    }
    static constexpr bool nocheck(size_y,size_x,size_y,size_y,size_x,size_x) {return true;}

  public:
    clique_cache(const cell_cache<topo_spec> &c, const bundle_cache<topo_spec> &b, size_t w) :
        clique_cache(c, b, w, nocheck) {}

    template<typename C>
    clique_cache(const cell_cache<topo_spec> &c, const bundle_cache<topo_spec> &b, size_t w, C &check) : 
        clique_cache(c, b, w, 0, 0, check) {}


    clique_cache(const cell_cache<topo_spec> &c, const bundle_cache<topo_spec> &b, size_t w, size_t bv, size_t bh) :
        clique_cache(c, b, w, bv, bh, nocheck) {}

    template<typename C>
    clique_cache(const cell_cache<topo_spec> &c, const bundle_cache<topo_spec> &b, size_t w, size_t bv, size_t bh, C &check) : 
            cells(c),
            bundles(b),
            width(w),
            bip_v(bv),
            bip_h(bh),
            mem(new size_t[memsize()]{}) {
        minorminer_assert(size_y(width) <= cells.topo.dim_y);
        minorminer_assert(size_x(width) <= cells.topo.dim_x);
        mem[0] = width;
        for(size_t i = 1; i < width; i++)
            mem[i] = mem[i-1] + memsize(i-1);
        compute_cache(check);
    }

    ~clique_cache() {
        if (mem != nullptr) {
            delete [] mem;
            mem = nullptr;
        }
    }

    maxcache get(size_t i) const {
        minorminer_assert(i < width);
        return maxcache(memrows(i), memcols(i), mem + mem[i]);
    }

    void print() {
        for(size_t i = 0; i < width-1; i++) {
            maxcache m = get(i);
            std::cout << mem[i] << ':' << memsize(i) << "?"<< std::endl;
            for(size_y y = 0; y < m.rows; y++) {
                for(size_x x = 0; x < m.cols; x++) {
                    std::cout << m.score(y, x) << '~' << (m.corners(y, x)) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
  private:
    struct horizontal_part {
        const bundle_cache<topo_spec> &bundles;
        size_y bip_height;
        size_x bip_width;
        horizontal_part(const bundle_cache<topo_spec> &b, size_y h, size_x w) : bundles(b), bip_height(h), bip_width(w) {}
        
        // compute the number of horizontal chains spanning the rectangle
        // [x0, x0+bip_width] \times [y0, y0+bip_height]
        size_t score(size_y y, size_x x0) {
            size_t s = 0;
            for (size_x x = x0; x <= x0 + bip_width; x++)
                s += bundles.get_line_score(1, x, y, y+bip_height);
            return s;
        }
        
        void inflate(size_y y, size_x x, vector<vector<size_t>> &emb) {
            bundles.inflate(1, y, y+bip_height, x, x+bip_width);
        }
    };

    struct vertical_part {
        const bundle_cache<topo_spec> &bundles;
        size_y bip_height;
        size_x bip_width;
        vertical_part(const bundle_cache<topo_spec> &b, size_y h, size_x w) : bundles(b), bip_height(h), bip_width(w) {}

        // compute the number of vertical chains spanning the rectangle
        // [x0, x0+bip_width] \times [y0, y0+bip_height]
        size_t score(size_y y0, size_x x) {
            size_t s = 0;
            for (size_y y = y0; y <= y0 + bip_height; y++)
                s += bundles.get_line_score(0, y, x, x+bip_width);
            return s;
        }

        void inflate(size_y y, size_x x, vector<vector<size_t>> &emb) {
            bundles.inflate(0, y, y+bip_height, x, x+bip_width);
        }
    };
    
    
    template<typename C>
    void compute_cache(C &check) {
        {
            size_y h = 1 + bip_h;
            size_x w = width - bip_h;
            auto zero = zerocache();
            // note to Jose:
            // when bip_h is nonzero, we need to replace this "zerocache" with
            // something that implements a maxcache-like interface, i.e. the
            // horizontal_part struct above (double-check orientation, I'm 
            // constantly getting this wrong).  You'll need to make sure that
            // the computed score is the right-sized rectangle corresponding
            // to the sum of scores computed therein.  Also make sure that 
            // the rectangle spanned by horizontal_part fits tightly with the
            // clique shape.
            // ~~~ busgraph_cache.draw_fragment_embedding is your friend ~~~
            if (bip_h)
                extend_cache(zero, h, w, check, corner::NE, corner::NW, corner::SW, corner::SE);
            else
                extend_cache(zero, h, w, check, corner::SW, corner::SE);            
        }
        for(size_t i = 1 + bip_h; i < width-1-bip_v; i++) {
            size_y h = i+1;
            size_x w = width-i;
            maxcache prev = get(coordinate_index(h-2_y));
            extend_cache(prev, h, w, check, corner::NE, corner::NW, corner::SW, corner::SE);
        }
        {
            size_y h = width-bip_v;
            size_x w = 1+bip_v;
            maxcache prev = get(coordinate_index(h-2_y));
            if (bip_v)
                extend_cache(prev, h, w, check, corner::NE, corner::NW, corner::SW, corner::SE);
            else
                extend_cache(prev, h, w, check, corner::NE, corner::SE);
        }
        // note to Jose: if bip_v is nonzero, it feels like we need one last 
        // cache extension here... but it might be easier to do that in 
        // `extract_solution` at slight cost to performance
    }

    template<typename T, typename C, typename ... Corners>
    inline void extend_cache(const T &prev, size_y h, size_x w, C &check, Corners ... corners) {
        maxcache next = get(coordinate_index(h-1_y));
        for(size_y y = 0; y <= cells.topo.dim_y-h; y++)
            for(size_x x = 0; x <= cells.topo.dim_x-w; x++)
                extend_cache(prev, next, y, y+h-1_y, x, x+w-1_x, check, corners...);
    }

    template<typename T, typename C, typename ... Corners>
    inline void extend_cache(const T &prev, maxcache &next,
                       size_y y0, size_y y1, size_x x0, size_x x1,
                       C &check, corner c, Corners ... corners) {
        extend_cache(prev, next, y0, y1, x0, x1, check, c);
        extend_cache(prev, next, y0, y1, x0, x1, check, corners...);
    }

    template<typename T, typename C>
    inline void extend_cache(const T &prev, maxcache &next,
                       size_y y0, size_y y1, size_x x0, size_x x1,
                       C &check, corner c) {
        size_y next_y, prev_y, yc; next_y = prev_y = yc = y0;
        size_x next_x, prev_x, xc; next_x = prev_x = xc = x0;
        corner skip_c;
        switch(c) {
            case corner::NW: next_x = x0+1_x; prev_y = y0+1_y; skip_c = corner::NWskip; break;
            case corner::SW: next_x = x0+1_x; yc = y1;         skip_c = corner::SWskip; break;
            case corner::NE: xc = x1;         prev_y = y0+1_y; skip_c = corner::NEskip; break;
            case corner::SE: xc = x1;         yc = y1;         skip_c = corner::SEskip; break;
            default: throw std::exception();
        }
        size_t score = prev.score(prev_y, prev_x);
        if(check(yc,xc,y0,y1,x0,x1))
            score += bundles.score(yc,xc,y0,y1,x0,x1);
        else
            c = skip_c;
        next.setmax(next_y, next_x, score, c);
    }

    corner inflate_first_ell(vector<vector<size_t>> &emb,
                             size_y &y, size_x &x, size_y h, size_x w, corner c) const {
        corner c0 = static_cast<corner>(1<< first_bit[c]);
        switch(c0) {
            case corner::NW: x--; bundles.inflate(y,  x,  y,y+h,x,x+w, emb); y++; break;
            case corner::SW: x--; bundles.inflate(y+h,x,  y,y+h,x,x+w, emb);      break;
            case corner::NE:      bundles.inflate(y,  x+w,y,y+h,x,x+w, emb); y++; break;
            case corner::SE:      bundles.inflate(y+h,x+w,y,y+h,x,x+w, emb);      break;
            case corner::NWskip: x--; y++; break;
            case corner::SWskip: x--;      break;
            case corner::NEskip:      y++; break;
            case corner::SEskip:           break;
            default: throw std::exception();
        }
        return c0;
    }

  public:
    bool extract_solution(vector<vector<size_t>> &emb) const {
        minorminer_assert(emb.size() == 0);
        size_y by;
        size_x bx;
        size_t bscore=0;
        maxcache scores = get(width-bip_v-1);
        // note to Jose: make a vertical_part struct here, and find the maximum
        // of `scores.score(y, x) + vpart.score(y, x)` -- again, make sure that
        // there are no off-by-one errors by looking for a tight fit.
        // ~~~ busgraph_cache.draw_fragment_embedding is your friend ~~~
        for(size_y y = 0; y < scores.rows; y++) {
            for(size_x x = 0; x < scores.cols; x++) {
                size_t s = scores.score(y, x);
                if (bscore < s) { 
                    bx = x; by = y; bscore = s;
                }
            }
        }
        if(bscore == 0) return false;
        // inflate the vertical_part here
        // and then, maybe, push an empty chain into the embedding to demark
        // where the vertical chains end
        corner bc = static_cast<corner>(scores.corners(by, bx));
        for(size_t i = width-bip_v-1; i-- > bip_h;) {
            inflate_first_ell(emb, by, bx, i+1, width-2-i, bc);
            bc = static_cast<corner>(get(i).corners(by, bx));
        }
        inflate_first_ell(emb, by, bx, bip_h, width-bip_h-1, bc);
        // again; consider a sentinal empty chain
        // then inflate the horizontal_part here.
        
        return true;
    }
};

template<typename topo_spec>
class clique_iterator {
    const cell_cache<topo_spec> &cells;
    const clique_cache<topo_spec> &cliq;
    size_t width;
    vector<std::tuple<size_y, size_x, corner>> basepoints;
    vector<std::tuple<size_t, size_y, size_x, corner>> stack;
    vector<vector<size_t>> emb;
    
  public:
    clique_iterator(const cell_cache<topo_spec> &c,
                    const clique_cache<topo_spec> &q) :
                    cells(c), cliq(q),
                    width(cliq.width), basepoints(0), stack(0) {

        //to prepare the iterator, we first compute the list of optimal 
        //basepoints -- which consist of a (y, x) pair to denote the rectangle
        //location, and a corner c to denote the orientation of the ell.  
        size_t score = 0;
        maxcache scores = cliq.get(width-1);
        for(size_y y = 0; y < scores.rows; y++)
            for(size_x x = 0; x < scores.cols; x++) {
                size_t s = scores.score(y, x);
                if(s < score) continue;
                else if (s > score) basepoints.clear();
                score = s;
                basepoints.emplace_back(y, x, scores.corners(y, x));
            }
    }

  private:
    bool advance() {
        //first, peel back the zeros (exhausted solutions) until we hit a
        //nonzero corner
        size_t n;
        size_y by;
        size_x bx;
        corner bc;
        while(stack.size()) {
            std::tie(n, by, bx, bc) = stack.back();                
            if(bc) {
                while(emb.size() > n)
                    emb.pop_back();
                break;
            } else {
                stack.pop_back();
            }
        }
        //now, regrow the stack -- if we peeled back the entire stack, we
        //didn't trim the embedding above and we throw the whole thing out
        if(stack.size() == 0) emb.clear();
        while(grow_stack());
        return (emb.size() != 0);
    }

    bool grow_stack() {
        //this grows a stack one step at a time.  the stack is inhomogeneous, so
        //it's a little bothersome to do so -- thus, this isn't a nice recursive
        //function -- we do one step at a time and return true if there's still
        //more work to do; false otherwise.

        size_t n, i;
        size_y by;
        size_x bx;
        corner lc, bc;
        if(stack.size() == 0) {
            //in this case, we've used up the last basepoint that was examined
            //go to the next basepoint (if available) and prime the stack and
            //the current embedding
            if(basepoints.size() == 0)
                return false;
            std::tie(by, bx, bc) = basepoints.back();
            basepoints.pop_back();
            cliq.inflate_first_ell(emb, by, bx, width-1, 0, bc);
            bc = cliq.get(width-2).corners(by, bx);
            //we'll grow the stack from here using the ell located at by, bx
            //with any of the corners in bc 
            stack.emplace_back(emb.size(), by, bx, bc);
            return true;
        } else if((i = stack.size()) < width) {
            //in this case, we're working on a nonempty stack -- so we grow the
            //clique from here
            std::tie(n, by, bx, bc) = stack.back();
            lc = cliq.inflate_first_ell(emb, by, bx, width-i-1, i, bc);
            //here we grow the stack with the first corner in bc -- mark it off
            //as used, so that the next time the stack reaches this depth we
            //won't repeat it
            std::get<3>(stack[i-1]) = static_cast<corner>(lc^bc);
            if(i < width-1) {
                bc = cliq.get(width-2-i).corners(by, bx);
                stack.emplace_back(emb.size(), by, bx, bc);
                return true;
            }
        }
        return false;
    }
  public:
    bool next(vector<vector<size_t>> &e) {
        e.clear();
        if(advance()) {
            for(auto &chain: emb)
                e.emplace_back(chain);
            return true;
        } else return false;
    }
};

const vector<vector<size_t>> empty_emb;

template<typename topo_spec>
class clique_yield_cache {
  private:
    const size_t length_bound;
    vector<size_t> clique_yield;
    vector<vector<vector<size_t>>> best_embeddings;

    size_t compute_length_bound(const zephyr_spec &topo) {
        return 4+topo.zdim;
    }

    size_t compute_length_bound(const pegasus_spec &topo) {
        return 5 + topo.pdim;
    }

    size_t compute_length_bound(const chimera_spec &topo) {
        return 2+coordinate_converter::min(topo.dim_y, topo.dim_x);
    }

  public:
    clique_yield_cache(const cell_cache<topo_spec> &cells) :
                       length_bound(compute_length_bound(cells.topo)),
                       clique_yield(length_bound, 0),
                       best_embeddings(length_bound, empty_emb) { compute_cache(cells); }

  private:
    size_t emb_max_length(const vector<vector<size_t>> &emb) const {
        size_t maxlen = 0;
        for(auto &chain: emb)
            maxlen = max(maxlen, chain.size());
        return maxlen;
    }

    void process_cliques(const clique_cache<topo_spec> &cliques) {
        vector<vector<size_t>> emb;
        if (cliques.extract_solution(emb)) {
            size_t real_len = emb_max_length(emb);
            if(clique_yield[real_len] < emb.size()) {
                clique_yield[real_len] = emb.size();
                best_embeddings[real_len] = emb;
            }
        }
    }

    void compute_cache_width_1(const cell_cache<topo_spec> &cells,
                               const bundle_cache<topo_spec> &bundles) {
        for(size_y y = 0; y < cells.topo.dim_y; y++)
            for(size_x x = 0; x < cells.topo.dim_x; x++) {
                size_t score = bundles.score(y,x,y,y,x,x);
                if(score > clique_yield[2]) {
                    vector<vector<size_t>> emb;
                    bundles.inflate(y,x,y,y,x,x, emb);
                    clique_yield[2] = score;
                    best_embeddings[2] = emb;
                    minorminer_assert(emb.size() == score);
                    minorminer_assert(emb_max_length(emb) == 2);
                }
                if(score == cells.topo.shore) return;
            }
    }

    void compute_cache_width_gt_1(const cell_cache<pegasus_spec> &cells,
                                  const bundle_cache<pegasus_spec> &bundles) {
        size_t maxw = coordinate_converter::min(cells.topo.dim_y, cells.topo.dim_x);

        for(size_t w = 2; w <= maxw; w++) {
            size_t min_length, max_length;
            get_length_range(cells.topo, w, min_length, max_length);
            {
                clique_cache<pegasus_spec> cliques(cells, bundles, w);
                process_cliques(cliques);
            }
            for(size_t len = min_length; len < max_length; len++) {
                auto check_length = [&bundles, len](size_y yc, size_x xc,
                                                    size_y y0, size_y y1,
                                                    size_x x0, size_x x1){
                    return bundles.length(yc,xc,y0,y1,x0,x1) <= len; 
                };
                clique_cache<pegasus_spec> cliques(cells, bundles, w, check_length);
                process_cliques(cliques);
            }
        }
    }
    
    void compute_cache_width_gt_1(const cell_cache<chimera_spec> &cells,
                                  const bundle_cache<chimera_spec> &bundles) {
        size_t maxw = coordinate_converter::min(cells.topo.dim_y, cells.topo.dim_x);

        for(size_t w = 2; w <= maxw; w++) {
            clique_cache<topo_spec> cliques(cells, bundles, w);
            process_cliques(cliques);
        }

    }

    void compute_cache_width_gt_1(const cell_cache<zephyr_spec> &cells,
                                  const bundle_cache<zephyr_spec> &bundles) {
        size_t maxw = coordinate_converter::min(cells.topo.dim_y, cells.topo.dim_x);

        for(size_t w = 2; w <= maxw; w++) {
            clique_cache<topo_spec> cliques(cells, bundles, w);
            process_cliques(cliques);
        }

    }

    void compute_cache(const cell_cache<zephyr_spec> &cells) {
        bundle_cache<zephyr_spec> bundles(cells);
        compute_cache_width_1(cells, bundles);
        compute_cache_width_gt_1(cells, bundles);
    }

    void compute_cache(const cell_cache<chimera_spec> &cells) {
        bundle_cache<chimera_spec> bundles(cells);
        compute_cache_width_1(cells, bundles);
        compute_cache_width_gt_1(cells, bundles);
    }
    
    void compute_cache(const cell_cache<pegasus_spec> &cells) {
        bundle_cache<pegasus_spec> bundles(cells);
        compute_cache_width_gt_1(cells, bundles);
    }

  public:
    const vector<vector<vector<size_t>>> &embeddings() {
        return best_embeddings;
    }
  
    static void get_length_range(const pegasus_spec &topo, size_t width, size_t &min_length, size_t &max_length) {
        max_length = 0;
        min_length = ~max_length;
        for(size_t i = 0; i < width; i++) {
            size_y h = width-i;
            size_x w = i+1;
            for(size_y y = 6; y < size_y(12); y++) {
                for(size_x x = 6; x < size_x(12); x++) {
                    size_t length;
                    length = topo.line_length(0, vert(x), vert(y), vert(y+h-1_y)) + topo.line_length(1, horz(y), horz(x), horz(x+w-1_x));
                    max_length = max(max_length, length);
                    min_length = min(min_length, length);
                    length = topo.line_length(0, vert(x), vert(y), vert(y+h-1_y)) + topo.line_length(1, horz(y+h-1_y), horz(x), horz(x+w-1_x));
                    max_length = max(max_length, length);
                    min_length = min(min_length, length);
                    length = topo.line_length(0, vert(x+w-1_x), vert(y), vert(y+h-1_y)) + topo.line_length(1, horz(y), horz(x), horz(x+w-1_x));
                    max_length = max(max_length, length);
                    min_length = min(min_length, length);
                    length = topo.line_length(0, vert(x+w-1_x), vert(y), vert(y+h-1_y)) + topo.line_length(1, horz(y+h-1_y), horz(x), horz(x+w-1_x));
                    max_length = max(max_length, length);
                    min_length = min(min_length, length);
                }
            }
        }
    }
};


}
