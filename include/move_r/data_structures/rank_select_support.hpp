#pragma once

#include <vector>
#include <iostream>
#include <filesystem>
#include <omp.h>
#include <sdsl/wt_gmr.hpp>

#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/hybrid_bit_vector.hpp>
#include <move_r/data_structures/interleaved_vectors.hpp>

/**
 * @brief a rank-select data structure using hybrid bit vectors (either sd_array or plain bit vector), or the gmr data structure
 * @tparam sym_t symbol type type
 * @tparam pos_t position type
 */
template <typename sym_t, typename pos_t = uint32_t>
class rank_select_support {
    protected:

    static_assert(std::is_same_v<pos_t,uint32_t> || std::is_same_v<pos_t,uint64_t>);
    
    static_assert(
        std::is_same_v<sym_t,char> ||
        std::is_same_v<sym_t,uint8_t> ||
        std::is_same_v<sym_t,uint16_t> ||
        std::is_same_v<sym_t,uint32_t> ||
        std::is_same_v<sym_t,uint64_t>
    );

    static constexpr bool str_input = std::is_same_v<sym_t,char>; // true <=> the input is a string
    static constexpr bool byte_alphabet = sizeof(sym_t) == 1; // true <=> use bit vectors
    static constexpr bool int_alphabet = !byte_alphabet; // true <=> use gmr

    using i_sym_t = std::conditional_t<str_input,uint8_t,sym_t>; // internal (unsigned) symbol type
    using hybrid_bv_t = hybrid_bit_vector<pos_t,true,false,true>; // hybrid bit vector type

    // ############################# VARIABLES FOR byte_alphabet = true #############################

    pos_t input_size = 0; // the size of the input
    std::vector<i_sym_t> alphabet; // the symbols occurring in the input sorted descendingly by their frequency
    std::vector<pos_t> freq; // stores at position v the frequency of v in the input

    /**
     * @brief hyb_bit_vecs[v] contains a hybrid bit vector [0..size-1] that marks (with ones) the occurrences of v
     *        in the input, for each value v occurring in the input
     */
    std::vector<hybrid_bv_t> hyb_bit_vecs;

    // ############################# VARIABLES FOR byte_alphabet = false #############################

    pos_t sigma = 0; // the number of distinct symbols in the input

    /**
     * @brief rank-select data structure for integer alphabets
     */
    sdsl::wt_gmr_rs<> gmr;

    /**
     * @brief [0..sigma-1] stores at position v the sum of the frequencies of all v' < v that occur in the input
     */
    interleaved_vectors<pos_t,pos_t> c_arr;
    
    /**
     * @brief [0..input_size-1] stores at position i the position of the j-th occurrence of v
     *        in the input, where c_arr[v] <= i < c_arr[v+1] and j = i - c_arr[v] + 1
     */
    interleaved_vectors<pos_t,pos_t> occs;

    // ######################################################################################

    pos_t symbol_idx(sym_t sym) const {
        if constexpr (str_input) {
            return char_to_uchar(sym);
        } else {
            return sym;
        }
    };

    /**
     * @brief builds the bit vectors
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    void build_bvs(const std::function<sym_t(pos_t)>& read, pos_t l, pos_t r, uint16_t num_threads) {
        input_size = r-l+1;
        std::vector<std::vector<pos_t>> occurrences_thr(num_threads,std::vector<pos_t>(256,0));

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i=l; i<=r; i++) {
            occurrences_thr[omp_get_thread_num()][symbol_idx(read(i))]++;
        }

        freq.resize(256);

        for (uint16_t cur_char=0; cur_char<256; cur_char++) {
            for (uint16_t i_p=0; i_p<num_threads; i_p++) {
                freq[cur_char] += occurrences_thr[i_p][cur_char];
            }

            if (freq[cur_char] != 0) {
                alphabet.emplace_back(cur_char);
            }
        }

        occurrences_thr.clear();
        occurrences_thr.shrink_to_fit();
        sigma = alphabet.size();
        uint16_t num_chars = sigma;
        std::sort(alphabet.begin(),alphabet.end(),[this](i_sym_t a, i_sym_t b){return freq[a] > freq[b];});
        std::vector<uint8_t> is_compressed(256,1);
        std::vector<sdsl::sd_vector_builder> sdv_builders(256);
        std::vector<sdsl::bit_vector> plain_bvs(256);

        #pragma omp parallel for num_threads(num_threads)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = alphabet[i];

            if (freq[cur_char] > input_size * hybrid_bit_vector<pos_t>::compression_threshold) {
                is_compressed[cur_char] = 0;
                plain_bvs[cur_char] = sdsl::bit_vector(input_size);
            } else {
                sdv_builders[cur_char] = sdsl::sd_vector_builder(input_size,freq[cur_char]);
            }
        }

        uint8_t cur_char;

        for (pos_t i=l; i<=r; i++) {
            cur_char = symbol_idx(read(i));

            if (is_compressed[cur_char]) {
                sdv_builders[cur_char].set(i-l);
            } else {
                plain_bvs[cur_char][i-l] = 1;
            }
        }

        hyb_bit_vecs.resize(256);

        #pragma omp parallel for num_threads(num_threads)
        for (uint16_t i=0; i<num_chars; i++) {
            uint8_t cur_char = symbol_idx(alphabet[i]);

            if (is_compressed[cur_char]) {
                hyb_bit_vecs[cur_char] = hybrid_bv_t(sdsl::sd_vector<>(sdv_builders[cur_char]));
            } else {
                hyb_bit_vecs[cur_char] = hybrid_bv_t(std::move(plain_bvs[cur_char]));
            }
        }
    }

    /**
     * @brief builds the gmr data structure
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param alphabet_size number of distinct symbols in vector (= maximum value in the input)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    void build_gmr(const std::function<sym_t(pos_t)>& read, pos_t alphabet_size, pos_t l, pos_t r) {
        sigma = alphabet_size;
        std::string vector_buffer_filename = "vec_" + random_alphanumeric_string(10);
        sdsl::int_vector_buffer<> vector_buffer(vector_buffer_filename,std::ios::out);
        freq.resize(alphabet_size,0);

        for (uint64_t i=l; i<=r; i++) {
            sym_t v = read(i);
            vector_buffer[i-l] = v;
            freq[v]++;
        }

        gmr = sdsl::wt_gmr_rs<>(vector_buffer,r-l+1);
        vector_buffer.close();
        std::filesystem::remove(vector_buffer_filename);
    }

    /**
     * @brief builds the c-array and occurrences
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    void build_select(const std::function<sym_t(pos_t)>& read, pos_t l, pos_t r, uint16_t num_threads) {
        uint8_t bytes_per_entry = (uint8_t)std::ceil(std::log2(r-l+1)/(double)8);
        c_arr = interleaved_vectors<pos_t,pos_t>({bytes_per_entry});
        c_arr.resize_no_init(sigma+1);
        c_arr.template set<0,pos_t>(0,0);

        for (pos_t v=1; v<=sigma; v++) {
            c_arr.template set<0,pos_t>(v,c_arr[v-1]+freq[v-1]);
        }

        freq.clear();
        freq.shrink_to_fit();
        std::vector<pos_t> occ_idx;
        no_init_resize(occ_idx,sigma);

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t v=0; v<sigma; v++) {
            occ_idx[v] = c_arr[v];
        }

        occs = interleaved_vectors<pos_t,pos_t>({bytes_per_entry});
        occs.resize_no_init(r-l+1);

        for (pos_t i=l; i<=r; i++) {
            pos_t v = symbol_idx(read(i));
            occs.template set<0,pos_t>(occ_idx[v],i-l);
            occ_idx[v]++;
        }
    }

    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @param v_s [0..input size-1] starting position of the v-interval
     * @param v_e [1..input size] end position of the v-interval
     * @return number of occurrences of v before index i
     */
    inline pos_t rank_v_interval(
        sym_t v, pos_t i,
        pos_t v_s,
        pos_t v_e
    ) const requires(int_alphabet) {
        pos_t pos;

        if (v_e-v_s <= 256) {
            pos = v_s;

            while (pos < v_e && occs[pos+1] < i) {
                pos++;
            }
        } else {
            pos = bin_search_max_lt<pos_t>(i,v_s,v_e-1,[this](pos_t x){return occs[x];});
        }

        return pos-v_s+1;
    }

    public:
    rank_select_support() = default;

    /**
     * @brief builds the data structure from the range [l,r] in an input (0 <= l <= r < string.size()),
     * else if l > r, then the the data structure is built for the whole input
     * @param string the input
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    rank_select_support(
        const std::string& string,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(byte_alphabet) {
        if (l > r) {
            l = 0;
            r = string.size()-1;
        }
        
        r = std::min<pos_t>(r,string.size()-1);
        build_bvs([&string](uint i){return string[i];},l,r,num_threads);
    }
    
    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    rank_select_support(
        const std::function<sym_t(pos_t)>& read,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(byte_alphabet) {
        build_bvs(read,l,r,num_threads);
    }

    /**
     * @brief builds the data structure by reading the input vector
     * @param vector input vector
     * @param alphabet_size number of distinct symbols in vector (= maximum value in vector)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     * @param num_threads the number of threads to use
     */
    rank_select_support(
        const std::vector<sym_t>& vector,
        pos_t alphabet_size,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(int_alphabet) {
        if (l > r) {
            l = 0;
            r = vector.size()-1;
        }
        
        r = std::min<pos_t>(r,vector.size()-1);
        build_gmr([&vector](uint i){return vector[i];},alphabet_size,l,r);
        build_select([&vector](uint i){return vector[i];},l,r,num_threads);
    }

    /**
     * @brief builds the data structure by reading the input using the function read
     * @param read function to read the input with; it is called with i in [l,r]
     * as a parameter and must return the value of the input at index i
     * @param alphabet_size number of distinct symbols in vector (= maximum value in vector)
     * @param l left range limit (l <= r)
     * @param r right range limit (l <= r)
     */
    rank_select_support(
        const std::function<sym_t(pos_t)>& read,
        pos_t alphabet_size,
        pos_t l = 1, pos_t r = 0,
        uint16_t num_threads = omp_get_max_threads()
    ) requires(int_alphabet) {
        build_gmr(read,alphabet_size,l,r);
        build_select(read,l,r,num_threads);
    }

    /**
     * @brief returns the value at index i in the input
     * @param i [0..size-1]
     * @return the value at index i in the input
     */
    inline sym_t operator[](pos_t i) const {
        if constexpr (byte_alphabet) {
            for (i_sym_t v : alphabet) {
                if (hyb_bit_vecs[v][i] == 1) {
                    return v;
                }
            }

            return 0;
        } else {
            return gmr[i];
        }
    }

    /**
     * @brief returns the size of the input
     * @return the size of the input
     */
    inline pos_t size() const {
        if constexpr (byte_alphabet) {
            return input_size;
        } else {
            return gmr.size();
        }
    }

    /**
     * @brief returns whether the input vector is empty
     * @return whether the input vector is empty
     */
    inline bool empty() const {
        return size() == 0;
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const {
        if constexpr (byte_alphabet) {
            uint64_t size = sizeof(pos_t)+alphabet.size()+sizeof(pos_t)*alphabet.size();

            for (uint16_t i=0; i<hyb_bit_vecs.size(); i++) {
                size += hyb_bit_vecs[i].size_in_bytes();
            }

            return size;
        } else {
            return sdsl::size_in_bytes(gmr)+
                c_arr.size_in_bytes()+
                occs.size_in_bytes();
        }
    }

    /**
     * @brief returns the number of distinct symbols in the input
     * @return the number of distinct symbols in the input
     */
    inline pos_t alphabet_size() const {
        return sigma;
    }

    inline bool contains(sym_t v) const {
        if constexpr (byte_alphabet) {
            return freq[symbol_idx(v)] != 0;
        } else {
            return v < sigma;
        }
    }

    /**
     * @brief returns the number of occurrences of the symbol v in the input (v must occur in the input)
     * @param v symbol that occurs in the input
     * @return the number of occurrences of v in the input
     */
    inline pos_t frequency(sym_t v) const {
        if constexpr (byte_alphabet) {
            return freq[symbol_idx(v)];
        } else {
            return c_arr[v+1]-c_arr[v];
        }
    }

    /**
     * @brief returns the number of occurrences of v before index i
     * @tparam use_gmr true <=> use the gmr data structure
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    template <bool use_gmr>
    inline pos_t rank(sym_t v, pos_t i) const requires(int_alphabet) {
        if constexpr (use_gmr) {
            return gmr.rank(i,v);
        } else {
            pos_t v_s = c_arr[v];
            pos_t v_e = c_arr[v+1]-1;

            if (occs[v_s] >= i) {
                return 0;
            }

            return rank_v_interval(v,i,v_s,v_e);
        }
    }

    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    inline pos_t rank(sym_t v, pos_t i) const {
        if constexpr (byte_alphabet) {
            return hyb_bit_vecs[symbol_idx(v)].rank_1(i);
        } else {
            pos_t v_s = c_arr[v];
            pos_t v_e = c_arr[v+1]-1;

            if (occs[v_s] >= i) {
                return 0;
            }

            if (v_e-v_s > 4096) {
                return gmr.rank(i,v);
            } else {
                return rank_v_interval(v,i,v_s,v_e);
            }
        }
    }
    
    /**
     * @brief returns the index of the i-th occurrence of v
     * @tparam use_gmr true <=> use the gmr data structure
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of v in the input]
     * @return index of the i-th occurrence of v
     */
    template <bool use_gmr>
    inline pos_t select(sym_t v, pos_t i) const requires(int_alphabet) {
        if constexpr (use_gmr) {
            return gmr.select(i,v);
        } else {
            return occs[c_arr[v]+i-1];
        }
    }
    
    /**
     * @brief returns the index of the i-th occurrence of v
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of v in the input]
     * @return index of the i-th occurrence of v
     */
    inline pos_t select(sym_t v, pos_t i) const {
        if constexpr (byte_alphabet) {
            return hyb_bit_vecs[symbol_idx(v)].select_1(i);
        } else {
            return select<false>(v,i);
        }
    }

    /**
     * @brief serializes the rank_select_support to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        if constexpr (byte_alphabet) {
            out.write((char*)&input_size,sizeof(pos_t));
            
            if (input_size > 0) {
                uint8_t num_chars = alphabet.size();
                out.write((char*)&num_chars,1);

                out.write((char*)&alphabet[0],alphabet.size());
                out.write((char*)&freq[0],256*sizeof(pos_t));

                for (uint16_t i=0; i<alphabet.size(); i++) {
                    hyb_bit_vecs[alphabet[i]].serialize(out);
                }
            }
        } else {
            pos_t vector_size = gmr.size();
            out.write((char*)&vector_size,sizeof(pos_t));
            out.write((char*)&sigma,sizeof(pos_t));
            
            if (vector_size > 0) {
                gmr.serialize(out);
                c_arr.serialize(out);
                occs.serialize(out);
            }
        }
    }

    /**
     * @brief loads the rank_select_support from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        if constexpr (byte_alphabet) {
            in.read((char*)&input_size,sizeof(pos_t));

            if (input_size > 0) {
                uint8_t num_chars;
                in.read((char*)&num_chars,1);

                alphabet.resize(num_chars);
                in.read((char*)&alphabet[0],num_chars);

                freq.resize(256);
                in.read((char*)&freq[0],256*sizeof(pos_t));

                hyb_bit_vecs.resize(256);

                for (uint16_t i=0; i<alphabet.size(); i++) {
                    hyb_bit_vecs[alphabet[i]].load(in);
                }
            }
        } else {
            pos_t vector_size;
            in.read((char*)&vector_size,sizeof(pos_t));
            in.read((char*)&sigma,sizeof(pos_t));

            if (vector_size > 0) {
                gmr.load(in);
                c_arr.load(in);
                occs.load(in);
            }
        }
    }

    std::ostream& operator>>(std::ostream& os) const {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is) {
        load(is);
        return is;
    }
};