cimport numpy as cnp

cdef extern from "lib/ssw.h":
    
    ctypedef struct s_profile:
        pass

    ctypedef struct s_align:
        cnp.uint16_t score1
        cnp.uint16_t score2
        cnp.int32_t ref_begin1
        cnp.int32_t ref_end1
        cnp.int32_t read_begin1
        cnp.int32_t read_end1
        cnp.int32_t ref_end2
        cnp.uint32_t* cigar
        cnp.int32_t cigarLen

    cdef s_profile* ssw_init(const cnp.int8_t* read,
                             const cnp.int32_t readLen,
                             const cnp.int8_t* mat,
                             const cnp.int32_t n,
                             const cnp.int8_t score_size)

    cdef void init_destroy(s_profile* p)

    cdef s_align* ssw_align(const s_profile* prof,
                            const cnp.int8_t* ref,
                            cnp.int32_t refLen,
                            const cnp.uint8_t weight_gapO,
                            const cnp.uint8_t weight_gapE,
                            const cnp.uint8_t flag,
                            const cnp.uint16_t filters,
                            const cnp.int32_t filterd,
                            const cnp.int32_t maskLen)

    cdef void align_destroy(s_align* a)

