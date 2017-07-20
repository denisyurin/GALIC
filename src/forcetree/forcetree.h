#ifndef FORCETREE_H
#define FORCETREE_H

#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif


/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000

#define MAX_TREE_LEVEL        30
#define MAX_TREE_ALLOC_FACTOR 30.0
#define MAX_IMPACT_BEFORE_OPTIMIZATION 1.03


#define BITFLAG_TOPLEVEL                    0
#define BITFLAG_DEPENDS_ON_LOCAL_MASS       1
#define BITFLAG_DEPENDS_ON_EXTERN_MASS      2
#define BITFLAG_INTERNAL_TOPLEVEL           6
#define BITFLAG_MULTIPLEPARTICLES           7
#define BITFLAG_NODEHASBEENKICKED           8
#define BITFLAG_CONTAINS_GAS                10


#define BITFLAG_MASK  ((1<< BITFLAG_CONTAINS_GAS) + (1 << BITFLAG_MULTIPLEPARTICLES))


static inline unsigned long long force_double_to_int(double d)
{
   union { double d; unsigned long long ull; } u;
   u.d=d;
   return (u.ull&0xFFFFFFFFFFFFFllu); 
}

static inline double force_int_to_double(unsigned long long x)
{
   union { double d; unsigned long long ull; } u;
   u.d = 1.0;
   u.ull |= x;
   return u.d;
}

int force_treebuild(int npart, int optimized_domain_mapping);
int force_treebuild_construct(int npart, int optimized_domain_mapping);
int force_treebuild_insert_single_point(int i, unsigned long long *intpos, int th, unsigned char level);
int force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z);
void force_insert_pseudo_particles(void);
#ifndef GPU_TREE
void force_update_node_recursive(int no, int sib, int father, int *last);
#else
int force_update_node_recursive(int no, int sib, int father, int *last, int depth);
#endif
void force_exchange_topleafdata(void);
void force_treeupdate_toplevel(int no, int topnode, int bits, int x, int y, int z);
void force_treeallocate(int maxpart, int maxindex);
void force_treefree(void);
void dump_particles(void);
int force_add_empty_nodes(void);
void force_short_range_init(void);
int force_treeevaluate(int target, int mode, int thread_id);
int force_treeevaluate_shortrange(int target, int mode, int thread_id, int measure_cost_flag);
int force_treeevaluate_ewald_correction(int i, int mode, int thread_id);
int force_treeevaluate_direct(int target, int mode);
void force_assign_cost_values(void);
void force_update_node_recursive_sse(int no, int sib, int father, int *last);
void force_optimize_domain_mapping(void);
double force_get_current_balance(double *impact);
void force_get_global_cost_for_leavenodes(int nexport);



#endif



