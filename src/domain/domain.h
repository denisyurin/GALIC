#ifndef ALLVARS_H
#include "../allvars.h"
#endif
#ifndef DOMAIN_H
#define DOMAIN_H


extern struct local_topnode_data
{
  peanokey Size;		/*!< number of Peano-Hilbert mesh-cells represented by top-level node */
  peanokey StartKey;		/*!< first Peano-Hilbert key in top-level node */
  long long Count;		/*!< counts the number of particles in this top-level node */
  int Daughter;			/*!< index of first daughter cell (out of 8) of top-level node */
  int Leaf;			/*!< if the node is a leaf, this gives its number when all leaves are traversed in Peano-Hilbert order */
  int Parent;
  int PIndex;			/*!< first particle in node */
}
 *topNodes, *branchNodes;	/*!< points to the root node of the top-level tree */

 struct domain_count_data
 {
   int task;
   int count;
   int origintask;
 };



extern struct domain_peano_hilbert_data
 {
   peanokey key;
   int index;
 }
  *mp;



extern struct trans_data
{
  MyIDType ID;
  int new_task;
  int new_index;
  int wrapped;
}
 *trans_table;

extern int N_trans;

extern int Nbranch;

extern double fac_load;


extern double totpartcount;

extern struct domain_cost_data
{
  int   no;
  int   Count;        /*!< a table that gives the total number of particles held by each processor */
}
 *DomainLeaveNode;



/*! toGo[partner] gives the number of particles on the current task that have to go to task 'partner'
 */
extern int *toGo;
extern int *toGet;
extern int *list_NumPart;
extern int *list_load;




int domain_check_for_local_refine_new(int i, MPI_Comm current_comm);
int domain_double_to_int(double d);
double domain_grav_tot_costfactor(int i);
double domain_hydro_tot_costfactor(int i);
void domain_init_sum_cost(void);
void domain_printf(char *buf);
void domain_report_balance(void);
int domain_sort_load(const void *a, const void *b);
int domain_compare_count(const void *a, const void *b);
int domain_sort_task(const void *a, const void *b);
void domain_post_checks(void);
void domain_prechecks(void);
void domain_insertnode(struct local_topnode_data *treeA, struct local_topnode_data *treeB, int noA, int noB);
void domain_add_cost(struct local_topnode_data *treeA, int noA, long long count, double cost, double sphcost);
int domain_compare_count(const void *a, const void *b);
void domain_rearrange_particle_sequence(void);
void domain_combine_topleaves_to_domains(int ncpu, int ndomain);
void domain_findSplit_load_balanced(int ncpu, int ndomain);
int domain_sort_loadorigin(const void *a, const void *b);
int domain_sort_segments(const void *a, const void *b);
void domain_combine_multipledomains(void);
void domain_allocate(void);
void domain_Decomposition(void);
int domain_check_memory_bound(void);
int domain_compare_key(const void *a, const void *b);
int domain_compare_key(const void *a, const void *b);
int domain_compare_toplist(const void *a, const void *b);
double domain_particle_costfactor(int i);
int domain_countToGo(void);
int domain_decompose(void);
int domain_determineTopTree(void);
void domain_exchange(void);
void domain_findExchangeNumbers(int task, int partner, int sphflag, int *send, int *recv);
void domain_findExtent(void);
void domain_findSplit(int cpustart, int ncpu, int first, int last);
void domain_findSplit_balanced(int cpustart, int ncpu, int first, int last);
void domain_free(void);
void domain_shiftSplit(void);
void domain_sumCost(void);
int domain_topsplit(int node, peanokey startkey);
int domain_topsplit_local(int node, peanokey startkey, int mode);
int domain_topsplit_special(void);
int domain_compare_key(const void *a, const void *b);
int domain_check_for_local_refine(int i, MPI_Comm comm, double work);
void domain_free_trick(void);
void domain_allocate_trick(void);
int domain_recursively_combine_topTree(int start, int ncpu);
void domain_walktoptree(int no);
void domain_optimize_domain_to_task_mapping(void);
int domain_compare_count(const void *a, const void *b);
void domain_allocate_lists(void);
void domain_free_lists(void);
void domain_pack_tree_branch(int no, int parent);
int domain_unpack_tree_branch(int no, int parent);
int domain_check_for_local_refine_alt(int i, int *current_taskset);
int domain_reduce_error_flag(int flag, int *current_taskset);
int domain_do_local_refine(int n, int **list);
void domain_preserve_relevant_topnode_data(void);
void domain_find_total_cost(void);
void domain_voronoi_dynamic_update_execute(void);
void domain_prepare_voronoi_dynamic_update(void);
void domain_voronoi_dynamic_flag_particles(void);
void domain_mark_in_trans_table(int i, int task);
void domain_exchange_and_update_DC(void);
int domain_compare_connection_ID(const void *a, const void *b);
int domain_compare_local_trans_data_ID(const void *a, const void *b);
int domain_compare_recv_trans_data_ID(const void *a, const void *b);
int domain_compare_recv_trans_data_oldtask(const void *a, const void *b);

void mysort_domain(void *b, size_t n, size_t s);

#endif
