#ifndef PROTO_H
#define PROTO_H

#include "allvars.h"
#include "forcetree/forcetree.h"

#include <math.h>
#include <stdlib.h>

#ifdef HAVE_HDF5
#include <hdf5.h>
#endif

int cmp_P_Rnd(const void *a, const void *b);
void shuffle_energies(int iter);
double parallel_sort(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *));
double parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar) (const void *, const void *), MPI_Comm comm);
void smooth_stack(double *data, int maxlevel);
double calc_stack_difference(double *d1, double *d2, int l, int i, int j, int maxlevel, double *ref1, double *ref2, double thresh, double *dist, int flag);

#ifdef VER_1_1
double calc_stack_difference_mod(double *d1, double *d2, int l, int i, int j, int maxlevel, double *ref1, double *ref2, double thresh, double *dist, int flag);
double calc_stack_sum(	double *ref, double *thr, int l, int i, int j, int maxlevel, double thresh, double *dist ); 
#endif


double calc_stack_difference_used(double *d1, double *d2, int l, int i, int j, int maxlevel, 
				  double *ref1, double *ref2, double *used1, double *used2,
				  double thresh, int flag);

double eval_smoothed_stack(double *din, int l, int i, int j, int maxlevel, double *ref, double thresh);
void calc_smoothed_stack(double *din, double *dout, int maxlevel, double *ref, double thresh);

double integrate_axisymmetric_jeans(double zstart, double zend, double R, int type);

double h_factor(double R, double z, int type);
double get_beta_of_type(double *pos, int type);

void free_allocated_memory(void);
void force_test(void);
void forcegrid_get_cell(double *pos, int *iR, int *iz, double *fR, double *fz);

double halo_get_potential(double *pos);
void halo_get_acceleration(double *pos, double *acc);
void halo_get_fresh_coordinate(double *pos);
double halo_generate_v(double rad);
double halo_get_potential_from_radius(double r);
double halo_get_density(double *pos);
double halo_get_mass_inside_radius(double r);
double halo_get_escape_speed(double *pos);
double halo_get_sigma2(double *pos);

void disk_get_fresh_coordinate(double *pos);
double disk_get_density(double *pos);
double disk_get_mass_inside_radius(double R);

double bugle_get_mass_inside_radius(double r);
void bulge_get_fresh_coordinate(double *pos);
double bulge_get_density(double *pos);
double bulge_get_mass_inside_radius(double r);
double bulge_get_escape_speed(double *pos);
double bulge_get_potential(double *pos);
double bulge_get_potential_from_radius(double r);
void bulge_get_acceleration(double *pos, double *acc);
double bulge_get_escape_speed(double *pos);
void output_rotcurve(void);

void densitygrid_sample_targetresponse(void);
void enable_core_dumps_and_fpu_exceptions(void);

double h_over_R(double R, double z, int type);


void line_search(void);
void calc_energy_grid_mass_maps(void);
void energygrid_get_cell(double *pos, int *iR, int *iz, double *fR, double *fz);
void calc_disp_components_for_particle(int n, double *v, double *vr2, double *vt2, double *vp2, double *vq2);

void structure_determination(void);
double structure_disk_angmomentum(void);
double structure_gc(double c);

double eval_fit(int n, double *vel, double *newdens, double *olddens);
#ifdef VER_1_1
double eval_fit_mod(int n, double *vel, double *newdens, double *olddens, double *egyROrbit_new, double *egyROrbit_old, 
																								  double *egyTOrbit_new, double *egyTOrbit_old,
																								  double *egyQOrbit_new, double *egyQOrbit_old,
																								  double *egyPOrbit_new, double *egyPOrbit_old );
#endif

double goldensection_search(int n, double ekin_a, double ekin_b, double ekin_c, double f_a, double f_b, double f_c, double *dir, double *egy, double *fnew, int *count);
double eval_fit_anisotropy(int, double alpha, double v, double *rad, double *perp);
void optimize(int n);
//void optimize_std(int n);
void free_all_response_fields(void);
void calc_all_response_fields(void);
void optimize_some_particles(void);

void forcegrid_allocate(void);
double forcegrid_get_potential(double *pos);
void forcegrid_get_acceleration(double *pos, double *acc);
double forcegrid_get_escape_speed(double *pos);

void forcedensitygrid_create(void);
void forcedensitygrid_calculate(void);

void densitygrid_allocate(void);
void densitygrid_get_cell(double *pos, int *iR, int *iz, double *fR, double *fz);
void forcedensitygrid_load(void);
void forcedensitygrid_save(void);

void commit_updates(void);
void init_updates(void);
void calc_global_fit(void);


void energygrid_allocate(void);


void reorient_particle_velocities(int iter);
void update_velocities(int iter);
void initialize_particles(void);

double get_density_of_type(double *pos, int type);
double get_vstream(double *pos, int type);
double get_z_disp_cylindrical(double *pos, int type);
double get_radial_disp_spherical(double *pos, int type);
void get_disp_rtp(double *pos, int type, double *disp_r, double *disp_t, double *disp_p, double *disp_q);
double get_r_disp_tilted(double *pos, int type);
double get_theta_disp_tilted(double *pos, int type);
double get_phi_disp(double *pos, int type);

void calculate_dispfield(void);
void calc_all_response_fields_and_gradients(void);
void log_message(int iter);
void calc_response_dispersion(void);
void allocate_memory(void);
void output_toomre_Q(void);
void add_to_energy_grid(double *pos, double mass, double vr2, double vt2, double vp2, double vq2,
                        double *egyMass, double *egyResponse_r, double *egyResponse_t, double *egyResponse_p, double *egyResponse_q);

double produce_orbit_response_field(double *pos, double *vel, int id, double *mfield, double mass, double timespan, int *orbitstaken);
#ifdef VER_1_1
double produce_orbit_response_field_mod(double *pos, double *vel, int id, double *mfield, double *egyfield_r, double *egyfield_t, double *egyfield_q, double *egyfield_p, double mass, double timespan, int *orbitstaken, int type);
#endif

void init(void);
void set_units(void);
void endrun(void);
void output_compile_time_options(void);
void set_softenings(void);

void read_parameter_file(char *fname);

void mpi_printf(const char *fmt, ...);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
void write_file(char *fname, int writeTask, int lastTask);
void get_dataset_name(enum iofields blocknr, char *buf);
void get_Tab_IO_Label(enum iofields blocknr, char *label);
int blockpresent(enum iofields blocknr, int write);
int get_particles_in_block(enum iofields blocknr, int *typelist);
int get_values_per_blockelement(enum iofields blocknr);
int get_datatype_in_block(enum iofields blocknr);
int get_bytes_per_blockelement(enum iofields blocknr, int mode);
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type);
void output_particles(int iter);
void output_density_field(int iter);
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);


void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file,
                                int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);
int dump_memory_table_buffer(char *p);
void mymalloc_init(void);

int permutation_compare(const void *a, const void *b);

double dabs(double a);
double dmax(double a, double b);
double dmin(double a, double b);
int imax(int a, int b);
int imin(int a, int b);
int get_part_count_this_task(int n);
size_t sizemax(size_t a, size_t b);
int my_ffsll(long long int i);
void reorder_particles(int *Id);
void gravity(void);

double second(void);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);
double timediff(double t0, double t1);

int get_thread_num(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
void peano_hilbert_order(void);
void peano_hilbert_key_inverse(peanokey key, int bits, int *x, int *y, int *z);
double mysort(void *base, size_t nel, size_t width, int (*compar) (const void *, const void *));

#endif

