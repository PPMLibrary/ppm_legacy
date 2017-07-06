#----------------------------------------------------------------------
#								      
#  Makefile : libppm.a     
#  Version  :							  
#  Author   :							  
#  Created  : 26/7 - 2004	  
#								      
#----------------------------------------------------------------------
#  Creator: makemake.pl 1.2
#  Rev: 2003/06/06 by Jens H. Walther  
#----------------------------------------------------------------------
SHELL       = /bin/sh
TARGET      = libppm.a
BIN	 = libppm.a
#----------------------------------------------------------------------
#  Compiler settings (Linux)
#----------------------------------------------------------------------
include Make.ifc
#include Make.pgf90
#include Make.g5
#include Make.mpiifc

#----------------------------------------------------------------------
#  Search path for RCS files					   
#----------------------------------------------------------------------
VPATH = ./RCS

#----------------------------------------------------------------------
#  Additional suffix rules					     
#----------------------------------------------------------------------
.SUFFIXES : .inc .inc,v .f,v .c,v
.f,v.f :
	 co $*.f

.c,v.c :
	 co $*.c

.inc,v.inc :
	 co $*.inc

#----------------------------------------------------------------------
#  Binary directory
#----------------------------------------------------------------------
bindir      = ../lib

#----------------------------------------------------------------------
#  Default target
#----------------------------------------------------------------------
all: $(TARGET)

#----------------------------------------------------------------------
#  Object files:						       
#----------------------------------------------------------------------
OBJS =\
	ppm_module_data.o\
	ppm_module_data_fmm.o\
	ppm_module_data_io.o\
	ppm_module_data_loadbal.o\
	ppm_module_data_mesh.o\
	ppm_module_data_mg.o\
	ppm_module_data_neighlist.o\
	ppm_module_data_part.o\
	ppm_module_data_ode.o\
	ppm_module_data_rmsh.o\
	ppm_module_data_tree.o\
	ppm_module_data_fieldsolver.o\
	ppm_module_data_gmm.o\
	ppm_module_write.o\
	ppm_module_util_time.o\
	ppm_module_substart.o\
	ppm_module_substop.o\
	ppm_module_log.o\
	ppm_module_error.o\
	ppm_module_alloc.o\
	ppm_module_check_topoid.o\
	ppm_module_check_meshid.o\
	ppm_module_util_cart2sph.o\
	ppm_module_util_commopt.o\
	ppm_module_util_fft_backward.o\
	ppm_module_util_fft_forward.o\
	ppm_module_util_invert_list.o\
	ppm_module_util_qsort.o\
	ppm_module_util_collapse_list.o\
	ppm_module_util_rank.o\
	ppm_module_util_sort.o\
	ppm_module_util_matinv_2x2.o\
	ppm_module_util_matinv_3x3.o\
	ppm_module_util_quadeq_real.o\
	ppm_module_util_cubeq_real.o\
	ppm_module_util_eigen_2sym.o\
	ppm_module_util_eigen_3sym.o\
	ppm_module_decomp_boxsplit.o\
	ppm_module_decomp_tree.o\
	ppm_module_tree_done.o\
	ppm_module_tree_alloc.o\
	ppm_module_tree_boxcost.o\
	ppm_module_tree_boxcut.o\
	ppm_module_tree_cutdir.o\
	ppm_module_tree_cutpos.o\
	ppm_module_tree_divcheck.o\
	ppm_module_tree.o\
	ppm_module_bem_basis.o\
	ppm_module_bem_quadrule_npoints.o\
	ppm_module_bem_quadrule.o\
	ppm_module_neighlist_MkNeighIdx.o\
	ppm_module_clist_destroy.o\
	ppm_module_neighlist_clist.o\
	ppm_module_neighlist_vlist.o\
	ppm_module_topo_check.o\
	ppm_module_topo_cost.o\
	ppm_module_topo_inquire.o\
	ppm_module_topo_store.o\
	ppm_module_topo_box2subs.o\
	ppm_module_map_part_global.o\
	ppm_module_map_part_partial.o\
	ppm_module_map_part_pop.o\
	ppm_module_map_part_push.o\
	ppm_module_map_part_remap.o\
	ppm_module_map_part_ring.o\
	ppm_module_map_part_ring_shift.o\
	ppm_module_map_part_send.o\
	ppm_module_map_part_ghost_get.o\
	ppm_module_map_part_get_sub.o\
	ppm_module_map_part_ghost_put.o\
	ppm_module_map_part_ghost.o\
	ppm_module_map_part.o\
	ppm_module_map_connect_distrib.o\
	ppm_module_map_connect_prune.o\
	ppm_module_map_connect_send.o\
	ppm_module_map_connect.o\
	ppm_module_mesh_alloc.o\
	ppm_module_mesh_block_intersect.o\
	ppm_module_map_field_ghost_get.o\
	ppm_module_map_field_ghost_init.o\
	ppm_module_map_field_ghost_put.o\
	ppm_module_map_field_pop_2d.o\
	ppm_module_map_field_pop_3d.o\
	ppm_module_map_field_push_2d.o\
	ppm_module_map_field_push_3d.o\
	ppm_module_map_field_send.o\
	ppm_module_map_field_ghost.o\
	ppm_module_map_field_global.o\
	ppm_module_map_field.o\
	ppm_module_comp_mesh_fd.o\
	ppm_module_comp_mesh_fft.o\
	ppm_module_comp_pp_cell.o\
	ppm_module_comp_pp_correct.o\
	ppm_module_comp_pp_doring.o\
	ppm_module_comp_pp_mk_table.o\
	ppm_module_comp_pp_ring.o\
	ppm_module_comp_pp_verlet.o\
	ppm_module_copy_image_subs.o\
	ppm_module_decomp_cartesian.o\
	ppm_module_decomp_null.o\
	ppm_module_decomp_pruned_cell.o\
	ppm_module_define_subs_bc.o\
	ppm_module_estimate_proc_speed.o\
	ppm_module_find_neigh.o\
	ppm_module_topo_metis_s2p.o\
	ppm_module_topo_subs2proc.o\
	ppm_module_mesh_on_subs.o\
	ppm_module_mesh_store.o\
	ppm_module_mesh_define.o\
	ppm_module_mesh_derive.o\
	ppm_module_mesh_finalize.o\
	ppm_module_mktopo.o\
	ppm_module_find_duplicates.o\
	ppm_module_get_cost.o\
	ppm_module_get_revision.o\
	ppm_module_io_unused_unit.o\
	ppm_module_io_close.o\
	ppm_module_io_delete.o\
	ppm_module_io_inquire.o\
	ppm_module_io_open.o\
	ppm_module_io_read_ascii.o\
	ppm_module_io_read_binary.o\
	ppm_module_io_set_unit.o\
	ppm_module_io_write_ascii.o\
	ppm_module_io_write_binary.o\
	ppm_module_doio.o\
	ppm_module_io.o\
	ppm_module_print_defines.o\
	ppm_module_set_decomp_cost.o\
	ppm_module_set_proc_speed.o\
	ppm_module_is_initialized.o\
	ppm_module_init.o\
	ppm_module_finalize.o\
	ppm_module_loadbal_inquire.o\
	ppm_module_impose_part_bc.o\
	ppm_module_mg_alloc.o\
	ppm_module_ode_alldone.o\
	ppm_module_ode_create_ode.o\
	ppm_module_ode_finalize.o\
	ppm_module_ode_init.o\
	ppm_module_ode_map_pop.o\
	ppm_module_ode_map_push.o\
	ppm_module_ode_start.o\
	ppm_module_ode_step.o\
	ppm_module_mg_finalize.o\
	ppm_module_mg_init.o\
	ppm_module_mg_smooth.o\
	ppm_module_mg_res.o\
	ppm_module_mg_prolong.o\
	ppm_module_mg_restrict.o\
	ppm_module_mg_core.o\
	ppm_module_mg_solv.o\
	ppm_module_rmsh_comp_weights.o\
	ppm_module_rmsh_create_part.o\
	ppm_module_interp_p2m.o\
	ppm_module_interp_m2p.o\
	ppm_module_rmsh_remesh.o\
	ppm_module_time.o\
	ppm_module_scale_domain.o\
	ppm_module_gmm_init.o\
	ppm_module_gmm_finalize.o\
	ppm_module_gmm_kickoff.o\
	ppm_module_gmm_march_fwd.o\
	ppm_module_gmm_march_bkwd.o\
	ppm_module_gmm_extend_fwd.o\
	ppm_module_gmm_extend_bkwd.o\
	ppm_module_gmm_march.o\
	ppm_module_gmm_reinitialize.o\
	ppm_module_gmm_cpt.o\
	ppm_module_gmm_extend.o\
	ppm_module_comp_part.o\
	ppm_module_comp_mesh.o\
	ppm_module_bem.o\
	ppm_module_fdsolver_init.o\
	ppm_module_fdsolver_fft_fd.o\
	ppm_module_fdsolver_fft_bd.o\
	ppm_module_fdsolver_poisson.o\
	ppm_module_fdsolver_map.o\
	ppm_module_fdsolver_solve.o\
	ppm_module_fdsolver_finalize.o\
	ppm_module_fieldsolver.o\
	ppm_module_topo.o\
	ppm_module_map.o\
	ppm_module_fmm_init.o\
	ppm_module_fmm_traverse.o\
	ppm_module_fmm_pretraverse.o\
	ppm_module_fmm_expchange.o\
	ppm_module_fmm_potential.o\
	ppm_module_fmm_finalize.o\
	ppm_module_user_util.o\
	ppm_module_fmm_expansion.o\
	ppm_module_loadbal.o\
	ppm_module_fmm.o\
	ppm_module_user_io.o\
	ppm_module_ode.o\
	ppm_module_mg.o\
	ppm_module_neighlist.o\
	ppm_module_rmsh.o\
	ppm_module_gmm.o\
	ppm_module_user.o

#----------------------------------------------------------------------
#  Dependencies:						       
#----------------------------------------------------------------------
ppm_module_alloc.o: ppm_module_alloc.f\
	./ppm_alloc_1d.f\
	./ppm_alloc_2d.f\
	./ppm_alloc_3d.f\
	./ppm_alloc_4d.f\
	./ppm_alloc_5d.f\
	./ppm_alloc_1dl.f\
	./ppm_alloc_2dl.f\
	./ppm_alloc_3dl.f\
	./ppm_alloc_4dl.f\
	./ppm_alloc_5dl.f\
	./ppm_alloc_1d.f\
	./ppm_alloc_2d.f\
	./ppm_alloc_3d.f\
	./ppm_alloc_4d.f\
	./ppm_alloc_5d.f\
	./ppm_alloc_1dl.f\
	./ppm_alloc_2dl.f\
	./ppm_alloc_3dl.f\
	./ppm_alloc_4dl.f\
	./ppm_alloc_5dl.f\
	./ppm_alloc_1d.f\
	./ppm_alloc_2d.f\
	./ppm_alloc_3d.f\
	./ppm_alloc_4d.f\
	./ppm_alloc_5d.f\
	./ppm_alloc_1dl.f\
	./ppm_alloc_2dl.f\
	./ppm_alloc_3dl.f\
	./ppm_alloc_4dl.f\
	./ppm_alloc_5dl.f\
	./ppm_alloc_1d.f\
	./ppm_alloc_2d.f\
	./ppm_alloc_3d.f\
	./ppm_alloc_4d.f\
	./ppm_alloc_5d.f\
	./ppm_alloc_1dl.f\
	./ppm_alloc_2dl.f\
	./ppm_alloc_3dl.f\
	./ppm_alloc_4dl.f\
	./ppm_alloc_5dl.f\
	./ppm_alloc_1d.f\
	./ppm_alloc_2d.f\
	./ppm_alloc_3d.f\
	./ppm_alloc_4d.f\
	./ppm_alloc_5d.f\
	./ppm_alloc_1dl.f\
	./ppm_alloc_2dl.f\
	./ppm_alloc_3dl.f\
	./ppm_alloc_4dl.f\
	./ppm_alloc_5dl.f\
	./ppm_alloc_1d.f\
	./ppm_alloc_2d.f\
	./ppm_alloc_3d.f\
	./ppm_alloc_4d.f\
	./ppm_alloc_5d.f\
	./ppm_alloc_1dl.f\
	./ppm_alloc_2dl.f\
	./ppm_alloc_3dl.f\
	./ppm_alloc_4dl.f\
	./ppm_alloc_5dl.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_alloc.f __ppm_module_alloc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_alloc.f -o ppm_module_alloc.o
ppm_module_bem_basis.o: ppm_module_bem_basis.f\
	./ppm_bem_basis.f\
	./ppm_bem_basis.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_bem_basis.f __ppm_module_bem_basis.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_bem_basis.f -o ppm_module_bem_basis.o
ppm_module_bem_quadrule.o: ppm_module_bem_quadrule.f\
	./ppm_bem_quadrule.f\
	./ppm_bem_quadrule.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_bem_quadrule.f __ppm_module_bem_quadrule.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_bem_quadrule.f -o ppm_module_bem_quadrule.o
ppm_module_bem_quadrule_npoints.o: ppm_module_bem_quadrule_npoints.f\
	./ppm_bem_quadrule_npoints.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_bem_quadrule_npoints.f __ppm_module_bem_quadrule_npoints.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_bem_quadrule_npoints.f -o ppm_module_bem_quadrule_npoints.o
ppm_module_comp_mesh_fd.o: ppm_module_comp_mesh_fd.f
	$(CPP) $(DEFINE)  -P ppm_module_comp_mesh_fd.f __ppm_module_comp_mesh_fd.f
	$(F77) $(FFLAGS)  -c __ppm_module_comp_mesh_fd.f -o ppm_module_comp_mesh_fd.o
ppm_module_comp_mesh_fft.o: ppm_module_comp_mesh_fft.f
	$(CPP) $(DEFINE)  -P ppm_module_comp_mesh_fft.f __ppm_module_comp_mesh_fft.f
	$(F77) $(FFLAGS)  -c __ppm_module_comp_mesh_fft.f -o ppm_module_comp_mesh_fft.o
ppm_module_comp_pp_cell.o: ppm_module_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f\
	./ppm_comp_pp_cell.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_comp_pp_cell.f __ppm_module_comp_pp_cell.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_comp_pp_cell.f -o ppm_module_comp_pp_cell.o
ppm_module_comp_pp_correct.o: ppm_module_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f\
	./ppm_comp_pp_correct.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_comp_pp_correct.f __ppm_module_comp_pp_correct.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_comp_pp_correct.f -o ppm_module_comp_pp_correct.o
ppm_module_comp_pp_doring.o: ppm_module_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f\
	./ppm_comp_pp_doring.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_comp_pp_doring.f __ppm_module_comp_pp_doring.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_comp_pp_doring.f -o ppm_module_comp_pp_doring.o
ppm_module_comp_pp_mk_table.o: ppm_module_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f\
	./ppm_comp_pp_mk_table.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_comp_pp_mk_table.f __ppm_module_comp_pp_mk_table.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_comp_pp_mk_table.f -o ppm_module_comp_pp_mk_table.o
ppm_module_comp_pp_ring.o: ppm_module_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f\
	./ppm_comp_pp_ring.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_comp_pp_ring.f __ppm_module_comp_pp_ring.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_comp_pp_ring.f -o ppm_module_comp_pp_ring.o
ppm_module_comp_pp_verlet.o: ppm_module_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f\
	./ppm_comp_pp_verlet.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_comp_pp_verlet.f __ppm_module_comp_pp_verlet.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_comp_pp_verlet.f -o ppm_module_comp_pp_verlet.o
ppm_module_copy_image_subs.o: ppm_module_copy_image_subs.f\
	./ppm_copy_image_subs.f\
	./ppm_copy_image_subs.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_copy_image_subs.f __ppm_module_copy_image_subs.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_copy_image_subs.f -o ppm_module_copy_image_subs.o
ppm_module_data_fieldsolver.o: ppm_module_data_fieldsolver.f
	$(CPP) $(DEFINE) $(INC) -P ppm_module_data_fieldsolver.f __ppm_module_data_fieldsolver.f
	$(F77) $(FFLAGS)  -c __ppm_module_data_fieldsolver.f -o ppm_module_data_fieldsolver.o
ppm_module_data_ode.o: ppm_module_data_ode.f
	$(F77) $(FFLAGS)  -c ppm_module_data_ode.f
ppm_module_data_gmm.o: ppm_module_data_gmm.f
	$(F77) $(FFLAGS)  -c ppm_module_data_gmm.f
ppm_module_data.o: ppm_module_data.f\
        ./ppm_param.h
	$(F77) $(FFLAGS)  -I./ -c ppm_module_data.f
ppm_module_data_fmm.o: ppm_module_data_fmm.f
	$(F77) $(FFLAGS)  -c ppm_module_data_fmm.f
ppm_module_data_io.o: ppm_module_data_io.f
	$(F77) $(FFLAGS)  -c ppm_module_data_io.f
ppm_module_data_loadbal.o: ppm_module_data_loadbal.f
	$(F77) $(FFLAGS)  -c ppm_module_data_loadbal.f
ppm_module_data_mesh.o: ppm_module_data_mesh.f
	$(F77) $(FFLAGS)  -c ppm_module_data_mesh.f
ppm_module_data_mg.o: ppm_module_data_mg.f
	$(CPP) $(DEFINE)  -P ppm_module_data_mg.f __ppm_module_data_mg.f
	$(F77) $(FFLAGS)  -c __ppm_module_data_mg.f -o ppm_module_data_mg.o
ppm_module_data_neighlist.o: ppm_module_data_neighlist.f
	$(F77) $(FFLAGS)  -c ppm_module_data_neighlist.f
ppm_module_data_part.o: ppm_module_data_part.f
	$(F77) $(FFLAGS)  -c ppm_module_data_part.f
ppm_module_data_rmsh.o: ppm_module_data_rmsh.f\
	./ppm_module_data.f
	$(F77) $(FFLAGS)  -c ppm_module_data_rmsh.f
ppm_module_data_tree.o: ppm_module_data_tree.f
	$(F77) $(FFLAGS)  -c ppm_module_data_tree.f
ppm_module_gmm_init.o: ppm_module_gmm_init.f\
	./ppm_gmm_init.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_init.f __ppm_module_gmm_init.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_init.f -o ppm_module_gmm_init.o
ppm_module_gmm_finalize.o: ppm_module_gmm_finalize.f\
	./ppm_gmm_finalize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_finalize.f __ppm_module_gmm_finalize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_finalize.f -o ppm_module_gmm_finalize.o
ppm_module_gmm_kickoff.o: ppm_module_gmm_kickoff.f\
	./ppm_gmm_kickoff.f\
	./ppm_gmm_lubksb.inc\
	./ppm_gmm_quasinewton.inc
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_kickoff.f __ppm_module_gmm_kickoff.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_kickoff.f -o ppm_module_gmm_kickoff.o
ppm_module_gmm_extend_bkwd.o: ppm_module_gmm_extend_bkwd.f\
	./ppm_gmm_extend_bkwd.f\
	./ppm_gmm_add_to_list.inc\
	./ppm_gmm_switches.inc\
	./ppm_gmm_getfdta.inc\
	./ppm_gmm_getdta.inc\
	./ppm_gmm_slvextn.inc
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_extend_bkwd.f __ppm_module_gmm_extend_bkwd.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_extend_bkwd.f -o ppm_module_gmm_extend_bkwd.o
ppm_module_gmm_extend_fwd.o: ppm_module_gmm_extend_fwd.f\
	./ppm_gmm_extend_fwd.f\
	./ppm_gmm_add_to_list.inc\
	./ppm_gmm_switches.inc\
	./ppm_gmm_getfdta.inc\
	./ppm_gmm_getdta.inc\
	./ppm_gmm_slvextn.inc
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_extend_fwd.f __ppm_module_gmm_extend_fwd.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_extend_fwd.f -o ppm_module_gmm_extend_fwd.o
ppm_module_gmm_march_bkwd.o: ppm_module_gmm_march_bkwd.f\
	./ppm_gmm_march_bkwd.f\
	./ppm_gmm_add_to_list.inc\
	./ppm_gmm_quadeq.inc\
	./ppm_gmm_switches.inc\
	./ppm_gmm_getfdta.inc\
	./ppm_gmm_slvupwnd.inc
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_march_bkwd.f __ppm_module_gmm_march_bkwd.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_march_bkwd.f -o ppm_module_gmm_march_bkwd.o
ppm_module_gmm_march_fwd.o: ppm_module_gmm_march_fwd.f\
	./ppm_gmm_march_fwd.f\
	./ppm_gmm_add_to_list.inc\
	./ppm_gmm_quadeq.inc\
	./ppm_gmm_switches.inc\
	./ppm_gmm_getfdta.inc\
	./ppm_gmm_slvupwnd.inc
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_march_fwd.f __ppm_module_gmm_march_fwd.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_march_fwd.f -o ppm_module_gmm_march_fwd.o
ppm_module_gmm_march.o: ppm_module_gmm_march.f\
	./ppm_gmm_march.f\
	./ppm_module_gmm_march_fwd.o\
	./ppm_module_gmm_march_bkwd.o\
	./ppm_module_gmm_extend_fwd.o\
	./ppm_module_gmm_extend_bkwd.o
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_march.f __ppm_module_gmm_march.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_march.f -o ppm_module_gmm_march.o
ppm_module_gmm_reinitialize.o: ppm_module_gmm_reinitialize.f\
	./ppm_gmm_reinitialize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_reinitialize.f __ppm_module_gmm_reinitialize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_reinitialize.f -o ppm_module_gmm_reinitialize.o
ppm_module_gmm_cpt.o: ppm_module_gmm_cpt.f\
	./ppm_gmm_cpt.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_cpt.f __ppm_module_gmm_cpt.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_cpt.f -o ppm_module_gmm_cpt.o
ppm_module_gmm_extend.o: ppm_module_gmm_extend.f\
	./ppm_gmm_extend.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_gmm_extend.f __ppm_module_gmm_extend.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_gmm_extend.f -o ppm_module_gmm_extend.o
ppm_module_decomp_boxsplit.o: ppm_module_decomp_boxsplit.f\
	./ppm_decomp_boxsplit.f\
	./ppm_decomp_boxsplit.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_decomp_boxsplit.f __ppm_module_decomp_boxsplit.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_decomp_boxsplit.f -o ppm_module_decomp_boxsplit.o
ppm_module_check_meshid.o: ppm_module_check_meshid.f\
	./ppm_check_meshid.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_check_meshid.f __ppm_module_check_meshid.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_check_meshid.f -o ppm_module_check_meshid.o
ppm_module_check_topoid.o: ppm_module_check_topoid.f\
	./ppm_check_topoid.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_check_topoid.f __ppm_module_check_topoid.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_check_topoid.f -o ppm_module_check_topoid.o
ppm_module_decomp_cartesian.o: ppm_module_decomp_cartesian.f\
	./ppm_decomp_cartesian.f\
	./ppm_decomp_cartesian.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_decomp_cartesian.f __ppm_module_decomp_cartesian.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_decomp_cartesian.f -o ppm_module_decomp_cartesian.o
ppm_module_decomp_null.o: ppm_module_decomp_null.f\
	./ppm_decomp_null.f\
	./ppm_decomp_null.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_decomp_null.f __ppm_module_decomp_null.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_decomp_null.f -o ppm_module_decomp_null.o
ppm_module_decomp_pruned_cell.o: ppm_module_decomp_pruned_cell.f\
	./ppm_decomp_pruned_cell.f\
	./ppm_decomp_pruned_cell.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_decomp_pruned_cell.f __ppm_module_decomp_pruned_cell.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_decomp_pruned_cell.f -o ppm_module_decomp_pruned_cell.o
ppm_module_decomp_tree.o: ppm_module_decomp_tree.f\
	./ppm_decomp_tree.f\
	./ppm_decomp_tree.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_decomp_tree.f __ppm_module_decomp_tree.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_decomp_tree.f -o ppm_module_decomp_tree.o
ppm_module_define_subs_bc.o: ppm_module_define_subs_bc.f\
	./ppm_define_subs_bc.f\
	./ppm_define_subs_bc.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_define_subs_bc.f __ppm_module_define_subs_bc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_define_subs_bc.f -o ppm_module_define_subs_bc.o
ppm_module_doio.o: ppm_module_doio.f\
	./ppm_doio.f\
	./ppm_doio.f\
	./ppm_doio.f\
	./ppm_doio.f\
	./ppm_doio.f\
	./ppm_doio.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_doio.f __ppm_module_doio.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_doio.f -o ppm_module_doio.o
ppm_module_error.o: ppm_module_error.f\
	./ppm_error.h\
	./ppm_error.f\
	./ppm_module_write.o
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_error.f __ppm_module_error.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_error.f -o ppm_module_error.o
ppm_module_estimate_proc_speed.o: ppm_module_estimate_proc_speed.f\
	./ppm_estimate_proc_speed.f\
	./ppm_estimate_proc_speed.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_estimate_proc_speed.f __ppm_module_estimate_proc_speed.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_estimate_proc_speed.f -o ppm_module_estimate_proc_speed.o
ppm_module_fdsolver_fft_fd.o: ppm_module_fdsolver_fft_fd.f\
	./ppm_fdsolver_fft_fd_2d.f\
	./ppm_fdsolver_fft_fd_2d.f\
	./ppm_fdsolver_fft_fd_3d.f\
	./ppm_fdsolver_fft_fd_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_fft_fd.f   __ppm_module_fdsolver_fft_fd.f  
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_fft_fd.f   -o ppm_module_fdsolver_fft_fd.o  
ppm_module_fdsolver_fft_bd.o: ppm_module_fdsolver_fft_bd.f\
	./ppm_fdsolver_fft_bd_2d.f\
	./ppm_fdsolver_fft_bd_2d.f\
	./ppm_fdsolver_fft_bd_3d.f\
	./ppm_fdsolver_fft_bd_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_fft_bd.f   __ppm_module_fdsolver_fft_bd.f  
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_fft_bd.f   -o ppm_module_fdsolver_fft_bd.o  
ppm_module_fdsolver_finalize.o: ppm_module_fdsolver_finalize.f\
	./ppm_fdsolver_finalize.f\
	./ppm_fdsolver_finalize.f\
	./ppm_fdsolver_finalize.f\
	./ppm_fdsolver_finalize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_finalize.f __ppm_module_fdsolver_finalize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_finalize.f -o ppm_module_fdsolver_finalize.o
ppm_module_fdsolver_init.o: ppm_module_fdsolver_init.f\
	./ppm_fdsolver_init.f\
	./ppm_fdsolver_init.f\
	./ppm_fdsolver_init.f\
	./ppm_fdsolver_init.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_init.f __ppm_module_fdsolver_init.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_init.f -o ppm_module_fdsolver_init.o
ppm_module_fdsolver_map.o: ppm_module_fdsolver_map.f\
	./ppm_fdsolver_map_2d.f\
	./ppm_fdsolver_map_3d.f\
	./ppm_fdsolver_map_2d.f\
	./ppm_fdsolver_map_3d.f\
	./ppm_fdsolver_map_2d.f\
	./ppm_fdsolver_map_3d.f\
	./ppm_fdsolver_map_2d.f\
	./ppm_fdsolver_map_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_map.f __ppm_module_fdsolver_map.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_map.f -o ppm_module_fdsolver_map.o
ppm_module_fdsolver_poisson.o: ppm_module_fdsolver_poisson.f\
	./ppm_fdsolver_poisson_2d.f\
	./ppm_fdsolver_poisson_3d.f\
	./ppm_fdsolver_poisson_2d.f\
	./ppm_fdsolver_poisson_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_poisson.f __ppm_module_fdsolver_poisson.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_poisson.f -o ppm_module_fdsolver_poisson.o
ppm_module_fdsolver_solve.o: ppm_module_fdsolver_solve.f\
	./ppm_fdsolver_solve_2d.f\
	./ppm_fdsolver_solve_3d.f\
	./ppm_fdsolver_solve_2d.f\
	./ppm_fdsolver_solve_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_fdsolver_solve.f __ppm_module_fdsolver_solve.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_fdsolver_solve.f -o ppm_module_fdsolver_solve.o
ppm_module_finalize.o: ppm_module_finalize.f\
	./ppm_finalize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_finalize.f __ppm_module_finalize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_finalize.f -o ppm_module_finalize.o
ppm_module_find_duplicates.o: ppm_module_find_duplicates.f\
	./ppm_find_duplicates.f\
	./ppm_find_duplicates.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_find_duplicates.f __ppm_module_find_duplicates.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_find_duplicates.f -o ppm_module_find_duplicates.o
ppm_module_find_neigh.o: ppm_module_find_neigh.f\
	./ppm_find_neigh.f\
	./ppm_find_neigh.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_find_neigh.f __ppm_module_find_neigh.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_find_neigh.f -o ppm_module_find_neigh.o
ppm_module_get_cost.o: ppm_module_get_cost.f\
	./ppm_get_cost.f\
	./ppm_get_cost.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_get_cost.f __ppm_module_get_cost.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_get_cost.f -o ppm_module_get_cost.o
ppm_module_get_revision.o: ppm_module_get_revision.f\
	./ppm_get_revision.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_get_revision.f __ppm_module_get_revision.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_get_revision.f -o ppm_module_get_revision.o
ppm_module_init.o: ppm_module_init.f\
	./ppm_init.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_init.f __ppm_module_init.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_init.f -o ppm_module_init.o
ppm_module_io_close.o: ppm_module_io_close.f\
	./ppm_io_close.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_close.f __ppm_module_io_close.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_close.f -o ppm_module_io_close.o
ppm_module_io_delete.o: ppm_module_io_delete.f\
	./ppm_io_delete.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_delete.f __ppm_module_io_delete.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_delete.f -o ppm_module_io_delete.o
ppm_module_io.o: ppm_module_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f\
	./ppm_io.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io.f __ppm_module_io.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io.f -o ppm_module_io.o
ppm_module_io_inquire.o: ppm_module_io_inquire.f\
	./ppm_io_inquire.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_inquire.f __ppm_module_io_inquire.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_inquire.f -o ppm_module_io_inquire.o
ppm_module_io_open.o: ppm_module_io_open.f\
	./ppm_io_open.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_open.f __ppm_module_io_open.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_open.f -o ppm_module_io_open.o
ppm_module_io_read_ascii.o: ppm_module_io_read_ascii.f\
	./ppm_io_read_ascii.f\
	./ppm_io_read_ascii.f\
	./ppm_io_read_ascii.f\
	./ppm_io_read_ascii.f\
	./ppm_io_read_ascii.f\
	./ppm_io_read_ascii.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_read_ascii.f __ppm_module_io_read_ascii.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_read_ascii.f -o ppm_module_io_read_ascii.o
ppm_module_io_read_binary.o: ppm_module_io_read_binary.f\
	./ppm_io_read_binary.f\
	./ppm_io_read_binary.f\
	./ppm_io_read_binary.f\
	./ppm_io_read_binary.f\
	./ppm_io_read_binary.f\
	./ppm_io_read_binary.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_read_binary.f __ppm_module_io_read_binary.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_read_binary.f -o ppm_module_io_read_binary.o
ppm_module_io_set_unit.o: ppm_module_io_set_unit.f\
	./ppm_io_set_unit.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_set_unit.f __ppm_module_io_set_unit.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_set_unit.f -o ppm_module_io_set_unit.o
ppm_module_io_unused_unit.o: ppm_module_io_unused_unit.f\
	./ppm_io_unused_unit.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_unused_unit.f __ppm_module_io_unused_unit.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_unused_unit.f -o ppm_module_io_unused_unit.o
ppm_module_io_write_ascii.o: ppm_module_io_write_ascii.f\
	./ppm_io_write_ascii.f\
	./ppm_io_write_ascii.f\
	./ppm_io_write_ascii.f\
	./ppm_io_write_ascii.f\
	./ppm_io_write_ascii.f\
	./ppm_io_write_ascii.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_write_ascii.f __ppm_module_io_write_ascii.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_write_ascii.f -o ppm_module_io_write_ascii.o
ppm_module_io_write_binary.o: ppm_module_io_write_binary.f\
	./ppm_io_write_binary.f\
	./ppm_io_write_binary.f\
	./ppm_io_write_binary.f\
	./ppm_io_write_binary.f\
	./ppm_io_write_binary.f\
	./ppm_io_write_binary.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_io_write_binary.f __ppm_module_io_write_binary.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_io_write_binary.f -o ppm_module_io_write_binary.o
ppm_module_is_initialized.o: ppm_module_is_initialized.f\
	./ppm_is_initialized.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_is_initialized.f __ppm_module_is_initialized.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_is_initialized.f -o ppm_module_is_initialized.o
ppm_module_loadbal_inquire.o: ppm_module_loadbal_inquire.f\
	./ppm_loadbal_inquire.f\
	./ppm_loadbal_inquire.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_loadbal_inquire.f __ppm_module_loadbal_inquire.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_loadbal_inquire.f -o ppm_module_loadbal_inquire.o
ppm_module_log.o: ppm_module_log.f\
	./ppm_log.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_log.f __ppm_module_log.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_log.f -o ppm_module_log.o
ppm_module_map_connect_distrib.o: ppm_module_map_connect_distrib.f\
	./ppm_map_connect_distrib.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_connect_distrib.f __ppm_module_map_connect_distrib.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_connect_distrib.f -o ppm_module_map_connect_distrib.o
ppm_module_map_connect.o: ppm_module_map_connect.f\
	./ppm_map_connect.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_connect.f __ppm_module_map_connect.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_connect.f -o ppm_module_map_connect.o
ppm_module_map_connect_prune.o: ppm_module_map_connect_prune.f\
	./ppm_map_connect_prune.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_connect_prune.f __ppm_module_map_connect_prune.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_connect_prune.f -o ppm_module_map_connect_prune.o
ppm_module_map_connect_send.o: ppm_module_map_connect_send.f\
	./ppm_map_connect_send.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_connect_send.f __ppm_module_map_connect_send.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_connect_send.f -o ppm_module_map_connect_send.o
ppm_module_map_field.o: ppm_module_map_field.f\
	./ppm_map_field.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field.f __ppm_module_map_field.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field.f -o ppm_module_map_field.o
ppm_module_map_field_ghost.o: ppm_module_map_field_ghost.f\
	./ppm_map_field_ghost.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_ghost.f __ppm_module_map_field_ghost.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_ghost.f -o ppm_module_map_field_ghost.o
ppm_module_map_field_ghost_get.o: ppm_module_map_field_ghost_get.f\
	./ppm_map_field_ghost_get.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_ghost_get.f __ppm_module_map_field_ghost_get.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_ghost_get.f -o ppm_module_map_field_ghost_get.o
ppm_module_map_field_ghost_init.o: ppm_module_map_field_ghost_init.f\
	./ppm_map_field_ghost_init.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_ghost_init.f __ppm_module_map_field_ghost_init.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_ghost_init.f -o ppm_module_map_field_ghost_init.o
ppm_module_map_field_ghost_put.o: ppm_module_map_field_ghost_put.f\
	./ppm_map_field_ghost_put.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_ghost_put.f __ppm_module_map_field_ghost_put.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_ghost_put.f -o ppm_module_map_field_ghost_put.o
ppm_module_map_field_global.o: ppm_module_map_field_global.f\
	./ppm_map_field_global.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_global.f __ppm_module_map_field_global.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_global.f -o ppm_module_map_field_global.o
ppm_module_map_field_pop_2d.o: ppm_module_map_field_pop_2d.f\
	./ppm_map_field_pop_2d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_pop_2d.f __ppm_module_map_field_pop_2d.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_pop_2d.f -o ppm_module_map_field_pop_2d.o
ppm_module_map_field_pop_3d.o: ppm_module_map_field_pop_3d.f\
	./ppm_map_field_pop_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_pop_3d.f __ppm_module_map_field_pop_3d.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_pop_3d.f -o ppm_module_map_field_pop_3d.o
ppm_module_map_field_push_2d.o: ppm_module_map_field_push_2d.f\
	./ppm_map_field_push_2d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_push_2d.f __ppm_module_map_field_push_2d.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_push_2d.f -o ppm_module_map_field_push_2d.o
ppm_module_map_field_push_3d.o: ppm_module_map_field_push_3d.f\
	./ppm_map_field_push_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_push_3d.f __ppm_module_map_field_push_3d.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_push_3d.f -o ppm_module_map_field_push_3d.o
ppm_module_map_field_send.o: ppm_module_map_field_send.f\
	./ppm_map_field_send.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_field_send.f __ppm_module_map_field_send.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_field_send.f -o ppm_module_map_field_send.o
ppm_module_map_part.o: ppm_module_map_part.f\
	./ppm_map_part.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part.f __ppm_module_map_part.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part.f -o ppm_module_map_part.o
ppm_module_map_part_ghost.o: ppm_module_map_part_ghost.f\
	./ppm_map_part_ghost.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_ghost.f __ppm_module_map_part_ghost.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_ghost.f -o ppm_module_map_part_ghost.o
ppm_module_map_part_get_sub.o: ppm_module_map_part_get_sub.f\
        ./ppm_map_part_get_sub.f
	$(CPP) $(DEFINE)  -I./ -P ppm_module_map_part_get_sub.f __ppm_module_map_part_get_sub.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_map_part_get_sub.f -o ppm_module_map_part_get_sub.o
ppm_module_map_part_ghost_get.o: ppm_module_map_part_ghost_get.f\
	./ppm_map_part_ghost_get.f\
	./ppm_map_part_ghost_get.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_ghost_get.f __ppm_module_map_part_ghost_get.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_ghost_get.f -o ppm_module_map_part_ghost_get.o
ppm_module_map_part_ghost_put.o: ppm_module_map_part_ghost_put.f\
	./ppm_map_part_ghost_put.f\
	./ppm_map_part_ghost_put.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_ghost_put.f __ppm_module_map_part_ghost_put.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_ghost_put.f -o ppm_module_map_part_ghost_put.o
ppm_module_map_part_global.o: ppm_module_map_part_global.f\
	./ppm_map_part_global.f\
	./ppm_map_part_global.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_global.f __ppm_module_map_part_global.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_global.f -o ppm_module_map_part_global.o
ppm_module_map_part_partial.o: ppm_module_map_part_partial.f\
	./ppm_map_part_partial.f\
	./ppm_map_part_partial.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_partial.f __ppm_module_map_part_partial.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_partial.f -o ppm_module_map_part_partial.o
ppm_module_map_part_pop.o: ppm_module_map_part_pop.f\
	./ppm_map_part_pop.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_pop.f __ppm_module_map_part_pop.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_pop.f -o ppm_module_map_part_pop.o
ppm_module_map_part_push.o: ppm_module_map_part_push.f\
	./ppm_map_part_push.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_push.f __ppm_module_map_part_push.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_push.f -o ppm_module_map_part_push.o
ppm_module_map_part_remap.o: ppm_module_map_part_remap.f\
	./ppm_map_part_remap.f\
	./ppm_map_part_remap.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_remap.f __ppm_module_map_part_remap.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_remap.f -o ppm_module_map_part_remap.o
ppm_module_map_part_ring.o: ppm_module_map_part_ring.f\
	./ppm_map_part_ring.f\
	./ppm_map_part_ring.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_ring.f __ppm_module_map_part_ring.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_ring.f -o ppm_module_map_part_ring.o
ppm_module_map_part_ring_shift.o: ppm_module_map_part_ring_shift.f\
	./ppm_map_part_ring_shift.f\
	./ppm_map_part_ring_shift.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_ring_shift.f __ppm_module_map_part_ring_shift.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_ring_shift.f -o ppm_module_map_part_ring_shift.o
ppm_module_map_part_send.o: ppm_module_map_part_send.f\
	./ppm_map_part_send.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_map_part_send.f __ppm_module_map_part_send.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_map_part_send.f -o ppm_module_map_part_send.o
ppm_module_mesh_alloc.o: ppm_module_mesh_alloc.f\
	./ppm_mesh_alloc_list.f\
	./ppm_mesh_alloc_equi.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_alloc.f __ppm_module_mesh_alloc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_alloc.f -o ppm_module_mesh_alloc.o
ppm_module_mesh_block_intersect.o: ppm_module_mesh_block_intersect.f\
	./ppm_mesh_block_intersect.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_block_intersect.f __ppm_module_mesh_block_intersect.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_block_intersect.f -o ppm_module_mesh_block_intersect.o
ppm_module_mesh_define.o: ppm_module_mesh_define.f\
	./ppm_mesh_define.f\
	./ppm_mesh_define.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_define.f __ppm_module_mesh_define.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_define.f -o ppm_module_mesh_define.o
ppm_module_mesh_derive.o: ppm_module_mesh_derive.f\
	./ppm_mesh_derive.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_derive.f __ppm_module_mesh_derive.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_derive.f -o ppm_module_mesh_derive.o
ppm_module_mesh_finalize.o: ppm_module_mesh_finalize.f\
	./ppm_mesh_finalize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_finalize.f __ppm_module_mesh_finalize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_finalize.f -o ppm_module_mesh_finalize.o
ppm_module_mesh_on_subs.o: ppm_module_mesh_on_subs.f\
	./ppm_mesh_on_subs.f\
	./ppm_mesh_on_subs.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_on_subs.f __ppm_module_mesh_on_subs.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_on_subs.f -o ppm_module_mesh_on_subs.o
ppm_module_mesh_store.o: ppm_module_mesh_store.f\
	./ppm_mesh_store.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mesh_store.f __ppm_module_mesh_store.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mesh_store.f -o ppm_module_mesh_store.o
ppm_module_topo_metis_s2p.o: ppm_module_topo_metis_s2p.f\
	./ppm_topo_metis_s2p.f\
	./ppm_topo_metis_s2p.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_topo_metis_s2p.f __ppm_module_topo_metis_s2p.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_topo_metis_s2p.f -o ppm_module_topo_metis_s2p.o
ppm_module_mg_alloc.o: ppm_module_mg_alloc.f\
	./ppm_mg_alloc_field.f\
	./ppm_mg_alloc_bc.f\
	./ppm_mg_alloc_field.f\
	./ppm_mg_alloc_bc.f\
	./ppm_mg_alloc_field.f\
	./ppm_mg_alloc_bc.f\
	./ppm_mg_alloc_field.f\
	./ppm_mg_alloc_bc.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_alloc.f __ppm_module_mg_alloc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_alloc.f -o ppm_module_mg_alloc.o
ppm_module_mktopo.o: ppm_module_mktopo.f\
	./ppm_topo_mkpart.f\
	./ppm_topo_mkfield.f\
	./ppm_topo_mkgeom.f\
	./ppm_topo_mktree.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mktopo.f __ppm_module_mktopo.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mktopo.f -o ppm_module_mktopo.o
ppm_module_clist_destroy.o: ppm_module_clist_destroy.f\
	./ppm_clist_destroy.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_clist_destroy.f __ppm_module_clist_destroy.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_clist_destroy.f -o ppm_module_clist_destroy.o
ppm_module_neighlist_clist.o: ppm_module_neighlist_clist.f\
	./ppm_neighlist_clist.f\
	./ppm_neighlist_clist.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_neighlist_clist.f __ppm_module_neighlist_clist.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_neighlist_clist.f -o ppm_module_neighlist_clist.o
ppm_module_neighlist_MkNeighIdx.o: ppm_module_neighlist_MkNeighIdx.f\
	./ppm_neighlist_MkNeighIdx.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_neighlist_MkNeighIdx.f __ppm_module_neighlist_MkNeighIdx.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_neighlist_MkNeighIdx.f -o ppm_module_neighlist_MkNeighIdx.o
ppm_module_neighlist_vlist.o: ppm_module_neighlist_vlist.f\
	./ppm_neighlist_vlist.f\
	./ppm_neighlist_vlist.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_neighlist_vlist.f __ppm_module_neighlist_vlist.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_neighlist_vlist.f -o ppm_module_neighlist_vlist.o
ppm_module_ode_alldone.o: ppm_module_ode_alldone.f\
	./ppm_ode_alldone.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_alldone.f __ppm_module_ode_alldone.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_alldone.f -o ppm_module_ode_alldone.o
ppm_module_ode_create_ode.o: ppm_module_ode_create_ode.f\
	./ppm_ode_create_ode.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_create_ode.f __ppm_module_ode_create_ode.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_create_ode.f -o ppm_module_ode_create_ode.o
ppm_module_impose_part_bc.o: ppm_module_impose_part_bc.f\
	./ppm_impose_part_bc.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_impose_part_bc.f __ppm_module_impose_part_bc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_impose_part_bc.f -o ppm_module_impose_part_bc.o
ppm_module_ode_finalize.o: ppm_module_ode_finalize.f\
	./ppm_ode_finalize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_finalize.f __ppm_module_ode_finalize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_finalize.f -o ppm_module_ode_finalize.o
ppm_module_ode_init.o: ppm_module_ode_init.f\
	./ppm_ode_init.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_init.f __ppm_module_ode_init.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_init.f -o ppm_module_ode_init.o
ppm_module_ode_map_pop.o: ppm_module_ode_map_pop.f\
	./ppm_ode_map_pop.f\
	./ppm_ode_map_pop.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_map_pop.f __ppm_module_ode_map_pop.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_map_pop.f -o ppm_module_ode_map_pop.o
ppm_module_ode_map_push.o: ppm_module_ode_map_push.f\
	./ppm_ode_map_push.f\
	./ppm_ode_map_push.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_map_push.f __ppm_module_ode_map_push.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_map_push.f -o ppm_module_ode_map_push.o
ppm_module_ode_start.o: ppm_module_ode_start.f\
	./ppm_ode_start.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_ode_start.f __ppm_module_ode_start.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_ode_start.f -o ppm_module_ode_start.o
ppm_module_ode_step.o: ppm_module_ode_step.f\
        ./ppm_ode_step.f\
        ./ppm_ode_step.f\
        ./ppm_ode_step.f\
        ./ppm_ode_step.f
	$(CPP) $(DEFINE)  -I./ -P ppm_module_ode_step.f __ppm_module_ode_step.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_ode_step.f -o ppm_module_ode_step.o
ppm_module_fmm_finalize.o: ppm_module_fmm_finalize.f\
	./ppm_fmm_finalize.f	
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_finalize.f __ppm_module_fmm_finalize.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_finalize.f -o ppm_module_fmm_finalize.o
ppm_module_fmm_init.o: ppm_module_fmm_init.f\
	./ppm_fmm_init.f	
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_init.f __ppm_module_fmm_init.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_init.f -o ppm_module_fmm_init.o
ppm_module_fmm_pretraverse.o: ppm_module_fmm_pretraverse.f\
	./ppm_fmm_pretraverse.f	
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_pretraverse.f __ppm_module_fmm_pretraverse.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_pretraverse.f -o ppm_module_fmm_pretraverse.o
ppm_module_fmm_traverse.o: ppm_module_fmm_traverse.f\
	./ppm_fmm_traverse.f	
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_traverse.f __ppm_module_fmm_traverse.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_traverse.f -o ppm_module_fmm_traverse.o
ppm_module_fmm_expchange.o: ppm_module_fmm_expchange.f\
	./ppm_fmm_expchange.f
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_expchange.f __ppm_module_fmm_expchange.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_expchange.f -o ppm_module_fmm_expchange.o
ppm_module_fmm_expansion.o: ppm_module_fmm_expansion.f\
	./ppm_fmm_expansion.f	
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_expansion.f __ppm_module_fmm_expansion.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_expansion.f -o ppm_module_fmm_expansion.o
ppm_module_fmm_potential.o: ppm_module_fmm_potential.f\
	./ppm_fmm_potential.f	
	$(CPP) $(DEFINE)  -I./ -P ppm_module_fmm_potential.f __ppm_module_fmm_potential.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_fmm_potential.f -o ppm_module_fmm_potential.o
ppm_module_mg_finalize.o: ppm_module_mg_finalize.f\
	./ppm_mg_finalize.f\
	./ppm_mg_finalize.f\
	./ppm_mg_finalize.f\
	./ppm_mg_finalize.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_finalize.f __ppm_module_mg_finalize.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_finalize.f -o ppm_module_mg_finalize.o
ppm_module_mg_init.o: ppm_module_mg_init.f\
	./ppm_mg_init.f\
	./ppm_mg_init.f\
	./ppm_mg_init.f\
	./ppm_mg_init.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_init.f __ppm_module_mg_init.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_init.f -o ppm_module_mg_init.o
ppm_module_mg_solv.o: ppm_module_mg_solv.f\
	./ppm_mg_solv.f\
	./ppm_mg_solv.f\
	./ppm_mg_solv.f\
	./ppm_mg_solv.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_solv.f __ppm_module_mg_solv.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_solv.f -o ppm_module_mg_solv.o
ppm_module_mg_core.o: ppm_module_mg_core.f\
	./ppm_mg_core.f\
	./ppm_mg_core.f\
	./ppm_mg_core.f\
	./ppm_mg_core.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_core.f __ppm_module_mg_core.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_core.f -o ppm_module_mg_core.o
ppm_module_mg_smooth.o: ppm_module_mg_smooth.f\
	./ppm_mg_smooth_coarse.f\
	./ppm_mg_smooth_fine.f\
	./ppm_mg_smooth_coarse.f\
	./ppm_mg_smooth_fine.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_smooth.f __ppm_module_mg_smooth.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_smooth.f -o ppm_module_mg_smooth.o
ppm_module_mg_res.o: ppm_module_mg_res.f\
	./ppm_mg_res_coarse.f\
	./ppm_mg_res_fine.f\
	./ppm_mg_res_coarse.f\
	./ppm_mg_res_fine.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_res.f __ppm_module_mg_res.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_res.f -o ppm_module_mg_res.o
ppm_module_mg_prolong.o: ppm_module_mg_prolong.f\
	./ppm_mg_prolong.f\
	./ppm_mg_prolong.f\
	./ppm_mg_prolong.f\
	./ppm_mg_prolong.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_prolong.f __ppm_module_mg_prolong.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_prolong.f -o ppm_module_mg_prolong.o
ppm_module_mg_restrict.o: ppm_module_mg_restrict.f\
	./ppm_mg_restrict.f\
	./ppm_mg_restrict.f\
	./ppm_mg_restrict.f\
	./ppm_mg_restrict.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_mg_restrict.f __ppm_module_mg_restrict.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_mg_restrict.f -o ppm_module_mg_restrict.o
ppm_module_interp_p2m.o: ppm_module_interp_p2m.f\
	./ppm_interp_p2m.f\
	./ppm_module_data_rmsh.f
	$(CPP) $(DEFINE) $(INC) -P ppm_module_interp_p2m.f __ppm_module_interp_p2m.f
	$(F77) $(FFLAGS) $(INC) -c __ppm_module_interp_p2m.f -o ppm_module_interp_p2m.o
ppm_module_interp_m2p.o: ppm_module_interp_m2p.f\
	./ppm_interp_m2p.f\
	./ppm_module_data_rmsh.f
	$(CPP) $(DEFINE) $(INC) -P ppm_module_interp_m2p.f __ppm_module_interp_m2p.f
	$(F77) $(FFLAGS) $(INC) -c __ppm_module_interp_m2p.f -o ppm_module_interp_m2p.o
ppm_module_rmsh_comp_weights.o: ppm_module_rmsh_comp_weights.f\
	./ppm_rmsh_comp_weights.f\
	./ppm_module_data_rmsh.f
	$(CPP) $(DEFINE) $(INC) -P ppm_module_rmsh_comp_weights.f __ppm_module_rmsh_comp_weights.f
	$(F77) $(FFLAGS) $(INC) -c __ppm_module_rmsh_comp_weights.f -o ppm_module_rmsh_comp_weights.o
ppm_module_rmsh_remesh.o: ppm_module_rmsh_remesh.f\
	./ppm_rmsh_remesh.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_rmsh_remesh.f __ppm_module_rmsh_remesh.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_rmsh_remesh.f  -o ppm_module_rmsh_remesh.o
ppm_module_rmsh_create_part.o: ppm_module_rmsh_create_part.f\
	./ppm_rmsh_create_part.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_rmsh_create_part.f __ppm_module_rmsh_create_part.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_rmsh_create_part.f   -o ppm_module_rmsh_create_part.o
ppm_module_print_defines.o: ppm_module_print_defines.f\
	./ppm_print_defines.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_print_defines.f __ppm_module_print_defines.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_print_defines.f -o ppm_module_print_defines.o
ppm_module_set_decomp_cost.o: ppm_module_set_decomp_cost.f\
	./ppm_set_decomp_cost.f\
	./ppm_set_decomp_cost.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_set_decomp_cost.f __ppm_module_set_decomp_cost.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_set_decomp_cost.f -o ppm_module_set_decomp_cost.o
ppm_module_set_proc_speed.o: ppm_module_set_proc_speed.f\
	./ppm_set_proc_speed.f\
	./ppm_set_proc_speed.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_set_proc_speed.f __ppm_module_set_proc_speed.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_set_proc_speed.f -o ppm_module_set_proc_speed.o
ppm_module_topo_subs2proc.o: ppm_module_topo_subs2proc.f\
	./ppm_topo_subs2proc.f\
	./ppm_topo_subs2proc.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_topo_subs2proc.f __ppm_module_topo_subs2proc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_topo_subs2proc.f -o ppm_module_topo_subs2proc.o
ppm_module_substart.o: ppm_module_substart.f\
	./substart.f\
	./substart.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_substart.f __ppm_module_substart.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_substart.f -o ppm_module_substart.o
ppm_module_substop.o: ppm_module_substop.f\
	./substop.f\
	./substop.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_substop.f __ppm_module_substop.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_substop.f -o ppm_module_substop.o
ppm_module_time.o: ppm_module_time.f\
	./ppm_time.f\
	./ppm_time.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_time.f __ppm_module_time.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_time.f -o ppm_module_time.o
ppm_module_scale_domain.o: ppm_module_scale_domain.f\
	./ppm_scale_domain.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_scale_domain.f __ppm_module_scale_domain.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_scale_domain.f -o ppm_module_scale_domain.o
ppm_module_topo_check.o: ppm_module_topo_check.f\
	./ppm_topo_check.f\
	./ppm_topo_check.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_topo_check.f __ppm_module_topo_check.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_topo_check.f -o ppm_module_topo_check.o
ppm_module_topo_cost.o: ppm_module_topo_cost.f\
	./ppm_topo_cost.f\
	./ppm_topo_cost.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_topo_cost.f __ppm_module_topo_cost.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_topo_cost.f -o ppm_module_topo_cost.o
ppm_module_topo_inquire.o: ppm_module_topo_inquire.f\
	./ppm_topo_inquire.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_topo_inquire.f __ppm_module_topo_inquire.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_topo_inquire.f -o ppm_module_topo_inquire.o
ppm_module_topo_store.o: ppm_module_topo_store.f\
	./ppm_topo_store.f\
	./ppm_topo_store.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_topo_store.f __ppm_module_topo_store.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_topo_store.f -o ppm_module_topo_store.o
ppm_module_topo_box2subs.o: ppm_module_topo_box2subs.f\
        ./ppm_topo_box2subs.f
	$(CPP) $(DEFINE)  -I./ -P ppm_module_topo_box2subs.f __ppm_module_topo_box2subs.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_topo_box2subs.f -o ppm_module_topo_box2subs.o
ppm_module_util_cart2sph.o: ppm_module_util_cart2sph.f\
        ./ppm_util_cart2sph.f
	$(CPP) $(DEFINE)  -I./ -P ppm_module_util_cart2sph.f __ppm_module_util_cart2sph.f
	$(F77) $(FFLAGS)  -I./ -c __ppm_module_util_cart2sph.f -o ppm_module_util_cart2sph.o
ppm_module_util_commopt.o: ppm_module_util_commopt.f\
	./ppm_util_commopt.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_commopt.f __ppm_module_util_commopt.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_commopt.f -o ppm_module_util_commopt.o
ppm_module_util_fft_backward.o: ppm_module_util_fft_backward.f\
	./ppm_util_fft_backward_2d.f\
	./ppm_util_fft_backward_3d.f\
	./ppm_util_fft_backward_2d.f\
	./ppm_util_fft_backward_3d.f\
	./ppm_util_fft_backward_2d.f\
	./ppm_util_fft_backward_3d.f\
	./ppm_util_fft_backward_2d.f\
	./ppm_util_fft_backward_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_fft_backward.f __ppm_module_util_fft_backward.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_fft_backward.f -o ppm_module_util_fft_backward.o
ppm_module_util_fft_forward.o: ppm_module_util_fft_forward.f\
	./ppm_util_fft_forward_2d.f\
	./ppm_util_fft_forward_3d.f\
	./ppm_util_fft_forward_2d.f\
	./ppm_util_fft_forward_3d.f\
	./ppm_util_fft_forward_2d.f\
	./ppm_util_fft_forward_3d.f\
	./ppm_util_fft_forward_2d.f\
	./ppm_util_fft_forward_3d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_fft_forward.f __ppm_module_util_fft_forward.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_fft_forward.f -o ppm_module_util_fft_forward.o
ppm_module_util_invert_list.o: ppm_module_util_invert_list.f\
	./ppm_util_invert_list.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_invert_list.f __ppm_module_util_invert_list.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_invert_list.f -o ppm_module_util_invert_list.o
ppm_module_util_collapse_list.o: ppm_module_util_collapse_list.f\
	./ppm_util_collapse_list.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_collapse_list.f __ppm_module_util_collapse_list.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_collapse_list.f -o ppm_module_util_collapse_list.o
ppm_module_util_rank.o: ppm_module_util_rank.f\
	./ppm_util_rank3d.f\
	./ppm_util_rank2d.f\
	./ppm_util_rank3d.f\
	./ppm_util_rank2d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_rank.f __ppm_module_util_rank.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_rank.f -o ppm_module_util_rank.o
ppm_module_util_qsort.o: ppm_module_util_qsort.f\
	./ppm_util_qsort.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_qsort.f __ppm_module_util_qsort.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_qsort.f -o ppm_module_util_qsort.o
ppm_module_util_sort.o: ppm_module_util_sort.f\
	./ppm_util_sort3d.f\
	./ppm_util_sort2d.f\
	./ppm_util_sort3d.f\
	./ppm_util_sort2d.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_sort.f __ppm_module_util_sort.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_sort.f -o ppm_module_util_sort.o
ppm_module_util_matinv_2x2.o: ppm_module_util_matinv_2x2.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_matinv_2x2.f __ppm_module_util_matinv_2x2.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_matinv_2x2.f -o ppm_module_util_matinv_2x2.o
ppm_module_util_matinv_3x3.o: ppm_module_util_matinv_3x3.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_matinv_3x3.f __ppm_module_util_matinv_3x3.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_matinv_3x3.f -o ppm_module_util_matinv_3x3.o
ppm_module_util_quadeq_real.o: ppm_module_util_quadeq_real.f\
	./ppm_util_quadeq_real.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_quadeq_real.f __ppm_module_util_quadeq_real.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_quadeq_real.f -o ppm_module_util_quadeq_real.o
ppm_module_util_cubeq_real.o: ppm_module_util_cubeq_real.f\
	./ppm_util_cubeq_real.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_cubeq_real.f __ppm_module_util_cubeq_real.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_cubeq_real.f -o ppm_module_util_cubeq_real.o
ppm_module_util_eigen_2sym.o: ppm_module_util_eigen_2sym.f\
	./ppm_util_eigen_2sym.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_eigen_2sym.f __ppm_module_util_eigen_2sym.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_eigen_2sym.f -o ppm_module_util_eigen_2sym.o
ppm_module_util_eigen_3sym.o: ppm_module_util_eigen_3sym.f\
	./ppm_util_eigen_3sym.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_eigen_3sym.f __ppm_module_util_eigen_3sym.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_eigen_3sym.f -o ppm_module_util_eigen_3sym.o
ppm_module_tree_done.o: ppm_module_tree_done.f\
	./ppm_tree_done.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_done.f __ppm_module_tree_done.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_done.f -o ppm_module_tree_done.o
ppm_module_tree_alloc.o: ppm_module_tree_alloc.f\
	./ppm_tree_alloc.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_alloc.f __ppm_module_tree_alloc.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_alloc.f -o ppm_module_tree_alloc.o
ppm_module_tree_boxcost.o: ppm_module_tree_boxcost.f\
	./ppm_tree_boxcost.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_boxcost.f __ppm_module_tree_boxcost.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_boxcost.f -o ppm_module_tree_boxcost.o
ppm_module_tree_boxcut.o: ppm_module_tree_boxcut.f\
	./ppm_tree_boxcut.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_boxcut.f __ppm_module_tree_boxcut.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_boxcut.f -o ppm_module_tree_boxcut.o
ppm_module_tree_cutdir.o: ppm_module_tree_cutdir.f\
	./ppm_tree_cutdir.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_cutdir.f __ppm_module_tree_cutdir.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_cutdir.f -o ppm_module_tree_cutdir.o
ppm_module_tree_cutpos.o: ppm_module_tree_cutpos.f\
	./ppm_tree_cutpos.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_cutpos.f __ppm_module_tree_cutpos.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_cutpos.f -o ppm_module_tree_cutpos.o
ppm_module_tree_divcheck.o: ppm_module_tree_divcheck.f\
	./ppm_tree_divcheck.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree_divcheck.f __ppm_module_tree_divcheck.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree_divcheck.f -o ppm_module_tree_divcheck.o
ppm_module_tree.o: ppm_module_tree.f\
	./ppm_tree.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_tree.f __ppm_module_tree.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_tree.f -o ppm_module_tree.o
ppm_module_util_time.o: ppm_module_util_time.f\
	./ppm_util_time.f\
	./ppm_util_time.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_util_time.f __ppm_module_util_time.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_util_time.f -o ppm_module_util_time.o
ppm_module_write.o: ppm_module_write.f\
	./ppm_write.f
	$(CPP) $(DEFINE)  $(INC) -P ppm_module_write.f __ppm_module_write.f
	$(F77) $(FFLAGS)  $(INC) -c __ppm_module_write.f -o ppm_module_write.o
ppm_module_comp_part.o: ppm_module_comp_part.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_comp_part.f -o ppm_module_comp_part.o
ppm_module_comp_mesh.o: ppm_module_comp_mesh.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_comp_mesh.f -o ppm_module_comp_mesh.o
ppm_module_bem.o: ppm_module_bem.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_bem.f -o ppm_module_bem.o
ppm_module_fieldsolver.o: ppm_module_fieldsolver.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_fieldsolver.f -o ppm_module_fieldsolver.o
ppm_module_map.o: ppm_module_map.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_map.f -o ppm_module_map.o
ppm_module_user_util.o: ppm_module_user_util.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_user_util.f -o ppm_module_user_util.o
ppm_module_loadbal.o: ppm_module_loadbal.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_loadbal.f -o ppm_module_loadbal.o
ppm_module_topo.o: ppm_module_topo.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_topo.f -o ppm_module_topo.o
ppm_module_user_io.o: ppm_module_user_io.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_user_io.f -o ppm_module_user_io.o
ppm_module_ode.o: ppm_module_ode.f
	$(F77) $(FFLAGS)  -I./ -c ppm_module_ode.f -o ppm_module_ode.o
ppm_module_fmm.o: ppm_module_fmm.f
	$(F77) $(FFLAGS)  -I./ -c ppm_module_fmm.f -o ppm_module_fmm.o
ppm_module_mg.o: ppm_module_mg.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_mg.f -o ppm_module_mg.o
ppm_module_neighlist.o: ppm_module_neighlist.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_neighlist.f -o ppm_module_neighlist.o
ppm_module_rmsh.o: ppm_module_rmsh.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_rmsh.f -o ppm_module_rmsh.o
ppm_module_gmm.o: ppm_module_gmm.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_gmm.f -o ppm_module_gmm.o
ppm_module_user.o: ppm_module_user.f
	$(F77) $(FFLAGS)  $(INC) -c ppm_module_user.f -o ppm_module_user.o

#----------------------------------------------------------------------
#  link								
#----------------------------------------------------------------------
$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) libppm.a $(OBJS)

#----------------------------------------------------------------------
#  Install							     
#----------------------------------------------------------------------
install: $(BIN)
	(cp -f $(BIN) $(bindir))
	(ranlib $(bindir)/$(BIN))

#----------------------------------------------------------------------
#  Clean							       
#----------------------------------------------------------------------
new: cleanall libppm.a
cleanall:
	 rm -f __*.f
	 rm -f $(OBJS)
	 rm -f $(TARGET)
	 rm -f $(BIN)
	 rm -f *.o
	 rm -f *.mod
	 rm -f *.lst
	 rm -f *.l
	 rm -f *.L

clean:
	 rm -f __*.f
	 rm -f *.lst


