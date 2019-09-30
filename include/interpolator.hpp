#ifndef INTERPOLATOR_HS
#define INTERPOLATOR_HS

#include<iostream>
#include<cstddef>
#include<GetPot>

// getfem includes
#include "getfem/getfem_generic_assembly.h"
#include "getfem/getfem_export.h"
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_derivatives.h"
#include "getfem/getfem_superlu.h"
#include "getfem/getfem_models.h"
#include "getfem/getfem_model_solvers.h" // for preconditioners
#include "getfem/getfem_import.h"


#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"

#include "defines.hpp"
#include "node.hpp"


class INTERPOLATOR{
public:
INTERPOLATOR(const GetPot& df, getfem::mesh_fem& mf_1, 
                               getfem::mesh_fem& mf_2, getfem::mesh_im& im_2, std::vector<scalar_type>& radius );
void build_operator(); // build the operators for coupling

inline sparse_matrix_type& get_Mbar(){return Mbar_;}
inline sparse_matrix_type& get_Mlin(){return Mlin_;}

private:
getfem::mesh_fem *mf_1_, *mf_2_;
getfem::mesh_im*  mim_2_;
sparse_matrix_type Mbar_,Mlin_;
std::vector<scalar_type>* radius_;
int NInt_;//Nb of discretisation point for 3D-1D interpolation
}; //end class interpolator

#endif
