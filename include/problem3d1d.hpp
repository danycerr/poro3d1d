#ifndef PROBLEM3D1D_HS
#define PROBLEM3D1D_HS


#include"biot.hpp"
#include"network.hpp"
#include"interpolator.hpp"

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


class PROBLEM_3d1d{
public:
   PROBLEM_3d1d(const GetPot& df);
   void assembly_mat();
   void assembly_rhs();
   void assembly_exchange_mat();
   void solve();
   void print(int iter=0);
private:

   BIOT biot_problem;
   
   NETWORK oned_problem;

   INTERPOLATOR op_3d1d_; // 3d1d operator

   bool coupling_=true;

   std::vector<scalar_type> sol_, sol_old_, rhs_;    /// solution and old solution
   sparse_matrix_type K_;                         /// iteration matrix 	
   
   
   getfem::size_type nb_dof_biot_; // nb dof biot problem
   getfem::size_type nb_dof_biot_p_; // nb dof biot problem
   getfem::size_type nb_dof_biot_u_; // nb dof biot problem
   getfem::size_type nb_dof_oned_; // nb dof 1d problem
};//end problem 3d1d


#endif
