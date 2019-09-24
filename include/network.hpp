#ifndef NETWORKHS
#define NETWORKHS

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

class NETWORK
{
public:
    vector_size_type nb_vertices_;

    NETWORK(const GetPot& );
    void assembly_mat();
    void assembly_rhs();
    void solve();
    void print();

private:
   //! Mesh for the vessel network  (1D)
   getfem::mesh mesh_;   
   getfem::mesh_im mim_;           /// the integration methods
   getfem::mesh_fem mf_;         /// the main mesh_fem, for the pressure solution
   getfem::mesh_fem mf_coef_;      /// the mesh_fem to represent pde coefficients 
   std::vector<getfem::node> BCList_;
   std::vector<scalar_type> BC_value_;
   std::vector<int> local_idx_bc_;
   size_type nb_branches_; // number branches in the network
   

   
   getfem::ga_workspace workspace_; /// generic workspace for the solution
   
   std::vector<scalar_type> sol_, sol_old_, rhs_;    /// solution and old solution
   sparse_matrix_type K_;                         /// iteration matrix 	
   
   
   void configure_wp(const GetPot& df);	
   void genBC();

   void build_mesh(getfem::mesh& mesh, const GetPot& df);
}; // end network class
#endif
