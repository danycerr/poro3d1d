#ifndef BIOTHS
#define BIOTHS

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

//namespaces and type defs
using bgeot::base_small_vector; /* special class for small (dim<16) vectors */
using bgeot::base_node; /* geometrical nodes (derived from base_small_vector)*/
using bgeot::scalar_type; /* = double */
using bgeot::size_type;   /* = unsigned long */
using bgeot::dim_type; 

typedef gmm::rsvector<scalar_type> sparse_vector_type;
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;
typedef gmm::col_matrix<sparse_vector_type> col_sparse_matrix_type;
typedef std::vector<scalar_type> plain_vector;

class BIOT
{
public:
    BIOT(const GetPot& );
    void assembly_mat();
    void assembly_rhs();
    void solve();
    void print(int iter=0);

    inline int ndof(){return mf_u_.nb_dof() + mf_p_.nb_dof();} // get dof number;

    inline int ndof_p(){return mf_p_.nb_dof();} // get dof number;
    inline int ndof_u(){return mf_u_.nb_dof();} // get dof number;

    inline sparse_matrix_type& get_iter_mat() {return K_;}
    inline std::vector<scalar_type>& get_iter_rhs() {return rhs_;}
    inline getfem::mesh_fem& get_fem_p() {return mf_p_;}
    inline getfem::mesh_fem& get_fem_u() {return mf_u_;}
    void   setsol(std::vector<scalar_type>& sol); // copy the solution
private:
   getfem::mesh mesh_;
   getfem::mesh_im mim_;           /// the integration methods
   getfem::mesh_fem mf_u_;         /// the main mesh_fem, for the displacement solution
   getfem::mesh_fem mf_p_;         /// the main mesh_fem, for the pressure solution
   getfem::mesh_fem mf_coef_;      /// the mesh_fem to represent pde coefficients 
   getfem::ga_workspace workspace_; /// generic workspace for the solution
   enum { BOTTOM = 1, TOP = 2, LEFT = 3, RIGHT =4, LEFTX = 5, RIGHTX =6}; // descriptor for zones

   int N_;                     /// problem dimension
   // material properties
   std::vector<scalar_type> tau_, mu_, bm_ ,lambda_, alpha_, permeability_,penalty_;

   std::vector<scalar_type> U_, U_old_, P_, P_old_, UP_, rhs_u_,rhs_p_, rhs_;    /// diplacement,  pressure, pressure_old
   sparse_matrix_type Kp_, Ku_, K_;                         /// iteration matrix 	

   void gen_bc(void);    /// create zones for boundary conditions
   void configure_wp(const GetPot& df);	

};
#endif
