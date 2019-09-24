#include "problem3d1d.hpp"

PROBLEM_3d1d::PROBLEM_3d1d(const GetPot& df):biot_problem(df),
                                             oned_problem(df)
{
   std::cout<<"PROBLEM_3d1d::PROBLEM_3d1d init vectors"<<std::endl;
   nb_dof_biot_ = biot_problem.ndof();
   nb_dof_oned_ = oned_problem.ndof();
   std::cout<<"nb dof biot "<< nb_dof_biot_<< " nb dof 1d "<< nb_dof_oned_<<std::endl;
   gmm::resize(rhs_,nb_dof_biot_ + nb_dof_oned_);
   gmm::resize(sol_, nb_dof_biot_ + nb_dof_oned_); gmm::clear(sol_);
   gmm::resize(sol_old_, nb_dof_biot_ + nb_dof_oned_); gmm::clear(sol_old_);
   std::cout<<"PROBLEM_3d1d:: init matricies"<<std::endl;
   gmm::resize(K_, nb_dof_biot_ + nb_dof_oned_, nb_dof_biot_ + nb_dof_oned_ ); gmm::clear(K_);
   
}

//====================================
void PROBLEM_3d1d::assembly_mat(){
   biot_problem.assembly_mat();
   gmm::copy(biot_problem.get_iter_mat(),gmm::sub_matrix(K_, 
                                             gmm::sub_interval(0, nb_dof_biot_), 
                                             gmm::sub_interval(0, nb_dof_biot_)));
   oned_problem.assembly_mat();
   gmm::copy(oned_problem.get_iter_mat(),gmm::sub_matrix(K_, 
                                             gmm::sub_interval(nb_dof_biot_, nb_dof_oned_), 
                                             gmm::sub_interval(nb_dof_biot_, nb_dof_oned_)));
   
}

//====================================
void PROBLEM_3d1d::assembly_rhs(){
   biot_problem.assembly_rhs();
   gmm::copy(biot_problem.get_iter_rhs(), gmm::sub_vector(rhs_,gmm::sub_interval(0, nb_dof_biot_)));
   oned_problem.assembly_rhs();
   gmm::copy(oned_problem.get_iter_rhs(), gmm::sub_vector(rhs_,gmm::sub_interval(nb_dof_biot_,nb_dof_oned_)));
}

//=====================================
void PROBLEM_3d1d::solve(){
   size_type restart = 50;
   gmm::identity_matrix PM; // no precond
   gmm::diagonal_precond<sparse_matrix_type> PR(K_);
   // gmm::mr_approx_inverse_precond<sparse_matrix_type> P(K, 10, 10E-17);
   // gmm::iteration iter(1.e-8);  // iteration object with the max residu
   //	iter.set_noisy(1);               // output of iterations (2: sub-iteration)
   //	iter.set_maxiter(1000); // maximum number of iterations
   //	gmm::gmres(K, UP, B, PR, restart, iter);
   scalar_type cond;
   gmm::SuperLU_solve(K_, sol_ , rhs_, cond);
   std::cout<<"condition number "<< cond<< std::endl;
   std::vector<scalar_type> buf(nb_dof_biot_);
   gmm::copy(gmm::sub_vector(sol_,gmm::sub_interval(0 , nb_dof_biot_)), buf);
   biot_problem.setsol(buf);
   buf.resize(nb_dof_oned_);
   gmm::copy(gmm::sub_vector(sol_,gmm::sub_interval(nb_dof_biot_ , nb_dof_oned_)), buf);
   oned_problem.setsol(buf);
   	
   // gmm::copy(gmm::sub_vector(UP_,gmm::sub_interval(nb_dof_u , nb_dof_p)), P_);	
   // gmm::copy(U_,U_old_);gmm::copy(P_,P_old_);
}
//====================================
void PROBLEM_3d1d::print(){
   biot_problem.print();
   oned_problem.print();
}
