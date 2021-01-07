#include<biot.hpp>

BIOT::BIOT(const GetPot& df): mim_(mesh_),
                              mf_u_(mesh_),
                              mf_p_(mesh_),
                              mf_coef_(mesh_),
                              tau_(1), mu_(1), bm_ (1),lambda_(1),
                              alpha_(1), permeability_(1),penalty_(1), beta_(1), fcoef_(1)
{
    std::cout<<"BIOT::BIOT build biot class"<<std::endl;
    bgeot::pgeometric_trans pgt = 
		bgeot::geometric_trans_descriptor(df("biot/mesh/mesh_type", "GT_PK(3,1)"));
    // Mesh generation
    N_ = pgt->dim(); 
    std::vector<size_type> nsubdiv(N_);
    std::fill(nsubdiv.begin(),nsubdiv.end(),df("biot/mesh/ndiv", 10));
    std::cout<< df("biot/mesh/noise", 0) <<std::endl;
    getfem::regular_unit_mesh(mesh_, nsubdiv, pgt, df("biot/mesh/noise", 0));
    ref_values.lref=df("biot/mesh/lref", 1.e+3);
    std::cout<<"end mesh generation"<<std::endl;

    // set  integration methods  
    getfem::pfem pf_u = getfem::fem_descriptor(df("biot/fem/fem_u", "FEM_PK(3,1)"));
    getfem::pfem pf_p = getfem::fem_descriptor(df("biot/fem/fem_p", "FEM_PK(3,1)"));
    getfem::pintegration_method ppi = getfem::int_method_descriptor(df("biot/fem/integr_method", "IM_TETRAHEDRON(6)"));
    //getfem::pintegration_method simp_ppi = getfem::int_method_descriptor(p_des.SIMPLEX_INTEGRATION);
    mim_.set_integration_method(mesh_.convex_index(), ppi);
    mf_u_.set_qdim(N_);
    mf_u_.set_finite_element(mesh_.convex_index(), pf_u); // finite element for displacement
    mf_p_.set_finite_element(mesh_.convex_index(), pf_p); //  finite element for pressure
    mf_coef_.set_finite_element(mesh_.convex_index(),getfem::classical_fem(pgt,0));// p0 for coefficient
    // ===boundary condition labelling
    gen_bc();
    // init matricies
    std::cout<<"init vectors"<<std::endl;
    getfem::size_type nb_dof_u = mf_u_.nb_dof();
    getfem::size_type nb_dof_p = mf_p_.nb_dof();
    gmm::resize(rhs_,nb_dof_u + nb_dof_p);
    gmm::resize(rhs_p_, nb_dof_p); 
    gmm::resize(rhs_u_, nb_dof_u);  // rhs fixed stress pressure displacement
    gmm::resize(U_, nb_dof_u); gmm::clear(U_);
    gmm::resize(U_old_, nb_dof_u); gmm::clear(U_old_);
    gmm::resize(P_, nb_dof_p); gmm::resize(P_old_, nb_dof_p);  
    std::fill(P_.begin(), P_.end(), 0); gmm::copy(P_,P_old_);
    gmm::resize(UP_,nb_dof_u + nb_dof_p); gmm::clear(UP_);
    std::cout<<"init matricies"<<std::endl;
    gmm::resize(Ku_, nb_dof_u , nb_dof_u ); gmm::clear(Ku_);
    gmm::resize(Kp_, nb_dof_p , nb_dof_p ); gmm::clear(Kp_);
    gmm::resize(K_,nb_dof_u + nb_dof_p ,nb_dof_u + nb_dof_p ); gmm::clear(K_);
    std::cout<<"number of dof "<< nb_dof_u<< " and "<<  nb_dof_p<<std::endl;
    std::cout <<"add variables to workspace"<<std::endl;
    workspace_.add_fem_variable("u", mf_u_, gmm::sub_interval(0, nb_dof_u), U_);
    workspace_.add_fem_variable("u_old", mf_u_, gmm::sub_interval(0, nb_dof_u), U_old_);
    workspace_.add_fem_variable("p", mf_p_, gmm::sub_interval(nb_dof_u, nb_dof_p), P_);
    workspace_.add_fem_variable("p_old", mf_p_, gmm::sub_interval(nb_dof_u, nb_dof_p), P_old_);
    configure_wp(df);
    
}/// end constructor


// ===========================================
// method for generation of bcs zones
// ===========================================
void BIOT::gen_bc()
{
   std::cout << "BIOT::gen_bc()"<< std::endl;
   getfem::mesh_region border_faces;
   getfem::outer_faces_of_mesh(mesh_, border_faces);
   for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
      assert(i.is_face());
      base_node un = mesh_.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);
      if ((un[N_-1] ) > 9.0E-1) { 
  	mesh_.region(TOP).add(i.cv(), i.f()); //all top iced
      } else if ((un[N_-1] ) < -9.0E-1) {  //the bottom surface is the most sharp
	mesh_.region(BOTTOM).add(i.cv(), i.f());
      } else if ((un[N_-2] ) < -1.0E-1) {
	mesh_.region(LEFT).add(i.cv(), i.f());
      } else if ((un[N_-2] ) > 1.0E-1) {
	mesh_.region(RIGHT).add(i.cv(), i.f());
      }	else if(N_=3){if ((un[N_-3] ) < -1.0E-1) {
			  mesh_.region(LEFTX).add(i.cv(), i.f());
			} else if ((un[N_-3] ) > 1.0E-1) {
			  mesh_.region(RIGHTX).add(i.cv(), i.f());
			}
		      }// end if n=3
    } // end for
}
// end bc generations
// =======================================================

void BIOT::configure_wp(const GetPot& df){
   std::cout << "BIOT::configure_wp Configuring workspace " << std::endl;
   ref_values.pref_out = df("biot/material/prefout", 0.9 );
   ref_values.pref_in = df("biot/material/prefin", 1. ); 
   double dp=ref_values.pref_out -ref_values.pref_in;
   ref_values.kref =  df("biot/material/k",1.);
   ref_values.tref =  ( ref_values.lref * ref_values.lref)/ref_values.kref;
   tau_[0] = ( df("time/dt", 1.0 )  ); // dt into  the workspace
   workspace_.add_fixed_size_constant("tau", tau_);
   double poisson=df("biot/material/poisson", 0.3 );
   double E=df("biot/material/E", 1. );
   double mu_s = E/( 2 * ( 1 + poisson) ) ;
   double lambda_l= E*poisson/ ( ( 1+poisson ) * (1 - 2 * poisson)) ;
   std::cout<<"mu "<<mu_s<< " lambda "<<  lambda_l<< " poisson "<< poisson<<" E "<< E <<std::endl;
   mu_[0] = mu_s;workspace_.add_fixed_size_constant("mu", mu_);
   //---------------------------------------------------------
   bm_[0] =  df("biot/material/biot_modulus",1.);
   workspace_.add_fixed_size_constant("bm", bm_);
   //---------------------------------------------------------
   lambda_[0] = lambda_l;
   workspace_.add_fixed_size_constant("lambda", lambda_);
   
   beta_[0] = ref_values.lref ;
   workspace_.add_fixed_size_constant("lref", beta_);
   //---------------------------------------------------------
   alpha_[0] = df("biot/material/alpha",1.);
   std::cout<<"alpha "<< alpha_[0]<<std::endl; 
   workspace_.add_fixed_size_constant("alpha", alpha_);
   //---------------------------------------------------------
   permeability_[0] =  df("biot/material/k",1.);
   workspace_.add_fixed_size_constant("permeability", permeability_);
   //---------------------------------------------------------
   penalty_[0] = 1.e+12; // 1---10
   workspace_.add_fixed_size_constant("penalty", penalty_);
   std::cout<<"BIOT::configure_wp end of workspace configuration"<<std::endl;
   
   fcoef_[0] = (2200*0.8 + 1000*0.2 -1000) *ref_values.lref*ref_values.lref;
   workspace_.add_fixed_size_constant("fr", fcoef_);
}


void BIOT::assembly_mat(){
   std::cout<<"BIOT::assembly_mat start assembling K"<<std::endl;
// ------------------ expressions --------------------------
   workspace_.add_expression("2*mu*Sym(Grad_u):Grad_Test_u + lambda*Div_u*Div_Test_u - lref*alpha*p.Div_Test_u ", mim_);
   workspace_.add_expression( "+(1/bm)*p.Test_p/permeability/tau + Grad_p.Grad_Test_p+ alpha*Test_p*Div_u/permeability/tau", mim_);
   // workspace_.assembly(2);
   // gmm::copy(workspace_.assembled_matrix(), K_);
   // workspace_.clear_expressions();
// bcs
   workspace_.add_expression("penalty/element_size*p*Test_p/permeability", mim_, TOP);
   workspace_.add_expression("-Grad_p.Normal*Test_p - Grad_Test_p.Normal*p ", mim_, TOP); 	
   // workspace_.add_expression("penalty/element_size*p*Test_p", mim_, LEFT);
   // workspace_.add_expression("penalty/element_size*p*Test_p", mim_, RIGHT);// 1 is the region	
   // workspace_.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim_, LEFT); 	
   // workspace_.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p", mim_, RIGHT); 
   // workspace_.add_expression("penalty/element_size*p*Test_p", mim_, LEFTX);
   // workspace_.add_expression("penalty/element_size*p*Test_p", mim_, RIGHTX);// 1 is the region	
   if(N_==3){	 	
   //   workspace_.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p ", mim_, LEFTX); 	
   //   workspace_.add_expression("-permeability*Grad_p.Normal*Test_p - permeability*Grad_Test_p.Normal*p", mim_, RIGHTX);
     }	
   workspace_.add_expression("penalty*u.Test_u" , mim_, BOTTOM); //neumann disp
   workspace_.assembly(2);
   gmm::add(workspace_.assembled_matrix(), K_);
   workspace_.clear_expressions();
   std::cout<<"BIOT::assembly_mat End assembling K"<<std::endl;
}



void BIOT::assembly_rhs(){
   std::cout<<"BIOT::assembly_mat start assembling rhs"<<std::endl;
   gmm::clear(rhs_);
   //======= RHS =====================
   if(N_==2) workspace_.add_expression("0.*0.001*(2200*0.8 + 1000*0.2 -1000 )*[0,-1].Test_u", mim_);
   if(N_==3) workspace_.add_expression("1.*fr*[0,0,-1].Test_u", mim_);
   workspace_.add_expression("+[+0.0e-6].Test_p + (1/bm)*p_old.Test_p/permeability/tau + alpha*Test_p*Div_u_old/permeability/tau", mim_);
   workspace_.set_assembled_vector(rhs_);
   workspace_.assembly(1);
   workspace_.clear_expressions();
   std::cout<<"BIOT::assembly_mat end assembling rhs"<<std::endl;

   // all bcs are homogeneus no rhs terms
}


//=====================================
void BIOT::solve(){
   size_type restart = 50;
   gmm::identity_matrix PM; // no precond
   gmm::diagonal_precond<sparse_matrix_type> PR(K_);
   // gmm::mr_approx_inverse_precond<sparse_matrix_type> P(K, 10, 10E-17);
   // gmm::iteration iter(1.e-8);  // iteration object with the max residu
   //	iter.set_noisy(1);               // output of iterations (2: sub-iteration)
   //	iter.set_maxiter(1000); // maximum number of iterations
   //	gmm::gmres(K, UP, B, PR, restart, iter);
   scalar_type cond;
   gmm::SuperLU_solve(K_, UP_ , rhs_, cond);
   std::cout<<"condition number "<< cond<< std::endl;
   getfem::size_type nb_dof_u = mf_u_.nb_dof();
   getfem::size_type nb_dof_p = mf_p_.nb_dof();
   gmm::copy(gmm::sub_vector(UP_,gmm::sub_interval(0, nb_dof_u)), U_);				
   gmm::copy(gmm::sub_vector(UP_,gmm::sub_interval(nb_dof_u , nb_dof_p)), P_);	
   gmm::copy(U_,U_old_);gmm::copy(P_,P_old_);
}
//====================================
void BIOT::setsol(std::vector<scalar_type>& sol){

   getfem::size_type nb_dof_u = mf_u_.nb_dof();
   getfem::size_type nb_dof_p = mf_p_.nb_dof();
   gmm::copy(gmm::sub_vector(sol,gmm::sub_interval(0, nb_dof_u)), U_);				
   gmm::copy(gmm::sub_vector(sol,gmm::sub_interval(nb_dof_u , nb_dof_p)), P_);	
   gmm::copy(U_,U_old_);gmm::copy(P_,P_old_);
}

//=====================================
void BIOT::print(int iter){
std::cout<< "BIOT::print() "<< "biot.vtk"<<std::endl;
getfem::vtk_export exp("biot."+std::to_string(iter)+".vtk");
exp.exporting(mf_u_);  	exp.write_point_data(mf_u_, U_, "u"); 
exp.write_point_data(mf_p_, P_, "p"); 
}
