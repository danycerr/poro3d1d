#include<network.hpp>

NETWORK::NETWORK (const GetPot& df): mim_(mesh_),
                                     mf_(mesh_),
                                     mf_coef_(mesh_),
                                     BC_value_(1)
{
  nb_vertices_.resize(0); nb_vertices_.clear();
  BCList_.resize(0);
  local_idx_bc_.resize(0);
  build_mesh(mesh_, df );
  std::cout<<"NETWORK::NETWORK  Configuring fem"<<std::endl;
  bgeot::pgeometric_trans pgt = 
		bgeot::geometric_trans_descriptor(df("1d/mesh/mesh_type", "GT_PK(1,1)"));
  getfem::pfem pf = getfem::fem_descriptor(df("1d/fem/fem", "FEM_PK(1,1)"));
  getfem::pintegration_method ppi = getfem::int_method_descriptor(df("biot/fem/integr_method", "IM_GAUSS1D(6)"));
  mim_.set_integration_method(mesh_.convex_index(), ppi);
  mf_.set_finite_element(mesh_.convex_index(), pf); // finite element for displacement
  mf_coef_.set_finite_element(mesh_.convex_index(),getfem::classical_fem(pgt,0));// p0 for coefficient
  std::cout<<"end NETWORK::NETWORK  Configuring fem"<<std::endl;
  genBC();
  // init matricies
  std::cout<<"NETWORK::init vectors"<<std::endl;
  getfem::size_type nb_dof = mf_.nb_dof();
  gmm::resize(rhs_,nb_dof);
  gmm::resize(sol_, nb_dof); gmm::clear(sol_);
  gmm::resize(sol_old_, nb_dof); gmm::clear(sol_old_);  
  std::fill(sol_.begin(), sol_.end(), 0); gmm::copy(sol_,sol_old_);
  gmm::resize(Y_, mf_coef_.nb_dof());std::fill(Y_.begin(), Y_.end(), 1.);// coefficient mass matrix
  std::cout<<"init matricies"<<std::endl;
  gmm::resize(K_, nb_dof, nb_dof ); gmm::clear(K_);
  std::cout<<"number of dof "<< nb_dof<<std::endl;
  std::cout <<"add variables to workspace"<<std::endl;
  workspace_.add_fem_variable("p", mf_, gmm::sub_interval(0, nb_dof), sol_);
  workspace_.add_fem_variable("p_old", mf_, gmm::sub_interval(0, nb_dof), sol_old_);
  configure_wp(df);
    
  
}//end constructor

void NETWORK::genBC(){
   getfem::mesh_region border_faces;
   getfem::outer_faces_of_mesh(mesh_, border_faces);
   for(int ibc=0; ibc<BCList_.size();ibc++){
      for (getfem::mr_visitor i(border_faces); !i.finished(); ++i) {
        for (int ipt =0;ipt<mesh_.points_of_convex(i.cv()).size();ipt++ ){ 
        double distance=0;
        for (int dim=0; dim <3; dim++) 
         {
          distance+= fabs( (mesh_.points_of_convex(i.cv())[ipt])[dim] - mesh_.points()[BCList_[ibc].idx][dim]);
         }
        // std::cout << "i= "<<ibc <<" local point" << ipt << " distance "<< distance<<std::endl;
        if(distance<1.e-8) {
            mesh_.region(100+BCList_[ibc].idx).add(i.cv(), i.f());
            local_idx_bc_.push_back(ipt);
            std::cout << "adding" <<  i.cv() << " to region "<< 100+BCList_[ibc].idx<<std::endl;
           }
       }
      }
   }

}

//=======================================================================
void  NETWORK::build_mesh(getfem::mesh& mesh, const GetPot& df){
   std::cout<< "NETWORK::build_mesh"<<std::endl;
   bgeot::pgeometric_trans pgt = 
		bgeot::geometric_trans_descriptor(df("1d/mesh/mesh_type", "GT_PK(1,1)"));
  
   mesh.clear();
   std::ifstream ifs(df("1d/mesh/input_fn", "segment.pts")); 
   // ===import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV);
   size_type Nb = 0; // nb of branches
   ifs.precision(16); ifs.seekg(0); ifs.clear();
   GMM_ASSERT1(bgeot::read_until(ifs, "BEGIN_LIST"), 
		"NETWORK::build_mesh This seems not to be a data file");

   size_type globalBoundaries = 0;
   while (bgeot::read_until(ifs, "BEGIN_ARC")) 
      {

         Nb++;nb_vertices_.emplace_back(0);
         std::vector<base_node> lpoints; 
         dal::dynamic_array<scalar_type> tmpv;
         std::string tmp, BCtype, value;
	 bool thend = false; 
	 size_type bcflag = 0;
	 size_type bcintI = 0, bcintF = 0;
	 getfem::node BCA, BCB;
         // Read an arc from data file and write to lpoints
	 while (!thend) {
	    bgeot::get_token(ifs, tmp, 1023);
	    if (tmp.compare("END_ARC") == 0) { 
               thend = true;
               }
            else if (ifs.eof()) {
               GMM_ASSERT1(0, "Unexpected end of stream");
               }
            else if (tmp.compare("BC") == 0) {
               bcflag++;bgeot::get_token(ifs, BCtype, 4);
               if (BCtype.compare("DIR") == 0) {
                  bgeot::get_token(ifs, value, 1023);
                  if (bcflag == 1) {
                     BCA.label = BCtype;BCA.value = stof(value); 
                     globalBoundaries++;
                     }
                   else if (bcflag == 2) {
                      BCB.label = BCtype; BCB.value = stof(value); 
                      globalBoundaries++;
                      }
                   else GMM_ASSERT1(0, "More than 2 BC found on this arc!");
               } // BCtype.compare DIR
               else if (BCtype.compare("MIX") == 0) {
               bgeot::get_token(ifs, value, 1023);
               if (bcflag == 1) {
                  BCA.label = BCtype; BCA.value = stof(value); 
                  globalBoundaries++;
                  }
               else if (bcflag == 2) {
                  BCB.label = BCtype; BCB.value = stof(value); 
                  globalBoundaries++;
                  }
               } // end bc type compare mix
               else if (BCtype.compare("INT") == 0) {
                  if (bcflag == 1) {
                     bcintI++;BCA.label = BCtype; 
		     }
                  else if (bcflag == 2) {
                     bcintF++;BCB.label = BCtype; 
		     }
		  else
		     GMM_ASSERT1(0, "More than 2 BC found on this arc!");
                  } // end of BCtype.compare("INT")
               else
                  GMM_ASSERT1(0, "Unknown Boundary condition");	  
               } // end of tmp.compare("BC") "BC" compare 
          else if (tmp.size() == 0) {
             GMM_ASSERT1(0, "Syntax error in file, at token '" 
                << tmp << "', pos=" << std::streamoff(ifs.tellg()));
          }

          else { // "point" case 
          nb_vertices_[Nb-1]++;
          int d = 0;
          while ( (isdigit(tmp[0]) != 0) || tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.'){ 
             tmpv[d++] = stof(tmp); 
             bgeot::get_token(ifs, tmp, 1023); 
             }
          if (d != 4) GMM_ASSERT1(0, "Points must have 3 coordinates - number of coordinates:" << d);
          base_node tmpn(tmpv[1], tmpv[2], tmpv[3]);
          lpoints.push_back(tmpn);
          if (tmp.compare("END_ARC") == 0) { thend = true; nb_vertices_[Nb-1]--; }
       } 				
    }  // end of inner while
    // Insert the arc into the 1D mesh and build a new region for the corresponding branch
    // Check validity of branch region
    GMM_ASSERT1(mesh.has_region(Nb-1)==0, "Overload in meshv region assembling!");
    for (size_type i=0; i<lpoints.size()-1; ++i ) {
       std::vector<size_type> ind(2);
       size_type ii = (i==0) ? 0 : i+1;
       size_type jj;
       if (ii == lpoints.size()-1) jj = 1;
       else if (ii == 0) jj = 2;
       else jj = ii+1;
       ind[0] = mesh.add_point(lpoints[ii]);
       ind[1] = mesh.add_point(lpoints[jj]);
       size_type cv;
       cv = mesh.add_convex(pgt, ind.begin());
       // Build branch regions
       mesh.region(Nb-1).add(cv);
       if ((bcflag>0) && (ii==0)&& (bcintI==0)) {
          BCA.idx = ind[0];BCList_.push_back(BCA);
          }
       if ((bcflag>1) && (jj==1) && (bcintF==0)) {
          BCB.idx = ind[1];BCList_.push_back(BCB);
          }
       } // end of inner for 
   } // end of outer while 	


   nb_branches_ = nb_vertices_.size();
   std::cout<< "Number of branches "<< nb_branches_<<std::endl;
   ifs.close();
   std::cout <<"radius import only constant"<<std::endl;
   radius_.resize(nb_branches_);
   for (int i=0; i< radius_.size(); i++) radius_[i]=df("1d/radius", 0.1);
   #ifdef MESH_INFO
   std::cout << "Mesh "<<std::endl;
   for (dal::bv_visitor i(mesh.points_index()); !i.finished(); ++i)
      std::cout << "Point of index " << i << " of the mesh: " << mesh.points()[i] << std::endl;
   #endif
   std::cout<<"end NETWORK::build_mesh"<<std::endl;
}

//=====================================================
void NETWORK::configure_wp(const GetPot& df){
 BC_value_[0] =  0.;
 workspace_.add_fixed_size_constant("bc_value", BC_value_);
}

//===========================================
void NETWORK::assembly_mat(){
   std::cout<<"NETWORK::assembly_mat start assembling K"<<std::endl;
   // ------------------ expressions --------------------------
   workspace_.add_expression( "+p.Test_p + Grad_p.Grad_Test_p", mim_);
   workspace_.assembly(2);
   gmm::copy(workspace_.assembled_matrix(), K_);
   workspace_.clear_expressions();
  
  std::cout<< BCList_<<std::endl;
  for(int i=0; i<BCList_.size();i++){
     if ((BCList_[i].label).compare("DIR")==0)
        {  int counter=0;
           for (dal::bv_visitor ipt(mf_.dof_on_region(100 + BCList_[i].idx)); !ipt.finished(); ++ipt){ 
           std::cout << counter << " " << ipt<<std::endl;
           if (local_idx_bc_[i]==counter) {
                  K_(ipt,ipt)+=1.e+8;std::cout << "adding bc on "<< ipt <<" at "<< local_idx_bc_[i] <<std::endl;
               }
           counter++;
           }
       }
   }
}

//===========================================
void NETWORK::assembly_rhs(){
   std::cout<<"NETWORK::assembly_mat start assembling rhs"<<std::endl;
   // ------------------ expressions --------------------------
   workspace_.add_expression( "+1*Test_p + p_old.Test_p", mim_);
   workspace_.set_assembled_vector(rhs_);
   workspace_.assembly(1);
   workspace_.clear_expressions();

  std::cout<< BCList_<<std::endl;
  for(int i=0; i<BCList_.size();i++){
     if ((BCList_[i].label).compare("DIR")==0)
        {  int counter=0;
           for (dal::bv_visitor ipt(mf_.dof_on_region(100 + BCList_[i].idx)); !ipt.finished(); ++ipt){ 
           std::cout << counter << " " << ipt<<std::endl;
           if (local_idx_bc_[i]==counter) {
                  rhs_[ipt]+=1.e+8*BCList_[i].value;std::cout << "adding rhs bc on "<< ipt <<" at "<< local_idx_bc_[i] <<std::endl;
               }
           counter++;
           }
       }
   }
}
  


//=====================================
void NETWORK::solve(){
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
   gmm::copy(sol_,sol_old_);
}

//=====================================
void NETWORK::print(){
std::cout<< "Network::print() "<< "nework.vtk"<<std::endl;
//getfem::vtk_export exp(p_des.datafilename + "." +  std::to_string(time) + ".vtk");
getfem::vtk_export exp("network.vtk");
 exp.exporting(mf_);  	exp.write_point_data(mf_, sol_, "pod"); 
// exp.write_point_data(mf_p_, P_, "p"); 
}

