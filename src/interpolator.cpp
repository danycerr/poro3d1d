#include"interpolator.hpp"

INTERPOLATOR::INTERPOLATOR(const GetPot& df, getfem::mesh_fem& mf_1
                                           , getfem::mesh_fem& mf_2, getfem::mesh_im& mim_2, std::vector<scalar_type>& radius  ):
                           mf_1_(&mf_1), mf_2_(&mf_2), mim_2_(&mim_2),radius_(&radius),
                           NInt_(df("interpolator/NInt", 50))

{
   std::cout<<"INTERPOLATOR::INTERPOLATOR Building interpolator class " << std::endl;
   std::cout<<"mf_1 nb dof "<<mf_1_->nb_dof()<<std::endl;
   
   // 2 to 1 (1d to 3d ) operators
   Mbar_.resize(mf_2_->nb_dof(), mf_1_->nb_dof());gmm::clear(Mbar_);
   Mlin_.resize(mf_2_->nb_dof(), mf_1_->nb_dof());gmm::clear(Mlin_);


}


void INTERPOLATOR::build_operator(){
   std::cout<<"INTERPOLATOR::build_operator Building operators 3D1D" << std::endl;

   gmm::clear(Mbar_); gmm::clear(Mlin_);
   const scalar_type Pi = 2*acos(0.0);
   size_type nb_dof_t     = mf_1_->nb_dof();
   size_type nb_dof_v     = mf_2_->nb_dof();
   // Linear interpolation map mf_1 --> mf_2
   getfem::interpolation(*mf_1_, *mf_2_, Mlin_);
   // Aux vectors for local interpolation
   std::vector<scalar_type> Pbari(NInt_); 
   std::vector<scalar_type> Pt(nb_dof_t); 
   size_type counter = 0;
   for (size_type i = 0; i < nb_dof_v; i++){
      counter++;
      if (counter*100 >= nb_dof_v) {
         counter = 0; 
	 std::cout << "*"; std::cout.flush();
        }
      // We need the following interpolation tool <getfem_interpolation.h>      
      getfem::mesh_trans_inv mti(mf_1_->linked_mesh());
   std::cout<<"End linear mesh_trans_inv"<< std::endl;
      // Build the list of point on the i-th circle:
      // ... first consider an orthonormal system v0, v1, v2:
      base_node v0;
      if (i==0) {
      v0 = mf_2_->point_of_basic_dof(i+1) - mf_2_->point_of_basic_dof(i);			
		} else {
      v0 = mf_2_->point_of_basic_dof(i) - mf_2_->point_of_basic_dof(i-1);
		}
      base_node v1(0.0, -v0[2], v0[1]);
      base_node v2(v0[1]*v0[1] +v0[2]*v0[2], -v0[0]*v0[1], -v0[0]*v0[2]);
      if (gmm::vect_norm2(v2) < 1.0e-8 * gmm::vect_norm2(v0)) {
         v1[0] = -v0[1]; v1[1] = v0[0]; v1[2] = 0.0;
         v2[0] = -v0[0]*v0[2]; v2[1] = -v0[1]*v0[2]; v2[2] = v0[0]*v0[0] +v0[1]*v0[1];
         }
      v1 = v1 / gmm::vect_norm2(v1);
      v2 = v2 / gmm::vect_norm2(v2);
      // ... then parametrize the circle:
      for (size_type j = 0; j < NInt_; j++){ 	
         mti.add_point( mf_2_->point_of_basic_dof(i) + 
            ((*radius_)[i])*(cos(2*Pi*j/NInt_)*v1 + sin(2*Pi*j/NInt_)*v2) );
         }
		// Get the local interpolation matrix Mbari
		sparse_matrix_type Mbari(NInt_, nb_dof_t); gmm::clear(Mbari);
		getfem::interpolation(*mf_1_, mti, Pt, Pbari, Mbari, 1);
		scalar_type sum_row = 0.0;
		for (size_type j=0; j < NInt_; ++j) {
			typename gmm::linalg_traits<sparse_matrix_type>::const_sub_row_type 
				row = mat_const_row(Mbari,j);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				it_nz = vect_const_begin(row);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				ite_nz = vect_const_end(row);
			for (; it_nz != ite_nz ; ++it_nz) {
				Mbar_(i, it_nz.index()) += (*it_nz);
				sum_row += (*it_nz);
			}
		}
		gmm::cout<<"=====dof_v"<<i<<"========sum_row"<<sum_row <<gmm::endl;
		typename gmm::linalg_traits<sparse_matrix_type>::sub_row_type 
			row = mat_row(Mbar_,i);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			it_nz = vect_begin(row);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			ite_nz = vect_end(row);
		for (; it_nz != ite_nz ; ++it_nz) {
			(*it_nz)/=sum_row;
		}

	} // end of outer for loop 
	std::cout << std::endl;

}



void INTERPOLATOR::extend2vector(getfem::mesh_fem& mf_1_vec, getfem::mesh_fem& mf_2_vec){
    
    std::cout<<"INTERPOLATOR::extend2vector extd interpolation 2 vectorial field"<<std::endl;
    mf_1_vec_ = &mf_1_vec;  // bulk
    mf_2_vec_ = &mf_2_vec; // network 
    
    Mbar_vec_.resize(mf_2_vec_->nb_dof(), mf_1_vec_->nb_dof());gmm::clear(Mbar_vec_);
    std::cout<<"size of Mbar_vec_ "<<mf_2_vec_->nb_dof()<<" x "<<  mf_1_vec_->nb_dof()<<std::endl;
    
    size_type nb_dof_t     = mf_1_->nb_dof();
    size_type nb_dof_v     = mf_2_->nb_dof();
    
    for (int idof=0; idof<nb_dof_v; idof++ ){
    bgeot::base_node pt_1d(mf_2_->point_of_basic_dof(idof));
        // evaluate corresponding 1d vec dof
    std::vector<int> vec_1d_dof;
    for(int idof_v=0; idof_v<mf_2_vec_->nb_dof(); idof_v++){
            bgeot::base_node pt_vec_1d(mf_2_vec_->point_of_basic_dof(idof_v));
            double distance = 0;
            for (int idim=0; idim<3; idim++) distance+=fabs(pt_1d[idim]-pt_vec_1d[idim]);
            if(distance<1.e-16){
                    vec_1d_dof.push_back(idof_v);
            }
    }

        
        
        
        
        std::vector<int> nonzero;
        for (int icol=0; icol<nb_dof_t; icol++){
            if(fabs(Mbar_(idof, icol))>1.e-16){
            nonzero.push_back(icol);
            }
        }
        // std::cout <<"dof "<<idof <<"non zero in "<< nonzero<<std::endl;
        for(int ibulkdof=0; ibulkdof<nonzero.size(); ibulkdof++){
                int bdof=nonzero[ibulkdof]; // target dof
                bgeot::base_node pt(mf_1_->point_of_basic_dof(bdof));
                std::vector<int> vec_dof;
                for(int ivecbulkdof=0; ivecbulkdof<mf_1_vec_->nb_dof(); ivecbulkdof++){
                    bgeot::base_node pt_vec(mf_1_vec_->point_of_basic_dof(ivecbulkdof));
                    double distance = 0;
                    for (int idim=0; idim<3; idim++) distance+=fabs(pt[idim]-pt_vec[idim]);
                    if(distance<1.e-16){
                        vec_dof.push_back(ivecbulkdof);
                    }
                }
            for(int idof_v_1d=0; idof_v_1d < vec_1d_dof.size(); idof_v_1d++ ){
                int ii= vec_1d_dof[idof_v_1d];
                int jj= vec_dof[idof_v_1d];
                Mbar_vec_(ii,jj) = Mbar_(idof,bdof);
            }
                
        }
    }
    
    std::cout<<" INTERPOLATOR::extend2vector end"<<std::endl;
}
