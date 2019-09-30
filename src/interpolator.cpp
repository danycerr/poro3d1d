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
