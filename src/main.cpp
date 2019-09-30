#include <iostream>
#include"biot.hpp"
#include"network.hpp"
#include"problem3d1d.hpp"

int main(int argc, char** argv)
{
   GetPot   cl(argc, argv);
   const std::string data_file_name = cl.follow("data", 2, "-f",
            "--file");

   GetPot dataFile(data_file_name.data());
   
   
   PROBLEM_3d1d bd(dataFile);
   bd.assembly_mat();
   bd.assembly_rhs();
   bd.assembly_exchange_mat();
   bd.solve();
   bd.print();

  if(0) { 
   NETWORK oned_problem(dataFile);
   oned_problem.assembly_mat();
   oned_problem.assembly_rhs();
   oned_problem.solve();
   oned_problem.print();
}

if(0) {
   BIOT biot_problem(dataFile);
   biot_problem.assembly_mat();
   biot_problem.assembly_rhs();
   biot_problem.solve();
   biot_problem.print();
      }   
}
