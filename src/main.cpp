#include <iostream>
#include"biot.hpp"
#include"network.hpp"

int main(int argc, char** argv)
{
   GetPot   cl(argc, argv);
   const std::string data_file_name = cl.follow("data", 2, "-f",
            "--file");

   GetPot dataFile(data_file_name.data());
   
   // init 3d problem
   BIOT biot_problem(dataFile);
   
   NETWORK oned_problem(dataFile);
   
   oned_problem.assembly_mat();
   oned_problem.assembly_rhs();
   oned_problem.solve();
   oned_problem.print();


if(0) {
   biot_problem.assembly_mat();
   biot_problem.assembly_rhs();
   biot_problem.solve();
   biot_problem.print();
      }   
}
