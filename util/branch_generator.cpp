#include<iostream>
#include <vector>
#include <iterator> // for ostream_iterator
#include <fstream>
#include <math.h>

class branch{
    std::vector<double> p1={0,0,0};
    std::vector<double> p2={0,0,0};
    int ndiv=0;
    void generate(void);
    void print(void);    
public:   
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    double val1=0;
    double val2=0;
    branch (std::vector<double>, std::vector<double>,int, double, double );
    };

branch::branch (std::vector<double> up1, std::vector<double> up2,int div,
                double valpt1,double valpt2): 
            p1(up1),p2(up2),ndiv(div), val1(valpt1), val2(valpt2)
    {
     generate();
    }
    
void branch::generate(){
    std::cout<< "Branch generation"<<std::endl;
    double dx= (p2[0]-p1[0])/ndiv;
    double dy= (p2[1]-p1[1])/ndiv;
    double dz= (p2[2]-p1[2])/ndiv;
    x.resize(ndiv+1);
    y.resize(ndiv+1);
    z.resize(ndiv+1);
    for(int i=0;i<ndiv+1;i++){
      x[i]=dx*i + p1[0];y[i]=dy*i+ p1[1];z[i]=dz*i+ p1[2];
      }
    }

void print(std::vector<branch> branches){
std::ofstream out;
out.open("./horizontal_single.pts");
out<<"BEGIN_LIST" <<std::endl;
for(int b=0;b<branches.size();b++){
out<<"BEGIN_ARC" <<std::endl;
out<<"BC DIR  " << branches[b].val1 <<std::endl;
out<<"BC MIX  " << branches[b].val2 <<std::endl;
out<<"\t" << 1<<"\t"<<branches[b].x[0]<< "\t" << branches[b].y[0]<< "\t" << branches[b].z[0]<<"\t"<<"begin"<<std::endl;
out<<"\t"  << 1<<"\t"
          << branches[b].x[branches[b].x.size()-1] << "\t"  
          << branches[b].y[branches[b].x.size()-1] << "\t" 
          << branches[b].z[branches[b].x.size()-1] << "\t"
          << "end" <<std::endl;
for (int i=1; i< branches[b].x.size()-1;i++)
            out<<"\t" << 1<<"\t"
                      << branches[b].x[i]<< "\t" 
                      << branches[b].y[i]<< "\t" 
                      << branches[b].z[i]<< "\t" 
                     << "point" <<std::endl;
out<<"END_ARC"<<std::endl;
}
out<<"END_LIST"<<std::endl;
out.close();
}

int main (){
    std::cout<<"Generator of branches"<<std::endl;
    
//     std::vector<double> p1={0.053,0.053,0.}; double pt1d=0.;
//     std::vector<double> p2={0.053,0.053,1.}; double pt2d=0.;
//     branch b1(p1,p2,1000, pt1d,pt2d);
    std::vector<branch> branches;
//     branches.push_back(b1);
    
  
//     radial wels
//     int ntheta=20;
//     double dtheta=360/ntheta;
//     double r=0.3, xc=p1[0], yc=p1[1], r2=0.2;    
//     for(int k=0;k<ntheta;k++){	
//      // std::vector<double> p3={r2*cos(k*dtheta)+xc,r2*sin(k*dtheta)+yc,0};double pt3d=0.77;
//      // std::vector<double> p4={r*cos(k*dtheta)+xc,r*sin(k*dtheta)+yc,0.1};double pt4d=0.;
// 
//      std::vector<double> p3={(double)(rand() % 100)/100.,0.66*((double)(rand() % 100)/100.) ,0};double pt3d=0.77;
//      std::vector<double> p4={(double)(rand() % 100)/100.,0.66*((double)(rand() % 100)/100.) ,0.1};double pt4d=0.;
// 
//     branch b2(p3,p4,10, pt3d,pt4d);
//     branches.push_back(b2);
// }
// //      assembly grid
//     int nrow=3; int ncol=3;
//     for (int ib=0; ib<3; ib++){
//         double dx=0.1/(nrow+1.);
//         double dy=0.1/(ncol+1.);
//         double x=0.1*ib; 
//         for(int ix=0;ix<nrow;ix++){
//         x+=dx;double y=0.;
//         for(int iy=0;iy<ncol;iy++){
//             y+=dy;
//             std::vector<double> p_low={x,y,0}; double p_low_v=30.;
//             std::vector<double> p_up={x,y,1};  double p_up_v=0.;
//             branch rod(p_low,p_up,31, p_low_v,p_up_v);
//             branches.push_back(rod);
//         }
//         }
//     }
// end assembly grid
    
    

    
//         // nacie esagonal grid
//     // central rod
//     double eps=1.e-9;double h=2.;int ndiv=10;
//     double p=8.39972*1.e-3;
//     //central rod
//     std::vector<double> p1={eps,eps,0.}; double pt1d=0.;
//     std::vector<double> p2={eps,eps,h }; double pt2d=0.;
//     branch b1(p1,p2,ndiv, pt1d,pt2d);
//     branches.push_back(b1);
//     for (int iring=0; iring<2;iring++)
//     {
//         int ntheta=6;
//         double dtheta=2*3.14/ntheta;
//         double r=(iring+1)*p;
//         double xc=eps, yc=-eps;    
//         for(int k=0;k<ntheta;k++){
//             std::vector<double> p3={r*cos(k*dtheta)+xc,r*sin(k*dtheta)+yc,0};double pt3d=0.;
//             std::vector<double> p4={r*cos(k*dtheta)+xc,r*sin(k*dtheta)+yc,h};double pt4d=0.;
//             branch b2(p3,p4,ndiv, pt3d,pt4d);
//             branches.push_back(b2);
//             //midlle points
//             for(int ipt=0; ipt<iring; ipt++){
//                 double loc_dx=r*(cos((k+1)*dtheta) - cos(k*dtheta))/(iring+1);
//                 double loc_dy=r*(sin((k+1)*dtheta) - sin(k*dtheta))/(iring+1);
//                 std::vector<double> p3_b={r*cos(k*dtheta)+xc + loc_dx*(ipt+1),
//                                           r*sin(k*dtheta)+yc + loc_dy*(ipt+1),0};
//                 std::vector<double> p4_b={r*cos(k*dtheta)+xc + loc_dx*(ipt+1),
//                                           r*sin(k*dtheta)+yc + loc_dy*(ipt+1),h};
//                 branch b_b(p3_b,p4_b,ndiv, pt3d,pt4d);
//                 branches.push_back(b_b);
//             }
//         }
//     }
    
    
    double eps=1.e-2;double h=1.;int ndiv=10;
    //central branch
    std::vector<double> p1={0,0.5+eps,0.5+2*eps}; double pt1d=10.;
    std::vector<double> p2={h,0.5 + eps,0.5 + 2*eps }; double pt2d=50.;
    branch b1(p1,p2,ndiv, pt1d,pt2d);
    branches.push_back(b1);
    print(branches); 
    
    std::cout<<"Endl of Generator of branches"<<std::endl;
    
    return 1;}
