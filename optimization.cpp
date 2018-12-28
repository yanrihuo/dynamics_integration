#include "utility/utility.h"
#include "utility/tic_toc.h"
#include <ceres/ceres.h>
#include "pose_local_parameterization.h"
#include "dynamics_base.h"
#include "dynamics_factor.h"
#include <unordered_map>
#include <queue>
#include <opencv2/core/eigen.hpp>
#include <math.h>
//#include <sstream>
//#include <iomanip>
using namespace std;
const double k1=30082*2;//qianlun
const double k2=31888*2;//houlun
const double m=1096;//zhiliang
const double l=2.3;
const double a=1.0377;//qianzhou
const double b=1.2623;//houzhou
const double K=m*(a/k2-b/k1)/(l*l);
int main(int argc, char** argv)
{
	Eigen::Vector3d v;
	vector <Eigen::Vector3d> vehicle_data;
	FILE *fp1;
	fp1 = fopen("CanData(copy).csv", "r");//
	fscanf(fp1, "%lf,%lf,%lf", &v[0], &v[1], &v[2]);
	while (1)
	{
		fscanf(fp1, "%lf,%lf,%lf", &v(0), &v(1), &v(2));
		//cout<<setprecision(11)<<v(0)<<endl;
		//v(1)=v(1)/3.6;
		//v(2)=v(2)*3.14159/180;
		vehicle_data.push_back(v);
		if (feof(fp1)) 
			break;
	}
	fclose(fp1);
	int m=vehicle_data.size()/5;
	double para_Pose[(m+1)][7]={0};
	double para_SpeedBias[(m+1)][9]={0};
	DynamicsBase *dynamics_integrations[(m+1)];
	ceres::Problem problem;
	ceres::LossFunction *loss_function;
	loss_function=new ceres::HuberLoss(1.0);
	Vector3d vel(0,0,0);
	Vector3d omiga(0,0,0);
	double beta;
	double t0=vehicle_data[0](0);
	double dt;
	//cout<<vehicle_data.size()<<endl;
	//cout<<sizeof(dynamics_integrations)<<endl;
	for(int i=0; i<vehicle_data.size(); i++)
	{
		
		beta=(1-m*vehicle_data[i](1)*vehicle_data[i](1)*a/(2*l*b*k2*3.6*3.6))*b*v(2)*3.14159/180/17.4
		/l/(1+K*vehicle_data[i](1)*vehicle_data[i](1)/3.6/3.6);
		omiga(0)=vehicle_data[i](1)/3.6*v(2)*3.14159/180/17.4/l/(1+K*vehicle_data[i](1)*vehicle_data[i](1)/3.6/3.6);
		vel(1)=vehicle_data[i](1)/3.6*sin(beta);
		vel(2)=vehicle_data[i](1)/3.6*cos(beta);
		dt=(vehicle_data[i](0)-t0)/1000000000;
		t0=vehicle_data[i](0);
		int j=i/5;
		//cout<<j<<endl;
		if(i%5==0) 
		{
			dynamics_integrations[j]=new DynamicsBase{vel,omiga};
			//cout<<j<<endl;
		}
		//cout<<i<<endl;
		cout<<"yujifen:"<<j<<endl;
		dynamics_integrations[j]->push_back(dt,vel,omiga);
	}


	for(int i=0; i<vehicle_data.size()/5; i++)
	{
		
		ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
		problem.AddParameterBlock(para_Pose[i], 7, local_parameterization);
		problem.AddParameterBlock(para_SpeedBias[i], 9);
	}
	for(int i=0; i<vehicle_data.size()/5; i++)
	{
		dynamicsFactor* dynamics_factor = new dynamicsFactor(dynamics_integrations[i]);
		problem.AddResidualBlock(dynamics_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[i+1], para_SpeedBias[i+1]);
	}


	ceres::Solver::Options options;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 20;
    //options.use_explicit_schur_complement = true;
    //options.minimizer_progress_to_stdout = true;
    //options.use_nonmonotonic_steps = true;
    /*
    if (marginalization_flag == MARGIN_OLD)
        options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0;
    else
        options.max_solver_time_in_seconds = SOLVER_TIME;
    */
    //TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    
    /*for (int i = 0; i < vehicle_data.size()/5 + 1; ++i)
    {
    	cout<<para_Pose[i]<<endl;
    }*/
    cout<<summary.BriefReport() <<endl;
}