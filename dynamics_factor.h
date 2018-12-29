#pragma once
#include <iostream>
#include <eigen3/Eigen/Dense>

#include "utility/utility.h"
//#include "parameters.h"
#include "dynamics_base.h"

#include <ceres/ceres.h>

class dynamicsFactor : public ceres::SizedCostFunction<6, 7, 7>
{
public:
	dynamicsFactor() = delete;
	dynamicsFactor(DynamicsBase* _dynamics_integration) :dynamics_integration(_dynamics_integration)
	{
	}
	//IMU对应的残差，对应ceres的结构，需要自己计算jacobian
	virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
	{

		Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
		Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

		//Eigen::Vector3d Vi(parameters[1][0], parameters[1][1], parameters[1][2]);
		//Eigen::Vector3d Bai(parameters[1][3], parameters[1][4], parameters[1][5]);
		//Eigen::Vector3d Bgi(parameters[1][6], parameters[1][7], parameters[1][8]);

		Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
		Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

		//Eigen::Vector3d Vj(parameters[3][0], parameters[3][1], parameters[3][2]);
		//Eigen::Vector3d Baj(parameters[3][3], parameters[3][4], parameters[3][5]);
		//Eigen::Vector3d Bgj(parameters[3][6], parameters[3][7], parameters[3][8]);

		//Eigen::Matrix<double, 15, 15> Fd;
		//Eigen::Matrix<double, 15, 12> Gd;

		//Eigen::Vector3d pPj = Pi + Vi * sum_t - 0.5 * g * sum_t * sum_t + corrected_delta_p;
		//Eigen::Quaterniond pQj = Qi * delta_q;
		//Eigen::Vector3d pVj = Vi - g * sum_t + corrected_delta_v;
		//Eigen::Vector3d pBaj = Bai;
		//Eigen::Vector3d pBgj = Bgi;

		//Vi + Qi * delta_v - g * sum_dt = Vj;
		//Qi * delta_q = Qj;

		//delta_p = Qi.inverse() * (0.5 * g * sum_dt * sum_dt + Pj - Pi);
		//delta_v = Qi.inverse() * (g * sum_dt + Vj - Vi);
		//delta_q = Qi.inverse() * Qj;


		Eigen::Map<Eigen::Matrix<double, 6, 1>> residual(residuals);
		residual = dynamics_integration->evaluate(Pi, Qi, Pj, Qj);

		Eigen::Matrix<double, 6, 6> sqrt_info = Eigen::LLT<Eigen::Matrix<double, 6, 6>>(dynamics_integration->covariance.inverse()).matrixL().transpose();
		//sqrt_info.setIdentity();
		residual = sqrt_info * residual;

		if (jacobians)
		{
			//cout<<"!!!"<<endl;
			double sum_dt = dynamics_integration->sum_dt;
			/*Eigen::Matrix3d dp_dba = dynamics_integration->jacobian.template block<3, 3>(O_P, O_BA);
			Eigen::Matrix3d dp_dbg = dynamics_integration->jacobian.template block<3, 3>(O_P, O_BG);

			Eigen::Matrix3d dq_dbg = dynamics_integration->jacobian.template block<3, 3>(O_R, O_BG);

			Eigen::Matrix3d dv_dba = dynamics_integration->jacobian.template block<3, 3>(O_V, O_BA);
			Eigen::Matrix3d dv_dbg = dynamics_integration->jacobian.template block<3, 3>(O_V, O_BG);*/

			if (dynamics_integration->jacobian.maxCoeff() > 1e8 || dynamics_integration->jacobian.minCoeff() < -1e8)
			{
				std::cout<<"numerical unstable in preintegration"<<std::endl;
				//std::cout << pre_integration->jacobian << std::endl;
				///                ROS_BREAK();
			}

			if (jacobians[0])
			{
				Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
				jacobian_pose_i.setIdentity();

				jacobian_pose_i.block<3, 3>(0, 0) = -Qi.inverse().toRotationMatrix();
				jacobian_pose_i.block<3, 3>(0, 3) = Utility::skewSymmetric(Qi.inverse() * (Pj - Pi ));
				jacobian_pose_i.block<3, 3>(3, 3) = -(Utility::Qleft(Qj.inverse() * Qi) * Utility::Qright(dynamics_integration->delta_q)).bottomRightCorner<3, 3>();
				//cout<<"qleft:"<<Utility::Qleft(Qj.inverse() * Qi)<<endl;
				jacobian_pose_i = sqrt_info * jacobian_pose_i;

				if (jacobian_pose_i.maxCoeff() > 1e8 || jacobian_pose_i.minCoeff() < -1e8)
				{
					std::cout<<"numerical unstable in preintegration"<<std::endl;
					//std::cout << sqrt_info << std::endl;
					return 0;
				}
				//cout<<"i:"<<jacobian_pose_i<<endl;
			}
			/*if (jacobians[1])
			{
				Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_i(jacobians[1]);
				jacobian_speedbias_i.setZero();
			}*/
			if (jacobians[1])
			{
				Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);
				jacobian_pose_j.setIdentity();

				jacobian_pose_j.block<3, 3>(0, 0) = Qi.inverse().toRotationMatrix();
				jacobian_pose_j.block<3, 3>(3, 3) = Utility::Qleft(dynamics_integration->delta_q.inverse() * Qi.inverse() * Qj).bottomRightCorner<3, 3>();
				jacobian_pose_j = sqrt_info * jacobian_pose_j;
				//cout<<"j:"<<jacobian_pose_j<<endl;
			}
			/*if (jacobians[3])
			{
				Eigen::Map<Eigen::Matrix<double, 15, 9, Eigen::RowMajor>> jacobian_speedbias_j(jacobians[3]);
				jacobian_speedbias_j.setZero();
			}*/
		}

		return true;
	}

	//bool Evaluate_Direct(double const *const *parameters, Eigen::Matrix<double, 15, 1> &residuals, Eigen::Matrix<double, 15, 30> &jacobians);

	//void checkCorrection();
	//void checkTransition();
	//void checkJacobian(double **parameters);
	DynamicsBase* dynamics_integration;

};