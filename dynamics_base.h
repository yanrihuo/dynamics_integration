#pragma once

#include "utility/utility.h"
//#include "parameters.h"

#include <ceres/ceres.h>
using namespace Eigen;
using namespace std;

class DynamicsBase
{
public:
	double dt;
	Eigen::Vector3d vel0, omiga0;
	Eigen::Vector3d vel1, omiga1;

	const Eigen::Vector3d linearized_vel, linearized_omiga;
	//Eigen::Vector3d linearized_ba, linearized_bg;

	Eigen::Matrix<double, 15, 15> jacobian, covariance;
	Eigen::Matrix<double, 15, 15> step_jacobian;
	Eigen::Matrix<double, 15, 18> step_V;
	Eigen::Matrix<double, 6, 6> noise;

	double sum_dt;
	Eigen::Vector3d delta_p;
	Eigen::Quaterniond delta_q;
	Eigen::Vector3d delta_v;

	std::vector<double> dt_buf;
	std::vector<Eigen::Vector3d> vel_buf;
	std::vector<Eigen::Vector3d> omiga_buf;
	DynamicsBase() = delete;
	DynamicsBase(const Eigen::Vector3d &_vel0, const Eigen::Vector3d &_omiga0)
		: vel0{ _vel0 }, omiga0{ _omiga0 }, linearized_vel{ _vel0 }, linearized_omiga{ _omiga0 },
		jacobian{ Eigen::Matrix<double, 15, 15>::Identity() }, covariance{ Eigen::Matrix<double, 15, 15>::Zero() },
		sum_dt{ 0.0 }, delta_p{ Eigen::Vector3d::Zero() }, delta_q{ Eigen::Quaterniond::Identity() }, delta_v{ Eigen::Vector3d::Zero() }

	{
		noise = Eigen::Matrix<double, 6, 6>::Zero();
		noise.block<3, 3>(0, 0) = (0.1 * 0.1) * Eigen::Matrix3d::Identity();
		noise.block<3, 3>(3, 3) = (0.1 * 0.1) * Eigen::Matrix3d::Identity();
	}

	void push_back(double dt, const Eigen::Vector3d &vel, const Eigen::Vector3d &omiga)
	{
		//cout<<dt<<endl;
		dt_buf.push_back(dt);
		//cout<<dt<<endl;
		vel_buf.push_back(vel);
		omiga_buf.push_back(omiga);

		propagate(dt, vel, omiga);
		cout<<delta_p<<endl<<endl<<endl;
	}

	/*void repropagate(const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg)
	{
		sum_dt = 0.0;
		acc_0 = linearized_acc;
		gyr_0 = linearized_gyr;
		delta_p.setZero();
		delta_q.setIdentity();
		delta_v.setZero();
		linearized_ba = _linearized_ba;
		linearized_bg = _linearized_bg;
		jacobian.setIdentity();
		covariance.setZero();
		for (int i = 0; i < static_cast<int>(dt_buf.size()); i++)
			propagate(dt_buf[i], acc_buf[i], gyr_buf[i]);
	}*/

	void dynamics_Integration(double _dt,
		const Eigen::Vector3d &_vel0, const Eigen::Vector3d &_omiga0,
		/*const Eigen::Vector3d &_vel1, const Eigen::Vector3d &_omiga1,*/
		const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
		Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
		bool update_jacobian)
	{
		//ROS_INFO("midpoint integration");
		Vector3d un_vel0 = delta_q * _vel0;
		//Vector3d un_omiga = 0.5 * (_omiga0 + _omiga1);
		result_delta_q = delta_q * Quaterniond(1, _omiga0(0) * _dt / 2, _omiga0(1) * _dt / 2, _omiga0(2) * _dt / 2);
		//Vector3d un_vel1 = result_delta_q * _vel1;
		//Vector3d un_vel = 0.5 * (un_vel0 + un_vel1);
		result_delta_p = delta_p + un_vel0 *_dt;

		//result_delta_v = delta_v + un_acc * _dt;
		//result_linearized_ba = linearized_ba;
		//result_linearized_bg = linearized_bg;

		if (update_jacobian)
		{
			Vector3d w_x = _omiga0;
			Vector3d v_0_x = _vel0;
			Matrix3d R_w_x, R_v_0_x;


			R_w_x << 0, -w_x(2), w_x(1),
				w_x(2), 0, -w_x(0),
				-w_x(1), w_x(0), 0;
			R_v_0_x << 0, -v_0_x(2), v_0_x(1),
				v_0_x(2), 0, -v_0_x(0),
				-v_0_x(1), v_0_x(0), 0;


			MatrixXd F = MatrixXd::Zero(15, 15);
			F.block<3, 3>(0, 0) = Matrix3d::Identity();
			F.block<3, 3>(0, 3) = -delta_q.toRotationMatrix()*R_v_0_x *_dt;
			F.block<3, 3>(3, 3) = Matrix3d::Identity() - R_w_x * _dt;

			//cout<<"A"<<endl<<A<<endl;

			MatrixXd V = MatrixXd::Zero(15, 6);
			V.block<3, 3>(0, 0) = delta_q.toRotationMatrix() * _dt;
			V.block<3, 3>(3, 3) = MatrixXd::Identity(3, 3) * _dt;

			//step_jacobian = F;
			//step_V = V;
			jacobian = F * jacobian;
			covariance = F * covariance * F.transpose() + V * noise * V.transpose();
		}

	}

	void propagate(double _dt, const Eigen::Vector3d &_vel0, const Eigen::Vector3d &_omiga0)
	{
		dt = _dt;
		vel0 = _vel0;
		omiga0 = _omiga0;
		Vector3d result_delta_p;
		Quaterniond result_delta_q;
		Vector3d result_delta_v;
		Vector3d result_linearized_ba;
		Vector3d result_linearized_bg;


		
		dynamics_Integration(_dt,vel0, omiga0, delta_p, delta_q, delta_v, 
			result_delta_p, result_delta_q, result_delta_v, 1);
		//checkJacobian(_dt, acc_0, gyr_0, acc_1, gyr_1, delta_p, delta_q, delta_v,
		//                    linearized_ba, linearized_bg);
		delta_p = result_delta_p;
		delta_q = result_delta_q;
		delta_v = result_delta_v;
		//linearized_ba = result_linearized_ba;
		//linearized_bg = result_linearized_bg;
		delta_q.normalize();
		sum_dt += dt;
		//acc_0 = acc_1;
		//gyr_0 = gyr_1;

	}

	Eigen::Matrix<double, 15, 1> evaluate(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi, const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
		const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj, const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj)
	{
		Eigen::Matrix<double, 15, 1> residuals;
		residuals.setZero();
		residuals.block<3, 1>(0, 0) = Qi.inverse() * (Pj - Pi) - delta_p;
		residuals.block<3, 1>(3, 0) = 2 * (delta_q.inverse() * (Qi.inverse() * Qj)).vec();
		return residuals;
	}
};