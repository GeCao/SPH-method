#pragma once
#include "vec3.hpp"
namespace SPH{
	class particle_node {
	public:
		P3 u; //粒子速度 | velocity of particle
		P3 place; //粒子位置 | place of particle
		double rho; //粒子密度 | density of particle
		double p; //压强 | Pressure
		P3 force_vis;
		P3 force_press;
		P3 acc_ext;
		char status; //粒子的状态属性（固体/液体/气体）| The status of particle (solid/liquid/gas)
		vector<int> neighbor_index;

		particle_node(double rho0, double p0, P3 acc_g):
			u(P3(0.0, 0.0, 0.0)),
			place(P3(0.0, 0.0, 0.0)),
			force_vis(P3(0.0, 0.0, 0.0)),
			force_press(P3(0.0, 0.0, 0.0)),
			acc_ext(acc_g),
			rho(rho0),
			p(p0),
			status('L')
		{
		}

		//改变粒子的位置信息
		void setplace(double new_x, double new_y, double new_z) {
			place = P3(new_x, new_y, new_z);
		}

		void setplace(P3 new_place) {
			place = new_place;
		}

		void set_status(char ch) {
			status = ch;
		}

		bool isWall() {
			return status == 'W';
		}
		bool isRigidBody() {
			return status == 'R';
		}
		bool isFluid() {
			return status == 'L';
		}
	};
}
