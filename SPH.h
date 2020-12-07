#pragma once
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <GL/glut.h>
#include "particle_node.h"
#include "vec3.hpp"
#include "RigidBody.h"
using namespace std;

namespace SPH {
	//盒子的大小：x=0,x=l,y=0,y=m.z=0,z=h,内部点为(l)*(m)*(h)个. | The size of Box: x=[0, l], y = [0, m], z = [0, h]
	static const int l = 30, m = 60, h = 40;
	static const double  Lx = 3.0, Ly = 6.0, Lz = 4.0;
	//定义pi的值 | Define the value of pi
	static const double pi = 3.14159265358979323846264338327950;
	static const P3 g_accelerate = P3(0.0, 0.0, -9.8); //重力加速度 | The gravity accelerate.

	class Sph {
	public:
		vector<particle_node> Nodes; //fluid nodes and bounadry nodes
		Wheel m_rigidbody;
		vector<P3> old_u;
		vector<int> idx_table[l + 1][m + 1][h + 1]; //一张索引表，存放了每个网格内部包含点的索引 \
			                                        //| A index table, which stored all the index of particles in corresponding grids.
		int m_grid_num[l + 1][m + 1][h + 1];
		int m_step; //当前已运行步数 | The step number
		int t_max;//最大的迭代计算步数 | The max steps this programm run.
		int m_ratio_h_dx; // = h/dx
		double m_h, m_dt, m_lambda, m_stiffness, m_KK, m_gamma;
		double dx, dy, dz;
		double m_Kpoly6;
		double a_Kpoly6;
		double vis_Kpoly6;
		double Vmax; //流场中的最大速度 | The max velocity of the fluid region.
		double radius; //绘制图像时的小球半径 | The radius of sphere while drawing the picture.
		double m_mass; //初始流体质量 | The initial mass of fluid particles.
		double rho_L; //初始的流体密度 | The initial density of fluid.
		double rho_R; //初始的刚体密度 | The initial density of rigid body.
		int total_num; //初始粒子总数为0 | The initial number of all the particles is 0.
		double p0; //初始压强 | The initial pressure of particles.
		double vis0; //初始粘度 | The initial viscosity of particles.
		std::string boundary_band; //including 6 char, denotes the boundary of Left/Right/Forward/Back/Down/Up\
								   // P: Periodical; W: Wall.

		bool if_dump; //是否要在运行过程中输出文件 | If we need to dump density files while we are running the code.
		int Nwri; // 每隔多少步输出一次文件 | Every "Nwri" steps passed, a density file will be dumped.

		Sph() :
			m_ratio_h_dx(2),
			dx(Lx / l),
			dy(Ly / m),
			dz(Lz / h),
			m_lambda(0.4),
			m_stiffness(1000.0),
			m_mass(1.0),
			m_gamma(1.0),
			radius(0.4),
			total_num(0),
			p0(0.0),
			vis0(0.002),
			Vmax(0.0),
			boundary_band("WWWWWW"),
			if_dump(false),
			Nwri(100),
			m_step(0),
			t_max(5000)
		{
			m_h = m_ratio_h_dx * dx;
			rho_L = m_mass / dx / dy / dz;
			rho_R = 300.0;
			m_dt = m_lambda * m_h / sqrt(m_stiffness);
			m_KK = m_stiffness * rho_L / m_gamma;
			m_Kpoly6 = 315 / ((64 * pi) * pow(m_h, 9));
			a_Kpoly6 = 45 / (pi * pow(m_h, 6));
			vis_Kpoly6 = 15.0 / (2.0 * pi * pow(m_h, 3));
			//m_rigidbody = Sphere(rho_R, Ly / 6.0, P3(Lx / 2.0, 2.0 * Ly / 3.0, Ly / 6.0 + 2 * dx));
			m_rigidbody = Wheel(rho_R, Lx / 6.0, Lx / 12.0, Lz / 5.0, P3(Lx / 2.0, Ly / 2.0, 0.0), 4);
			for (int i = 0; i <= l; i++) {
				for (int j = 0; j <= m; j++) {
					for (int r = 0; r <= h; r++) {
						m_grid_num[i][j][r] = 0;
					}
				}
			}

			initialization();
		}

		void initialization();

		void compute_press_accelerate(int k);

		void compute_vis_accelerate(int k);

		void Compute_Rho(int k);

		void Set_BoundaryNodes();

		void add_node(int i, int j, int r, char status);

		void dump_file(string file_name);
		bool read_file(string file_name);

		void drawing_liquid(void);

		double comp_error();

		void record_velocity();

		void compute(void);

		void step();
		void OnlyReadFileStep();

		int getstep() { return m_step; }
		int getMaxstep() { return t_max; }

		double Muller03Kernel_Basic(P3 origin_point, P3 detect_point) {
			double r_square_norm = (detect_point - origin_point).squarenorm();
			if (m_h * m_h > r_square_norm) { return m_Kpoly6 * pow(m_h * m_h - r_square_norm, 3); }
			else { return 0.0; }
		}

		P3 Gradient_Muller03Kernel_Basic(P3 origin_point, P3 detect_point) {
			double r_square_norm = (detect_point - origin_point).squarenorm();
			if (m_h * m_h > r_square_norm) {
				return (detect_point - origin_point) * (6.0 * m_Kpoly6 * pow(m_h * m_h - r_square_norm, 2));
			}
			else { return P3(0.0, 0.0, 0.0); }
		}

		double Muller03Kernel_Pressure(P3 origin_point, P3 detect_point) {
			double r_norm = (origin_point - detect_point).norm();
			if (m_h > r_norm) { return a_Kpoly6 / 3.0 * pow((m_h - r_norm), 3); }
			else { return 0.0; }
		}

		P3 Gradient_Muller03Kernel_Pressure(P3 origin_point, P3 detect_point) {
			double r_norm = (origin_point - detect_point).norm();
			if (r_norm > m_h) { return P3(0.0, 0.0, 0.0); }
			return (detect_point - origin_point).normalized() * (a_Kpoly6 * pow(m_h - r_norm, 2));
		}

		P3 Gradient_Muller03Kernel_Vis(P3 origin_point, P3 detect_point) {
			double r_norm = (origin_point - detect_point).norm();
			if (r_norm > m_h) { return P3(0.0, 0.0, 0.0); }
			return (detect_point - origin_point) * (vis_Kpoly6 * \
				(-1.5 * r_norm / pow(m_h, 3) + 2.0 / m_h / m_h - 0.5 * m_h / pow(r_norm, 3)));
		}
	};
}