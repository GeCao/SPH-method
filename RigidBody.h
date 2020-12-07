#pragma once
#include "vec3.hpp"
#include "mat3.hpp"
#include "particle_node.h"
class RigidBody {
public:
	P3 b, v; //质心位置 & 线速度
	//vector<double> Quaternion; //四元数，代替旋转矩阵
	mat3 Rotation;
	P3 angular_velocity; //角速度

	P3 linear_momentum, linear_force; //线动量 & 施加力
	P3 angular_momentum, torque; //角动量 & Torque

	mat3 Ib_inv;

	RigidBody() :
		b(P3()),
		v(P3()),
		Rotation(mat3()),
		angular_velocity(P3()),
		linear_momentum(P3()),
		linear_force(P3()),
		angular_momentum(P3()),
		torque(P3()),
		Ib_inv(mat3())
	{
		
	}

	virtual void update(double m_dt) = 0;
};

class Wheel :public RigidBody {
public:
	std::vector<int> shape; //Store the index rigid body nodes
	P3 m_wheel_center; //The center of wheel.
	double m_wheel_radius_outsize; //The radius of wheel.(Out circle)
	double m_wheel_radius_insize;//The radius of wheel.(In circle)
	double m_wheel_height;//The height of wheel
	int m_leafnum; //The number of leaves in the wheel.
	double swirl_velocity;
	double m_Mass;

	Wheel() {
		m_wheel_center = P3();
		m_wheel_radius_outsize = 0.0;
		m_wheel_radius_insize = 0.0;
		m_wheel_height = 0.0;
		m_leafnum = 0;
		m_Mass = 1.0;
		shape.clear();
		Ib_inv = mat3();

		b = m_wheel_center;
	}

	Wheel(double rho_R, double radius_outsize, double radius_insize, double height, P3 center, int leaf_num) {
		m_wheel_center = center;
		m_wheel_radius_outsize = radius_outsize;
		m_wheel_radius_insize = radius_insize;
		m_wheel_height = height;
		m_leafnum = leaf_num;
		m_Mass = rho_R * height * 4.0 * 3.1415926 * \
			(m_wheel_radius_insize * m_wheel_radius_insize + m_wheel_radius_outsize * m_wheel_radius_outsize) / 2.0;
		shape.clear();
		Ib_inv = mat3();

		b = m_wheel_center;
	}

	virtual void update(double m_dt) { //Semi-implicit
		swirl_velocity = 3.1415926 / 18.0 / m_dt;
	}
};

class Sphere :public RigidBody {
public:
	std::vector<int> shape; //Store the index rigid body nodes
	P3 m_sphere_center; //The center of sphere.
	double m_sphere_radius; //The radius of sphere.
	double m_Mass;

	Sphere(){
		m_sphere_radius = 0.0;
		m_Mass = 1.0;
		m_sphere_center = P3();
		shape.clear();
		Ib_inv = mat3() * (0.4 * m_Mass * m_sphere_radius * m_sphere_radius);

		b = m_sphere_center;
	}

	Sphere(double rho_R, double radius, P3 center) {
		m_sphere_radius = radius;
		m_Mass = rho_R * (4.0 * 3.1415926 / 3.0 * pow(radius, 3));
		m_sphere_center = center;
		shape.clear();
		Ib_inv = mat3() * (0.4 * m_Mass * m_sphere_radius * m_sphere_radius);

		b = m_sphere_center;
	}

	virtual void update(double m_dt) { //Semi-implicit
		mat3 I_inv = Rotation.multi(Ib_inv).multi(Rotation);

		//已经更新的变量：linear force & torque: 2 variables done!
		linear_momentum += linear_force * m_dt; // 3 variables done!
		angular_momentum += torque * m_dt; // 4 variables done!
		cout << "forcces: " << linear_force << std::endl;

		torque -= angular_velocity.cross_product(I_inv.multi(angular_velocity)); //加入科里奥利力
		angular_velocity += I_inv.multi(torque) * m_dt; // 5 variables done!

		v = linear_momentum / m_Mass; // 6 variables done!
		b += v * m_dt; // 7 variables done!
		m_sphere_center = b;

		mat3 Omega(std::vector<P3>{
			P3(0.0, -angular_velocity.getz(), angular_velocity.gety()),
		    P3(angular_velocity.getz(), 0.0, -angular_velocity.getx()),
		    P3(-angular_velocity.gety(), angular_velocity.getx(), 0.0)});
		Rotation += Omega.multi(Rotation) * m_dt;
	}
};