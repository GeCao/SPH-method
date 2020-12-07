#pragma once
#ifndef __VEC3_H__
#define __VEC3_H__ 
using namespace std;
class P3{
private:
	double x, y, z;
public:

	//初始化函数
	P3(double x_0 = 0, double y_0 = 0, double z_0 = 0) : x(x_0), y(y_0), z(z_0){}

	//算子重载：
	//===================================
	//====        "+-" for P3        ====
	//===================================
	P3 operator+(const P3& a)const {
		return P3(a.x + x, a.y + y, a.z + z);
	}
	P3 operator-(const P3& a)const {
		return P3(x - a.x, y - a.y, z - a.z);
	}
	P3& operator+=(const P3&a) { return *this = *this + a; }
	P3& operator-=(const P3&a) { return *this = *this - a; }
	P3& operator*=(double p) { return *this = *this * p; }
	P3& operator/=(double p) { return *this = *this / p; }
	//===================================
	//====    "+-*/" for constant    ====
	//===================================
	P3 operator+(const double a)const {
		return P3(a + x, a + y, a + z);
	}
	P3 operator-(const double a)const {
		return P3(x - a, y - a, z - a);
	}
	P3 operator*(const double a)const {
		return P3(a * x, a * y, a * z);
	}
	P3 operator/(const double a)const {
		if (a == 0){ std::cout << "There's a problem that a number divides 0!" << std::endl; exit(1); }
		return P3(x / a, y / a, z / a);
	}
	//===================================
	//====         By reverse        ====
	//===================================
	P3 operator-()const {
		return P3(-x, -y, -z);
	}


	friend ostream& operator<<(ostream &os, const P3& a)
	{
		os << "Point/Vector::(" << a.x << ", " << a.y << ", " << a.z << " )";
		return os;
	}


	bool operator==(const P3&a) const { return x == a.x && y == a.y && z == a.z; }
	bool operator!=(const P3&a) const { return x != a.x || y != a.y || z != a.z; }

	//函数部分
	double maxCoeff() const {//返回三维向量中的最大值
		return x>y&&x>z ? x : y>z ? y : z; 
	}
	double maxAbsCoeff() const {//返回三维向量中的最大绝对值
		return abs(x) > abs(y) && abs(x) > abs(z) ? abs(x) : abs(y) > abs(z) ? abs(y) : abs(z);
	}
	double dot(const P3& a, const P3& b)const {//向量的点乘
		return (a.x*b.x + a.y*b.y + a.z*b.z);
	}
	double dot(const P3& a)const {//向量的点乘
		return (a.x*x + a.y*y + a.z*z);
	}
	P3 mult(const P3&a) const { return P3(x*a.x, y*a.y, z*a.z); }
	double norm()const {//向量的模长
		return sqrt(x * x + y * y + z * z);
	}
	double squarenorm()const {//向量模长的平方
		return x * x + y * y + z * z;
	}
	P3 normalized() {//向量的正则化
		double n_norm = this->norm();
		x = x / n_norm; y = y / n_norm; z = z / n_norm;
		return *this;
	}
	P3 cross_product(const P3& a)const {//向量的叉乘
		return P3(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	}
	P3 clip(double r0 = 0, double r1 = 1) const { //进行范围内的切片
		return P3(x>r1 ? r1 : x<r0 ? r0 : x,
			y>r1 ? r1 : y<r0 ? r0 : y,
			z>r1 ? r1 : z<r0 ? r0 : z); 
	}
	double getx() const { return x; }
	double gety() const { return y; }
	double getz() const { return z; }

	double setx(double x_) { x = x_; }
	double sety(double y_) { y = y_; }
	double setz(double z_) { z = z_; }
};

#endif // __VEC3_H__