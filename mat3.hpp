#pragma once
#ifndef __MAT3_H__
#define __MAT3_H__ 
#include "vec3.hpp"
class mat3
{
public:
	mat3() {
		m_matrix = std::vector<P3>{ P3(1.0, 0.0, 0.0), P3(0.0, 1.0, 0.0), P3(0.0, 0.0, 1.0) };
	}
	mat3(std::vector<P3> other_mat) {
		m_matrix = other_mat;
	}

	mat3 Identity() const {
		return mat3(std::vector<P3>{P3(1.0, 0.0, 0.0), P3(0.0, 1.0, 0.0), P3(0.0, 0.0, 1.0)});
	}
	mat3 Zero() const {
		return mat3(std::vector<P3>(3, P3(0.0, 0.0, 0.0)));
	}

	P3 multi(P3 other_vector) const {
		return P3(m_matrix[0].dot(other_vector),
			m_matrix[1].dot(other_vector),
			m_matrix[2].dot(other_vector));
	}
	mat3 T() const {
		std::vector<P3> new_matrix(3, P3());
		for (int i = 0; i < 3; i++) {
			new_matrix[i] = this->col(i);
		}
		return mat3(new_matrix);
	}
	mat3 multi(mat3 other_mat) const {
		std::vector<P3> new_matrix(3, P3());
		for (int i = 0; i < 3; i++) {
			new_matrix[i] = P3(m_matrix[i].dot(other_mat.col(0)),
				m_matrix[i].dot(other_mat.col(1)),
				m_matrix[i].dot(other_mat.col(2)));
		}
		return mat3(new_matrix);
	}

	mat3 operator*(const double a) const {
		std::vector<P3> new_matrix(3, P3());
		for (int i = 0; i < 3; i++) {
			new_matrix[i] = m_matrix[i] * a;
		}
		return mat3(new_matrix);
	}
	mat3 operator+(const mat3 a)const {
		std::vector<P3> new_matrix(3, P3());
		for (int i = 0; i < 3; i++) {
			new_matrix[i] = m_matrix[i] + a.row(i);
		}
		return mat3(new_matrix);
	}
	mat3& operator+=(const mat3& a) {
		return *this = *this + a; 
	}

	P3 row(int i) const {
		if (i >= 3 || i < 0) {
			std::cout << "Wrong row index to get the row from matrix!" << std::endl;
			system("pause");
			exit(-1);
		}
		return m_matrix[i];
	}
	P3 col(int j) const {
		if (j >= 3 || j < 0) {
			std::cout << "Wrong row index to get the row from matrix!" << std::endl;
			system("pause");
			exit(-1);
		}
		switch (j)
		{
		case 0:
			return P3(m_matrix[0].getx(), m_matrix[1].getx(), m_matrix[2].getx());
		case 1:
			return P3(m_matrix[0].gety(), m_matrix[1].gety(), m_matrix[2].gety());
		case 2:
			return P3(m_matrix[0].getz(), m_matrix[1].getz(), m_matrix[2].getz());
		default:
			std::cout << "Wrong row index to get the row from matrix!" << std::endl;
			system("pause");
			exit(-1);
		}
	}

private:
	/*
	* [0.x, 0.y, 0.z]
	* [1.x, 1.y, 1.z]
	* [2.x, 2.y, 2.z]
	*/
	std::vector<P3> m_matrix;
};
#endif