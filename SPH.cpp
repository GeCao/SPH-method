#include "SPH.h"

namespace SPH {
	void Sph::add_node(int i, int j, int r, char status) {
		if (m_grid_num[i][j][r] == 1) {
			if (status != 'L') {
				if (!Nodes[idx_table[i][j][r][0]].isRigidBody() && status == 'R') {
					m_rigidbody.shape.push_back(idx_table[i][j][r][0]);
				}
				Nodes[idx_table[i][j][r][0]].set_status(status);
			}
			return;
		}
		particle_node temp_node(rho_L, p0, g_accelerate);

		idx_table[i][j][r].push_back(total_num++); //Push index into index table
		temp_node.setplace(i * dx, j * dy, r * dz);
		temp_node.set_status(status);
		Nodes.push_back(temp_node);//Push new node into node list.
		if(temp_node.isRigidBody()){ m_rigidbody.shape.push_back(Nodes.size() - 1); }
		m_grid_num[i][j][r] = 1;
	}

	void Sph::Set_BoundaryNodes() {
		//Special Zone: Wall
		if (boundary_band[0] == 'W') { //Left
			for (int j = 0; j <= m; j++)
				for (int r = 0; r <= h; r++) {
					add_node(0, j, r, 'W'); add_node(1, j, r, 'W');
				}
		}
		if (boundary_band[1] == 'W') { //Right
			for (int j = 0; j <= m; j++)
				for (int r = 0; r <= h; r++) {
					add_node(l, j, r, 'W'); add_node(l - 1, j, r, 'W');
				}
		}
		if (boundary_band[2] == 'W') { //Forward
			for (int i = 0; i <= l; i++)
				for (int r = 0; r <= h; r++) {
					add_node(i, 0, r, 'W'); add_node(i, 1, r, 'W');
				}
		}
		if (boundary_band[3] == 'W') { //Back
			for (int i = 0; i <= l; i++)
				for (int r = 0; r <= h; r++) {
					add_node(i, m, r, 'W'); add_node(i, m - 1, r, 'W');
				}
		}
		if (boundary_band[4] == 'W') { //Down
			for (int i = 0; i <= l; i++)
				for (int j = 0; j <= m; j++) {
					add_node(i, j, 0, 'W'); add_node(i, j, 1, 'W');
				}
		}
		if (boundary_band[5] == 'W') { //Up
			for (int i = 0; i <= l; i++)
				for (int j = 0; j <= m; j++) {
					add_node(i, j, h, 'W'); add_node(i, j, h - 1, 'W');
				}
		}
	}

	void Sph::initialization() {
		particle_node temp_node(rho_L, p0, g_accelerate);
		int i = 0, j = 0, r = 0;
		std::cout << "[Init]: Intitialization begin" << std::endl;

		//How to set boundary?
		// if m = 30 (in y direction), then the Boundary is set in idx = 0, 30, and the(1, 29) is fluid.
		Set_BoundaryNodes();
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				for (r = 0; r <= h; r++) {
					/*
					if (pow(i * dx - m_rigidbody.m_sphere_center.getx(), 2) +
						pow(j * dy - m_rigidbody.m_sphere_center.gety(), 2) +
						pow(r * dz - m_rigidbody.m_sphere_center.getz(), 2) < pow(m_rigidbody.m_sphere_radius, 2)) {
						add_node(i, j, r, 'R');
					}
					*/
					if (r * dz > 2 * dz + m_rigidbody.m_wheel_height) { continue; } //The 2*dy is because the boundary layer.
					if (r * dz < 2 * dz) { continue; }
					double local_x = i * dx - m_rigidbody.m_wheel_center.getx();
					double local_y = j * dy - m_rigidbody.m_wheel_center.gety();
					double r_norm2 = local_x * local_x + local_y * local_y;
					if (r_norm2 > m_rigidbody.m_wheel_radius_outsize * m_rigidbody.m_wheel_radius_outsize) { continue; }
					else if (r_norm2 <= m_rigidbody.m_wheel_radius_insize * m_rigidbody.m_wheel_radius_insize) {
						add_node(i, j, r, 'R');
					}
					else {
						double dtheta = 3.1415926 / m_rigidbody.m_leafnum;
						double r_norm = sqrt(r_norm2);
						double costheta = local_x / r_norm, sintheta = local_y / r_norm;
						double theta = acos(costheta);
						if (sintheta < 0.0) { theta = 2 * 3.1415926 - theta; }
						int dnum = (int)((theta + 2 * 3.1415926) / dtheta);
						if (dnum % 2 == 0) {
							add_node(i, j, r, 'R');
						}
					}
				}
			}
		}
		for (i = 0; i <= l; i++) {
			for (j = 0; j < m / 3; j++) {
				for (r = 0; r < 3 * h / 4; r++) {
					add_node(i, j, r, 'L');
				}
			}
		}

		std::cout << "The total number of particles: " << total_num << std::endl;
		old_u = std::vector<P3>(total_num, P3(0.0, 0.0, 0.0));
	}

	void Sph::record_velocity() {
		int k = 0;
#pragma omp parallel for private(k)
		for (k = 0; k < total_num; k++) {
			old_u[k] = Nodes[k].u;
		}
	}

	double Sph::comp_error()
	{
		int k = 0;
		double final_error = 0.0;

#pragma omp parallel for private(k)
		for (k = 0; k < total_num; k++) {
			final_error += (old_u[k] - Nodes[k].u).squarenorm();
		}
		return sqrt(final_error);
	}

	void Sph::compute_vis_accelerate(int k) {
		double Vb = 0.0;
		for (int i = 0; i < Nodes[k].neighbor_index.size(); i++) {
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if (!Nodes[neighbor_ID].isFluid()) { //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(Nodes[k].place, Nodes[neighbor_ID].place);
			}
		}
		Vb = 1.0 / Vb;

		for (int i = 0; i < Nodes[k].neighbor_index.size(); i++) {
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if (neighbor_ID == k) { continue; }

			double r_sqnorm = (Nodes[k].place - Nodes[neighbor_ID].place).squarenorm();
			if (r_sqnorm > m_h * m_h) { continue; }

			P3 Gradient_Wij = Gradient_Muller03Kernel_Pressure(Nodes[k].place, Nodes[neighbor_ID].place);
			double u_dot_r = (Nodes[k].u - Nodes[neighbor_ID].u).dot(Nodes[k].place - Nodes[neighbor_ID].place);
			if (Nodes[neighbor_ID].isFluid()) {
				//fluid - fluid
				Nodes[k].force_vis += Gradient_Wij * ((vis0) * 2.0 * (3 + 2.0) * (m_mass*m_mass / Nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
			}
			else {
				//fluid - obstacle
				double mass_rigid = m_mass;
				P3 temp_force = Gradient_Wij * ((vis0) * 2.0 * (3 + 2.0) * (mass_rigid * m_mass / Nodes[neighbor_ID].rho) * u_dot_r / r_sqnorm);
				Nodes[k].force_vis += temp_force;
				Nodes[neighbor_ID].force_vis += -temp_force;
			}
		}
	}

	void Sph::compute_press_accelerate(int k) {
		if (Nodes[k].p < 0.0) { return; }

		double Vb = 0.0;
		for (int i = 0; i < Nodes[k].neighbor_index.size(); i++) {
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if (!Nodes[neighbor_ID].isFluid()) { //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(Nodes[k].place, Nodes[neighbor_ID].place);
			}
		}
		Vb = 1.0 / Vb;

		for (int i = 0; i < Nodes[k].neighbor_index.size(); i++) {
			int neighbor_ID = Nodes[k].neighbor_index[i];

			if (neighbor_ID == k) { continue; }
			if ((Nodes[k].place - Nodes[neighbor_ID].place).squarenorm() > m_h * m_h) { continue; }

			P3 Gradient_Wij = Gradient_Muller03Kernel_Pressure(Nodes[k].place, Nodes[neighbor_ID].place);
			if (Nodes[neighbor_ID].isWall()) {
				//fluid - Obstacle
				double mass_rigid = m_mass;
				P3 temp_force = Gradient_Wij * (-2.0 * m_mass * mass_rigid * Nodes[k].p / Nodes[k].rho / Nodes[k].rho);
				Nodes[k].force_press += temp_force;
				Nodes[neighbor_ID].force_press += -temp_force;
			}
			else {
				//fluid - fluid
				Nodes[k].force_press += Gradient_Wij * (-m_mass * m_mass * \
					(Nodes[k].p / Nodes[k].rho / Nodes[k].rho + \
						Nodes[neighbor_ID].p / Nodes[neighbor_ID].rho / Nodes[neighbor_ID].rho));
			}
		}
	}

	void Sph::Compute_Rho(int k) {
		if (Nodes[k].isWall()) { return; }
		Nodes[k].rho = 0.0;
		double Vb = 0.0;
		for (int i = 0; i < Nodes[k].neighbor_index.size(); i++) {
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if (!Nodes[neighbor_ID].isFluid()) { //对运动刚体和固定边界的密度插值修正: Akinci et al. (2012).
				Vb += Muller03Kernel_Basic(Nodes[k].place, Nodes[neighbor_ID].place);
			}
		}
		Vb = 1.0 / Vb;
		for (int i = 0; i < Nodes[k].neighbor_index.size(); i++) {
			int neighbor_ID = Nodes[k].neighbor_index[i];
			if ((Nodes[k].place - Nodes[neighbor_ID].place).squarenorm() > m_h * m_h) { continue; }

			P3 Grad_W = Gradient_Muller03Kernel_Basic(Nodes[k].place, Nodes[neighbor_ID].place);
			if (!Nodes[neighbor_ID].isFluid()) {
				double mass_rigid = m_mass;
				Nodes[k].rho += mass_rigid * Muller03Kernel_Basic(Nodes[k].place, Nodes[neighbor_ID].place);
				Nodes[k].rho += mass_rigid * m_dt * (Nodes[k].u - Nodes[neighbor_ID].u).dot(Grad_W);
			}
			else {
				Nodes[k].rho += m_mass * Muller03Kernel_Basic(Nodes[k].place, Nodes[neighbor_ID].place);
				Nodes[k].rho += m_mass * m_dt * (Nodes[k].u - Nodes[neighbor_ID].u).dot(Grad_W);
			}
		}
	}

	void Sph::compute(void)
	{
		int i = 0, j = 0, r = 0, k = 0;

		//step0:更新表：
#pragma omp parallel for private(i,j,r)
		for (i = 0; i <= l; i++) {
			for (j = 0; j <= m; j++) {
				for (r = 0; r <= h; r++) {
					idx_table[i][j][r].clear();
					m_grid_num[i][j][r] = 0;
				}
			}
		}
		for (k = 0; k < total_num; k++)
		{
			i = (int)(Nodes[k].place.getx() / dx);
			j = (int)(Nodes[k].place.gety() / dy);
			r = (int)(Nodes[k].place.getz() / dz);
			if (i < 0 || i > l) { std::cout << "A particle flied out from x direction:" << Nodes[k].status << std::endl; continue; }
			if (j < 0 || j > m) { std::cout << "A particle flied out from y direction:" << Nodes[k].status << std::endl; continue; }
			if (r < 0 || r > h) { std::cout << "A particle flied out from z direction:" << Nodes[k].status << std::endl; continue; }
			idx_table[i][j][r].push_back(k);
			m_grid_num[i][j][r] += 1;
		}
		for (k = 0; k < total_num; k++) {
			Nodes[k].force_vis = P3(0.0, 0.0, 0.0);
			Nodes[k].neighbor_index.clear();
			int comp_i = (int)(Nodes[k].place.getx() / dx);
			int comp_j = (int)(Nodes[k].place.gety() / dy);
			int comp_r = (int)(Nodes[k].place.getz() / dz);
			for (i = ((comp_i >= m_ratio_h_dx) ? (comp_i - m_ratio_h_dx) : 0); i <= ((comp_i <= l - m_ratio_h_dx) ? (comp_i + m_ratio_h_dx) : l); i++) {
				for (j = ((comp_j >= m_ratio_h_dx) ? (comp_j - m_ratio_h_dx) : 0); j <= ((comp_j <= m - m_ratio_h_dx) ? (comp_j + m_ratio_h_dx) : m); j++) {
					for (r = ((comp_r >= m_ratio_h_dx) ? (comp_r - m_ratio_h_dx) : 0); r <= ((comp_r <= h - m_ratio_h_dx) ? (comp_r + m_ratio_h_dx) : h); r++) {
						for (int comp_index = 0; comp_index < m_grid_num[i][j][r]; comp_index++) {
							int neighbor_ID = idx_table[i][j][r][comp_index];
							double r_sqnorm = (Nodes[k].place - Nodes[neighbor_ID].place).squarenorm();
							if (r_sqnorm <= m_h * m_h) { Nodes[k].neighbor_index.push_back(neighbor_ID); }
						}
					}
				}
			}
		}

		//step1:计算粒子的Non-Pressure加速度：
#pragma omp parallel for private(k)
		for (k = 0; k < total_num; k++) {
			//Nodes[k].force_vis = P3(0.0, 0.0, 0.0); execute in step 0.
			//为什么只计算fluid的pressure加速度：因为边界条件不移动，没有必要；而刚体可以通过反作用力来直接给出
			if (Nodes[k].isWall()) { continue; }
			compute_vis_accelerate(k);
		}

		//step2:更新粒子的速度
#pragma omp parallel for private(k)
		for (k = 0; k < total_num; k++) {
			if (Nodes[k].isWall()) { continue; }
			Nodes[k].u += (Nodes[k].force_vis / m_mass + Nodes[k].acc_ext) * m_dt;
		}

		int iter_num = 0;
		double rho_error = 0.0;
		while (++iter_num <= 10) {
			//Step3: 根据光滑核函数计算粒子的插值密度：
#pragma omp parallel for private(k)
			for (k = 0; k < total_num; k++) {
				Compute_Rho(k);
			}

			//Step4: 计算粒子的压强：
#pragma omp parallel for private(k)
			for (k = 0; k < total_num; k++) {
				Nodes[k].force_press = P3(0.0, 0.0, 0.0);
				Nodes[k].p = m_KK * (pow(Nodes[k].rho / rho_L, (int)m_gamma) - 1.0);
				if (Nodes[k].p < 0.0) { Nodes[k].p = 0.0; }
			}

			//step5:计算粒子的Pressure加速度：
#pragma omp parallel for private(k)
			for (k = 0; k < total_num; k++) {
				//Nodes[k].acc_press = P3(0.0, 0.0, 0.0); execute in step 4.
				//为什么只计算fluid的pressure加速度：因为边界条件不移动，没有必要；而刚体可以通过反作用力来直接给出
				if (!Nodes[k].isFluid()) { continue; }
				compute_press_accelerate(k);
			}

			//step6:更新粒子的速度
#pragma omp parallel for private(k)
			for (k = 0; k < total_num; k++) {
				if (!Nodes[k].isFluid()) { continue; }
				Nodes[k].u += Nodes[k].force_press * (m_dt / m_mass);
			}

			//Check if we can break out the iteration:
			rho_error = 0.0;
#pragma omp parallel for private(k)
			for (k = 0; k < total_num; k++) {
				rho_error += (Nodes[k].rho - rho_L) * (Nodes[k].rho - rho_L);
			}
			rho_error = sqrt(rho_error);
			if (rho_error < 0.1) { break; }
		}
		std::cout << "rho_error: " << rho_error << std::endl;

		//step6.5:更新刚体运动
		/*
		m_rigidbody.linear_force = P3(0.0);
		m_rigidbody.torque = P3(0.0);
#pragma omp parallel for private(k,i)
		for (i = 0; i < m_rigidbody.shape.size(); i++) {
			k = m_rigidbody.shape[i];
			m_rigidbody.linear_force += Nodes[k].force_press + Nodes[k].force_vis;
			m_rigidbody.torque += (Nodes[k].place - m_rigidbody.m_sphere_center).cross_product(m_rigidbody.linear_force);
		}
		P3 sphere_center = m_rigidbody.b;
		m_rigidbody.update(m_dt);
		std::cout << "Sphere center: " << m_rigidbody.m_sphere_center << std::endl;
#pragma omp parallel for private(k,i)
		for (i = 0; i < m_rigidbody.shape.size(); i++) {
			k = m_rigidbody.shape[i];
			Nodes[k].u = m_rigidbody.v; // +m_rigidbody.angular_velocity.cross_product(Nodes[k].place - sphere_center);
			Nodes[k].place += Nodes[k].u * m_dt;
		}
		*/
		m_rigidbody.update(m_dt);
		double dTheta_dt = m_rigidbody.swirl_velocity;
		double dTheta = dTheta_dt * m_dt;
		double cosdTheta = cos(dTheta), sindTheta = sin(dTheta);
#pragma omp parallel for private(k,i)
		for (i = 0; i < m_rigidbody.shape.size(); i++) {
			k = m_rigidbody.shape[i];
			double local_x = Nodes[k].place.getx() - m_rigidbody.m_wheel_center.getx();
			double local_y = Nodes[k].place.gety() - m_rigidbody.m_wheel_center.gety();
			double r_norm = sqrt(local_x * local_x +local_y * local_y);
			double costheta = local_x / r_norm, sintheta = local_y / r_norm;
			//update:
			Nodes[k].setplace((costheta* cosdTheta - sintheta * sindTheta) * r_norm + m_rigidbody.m_wheel_center.getx(),
				(sintheta* cosdTheta + costheta * sindTheta) * r_norm + m_rigidbody.m_wheel_center.gety(),
				Nodes[k].place.getz());
		}


		//step7:更新粒子的位置：
#pragma omp parallel for private(k)
		for (k = 0; k < total_num; k++) {
			if (!Nodes[k].isFluid()) { continue; }
			Nodes[k].place += Nodes[k].u * m_dt;
			if (boundary_band[0] == 'P') { //Left
				if (Nodes[k].place.getx() < 0) { Nodes[k].place += P3(Lx, 0.0, 0.0); }
			}
			if (boundary_band[1] == 'P') { //Right
				if (Nodes[k].place.getx() > Lx) { Nodes[k].place += P3(-Lx, 0.0, 0.0); }
			}
			if (boundary_band[2] == 'P') { //Forward
				if (Nodes[k].place.gety() < 0) { Nodes[k].place += P3(0.0, Ly, 0.0); }
			}
			if (boundary_band[3] == 'P') { //Back
				if (Nodes[k].place.gety() > Ly) { Nodes[k].place += P3(0.0, -Ly, 0.0); }
			}
			if (boundary_band[4]== 'P') { //Down
				if (Nodes[k].place.getz() < 0) { Nodes[k].place += P3(0.0, 0.0, Lz); }
			}
			if (boundary_band[5] == 'P') { //Up
				if (Nodes[k].place.getz() > Lz) { Nodes[k].place += P3(0.0, 0.0, -Lz); }
			}
		}

		Vmax = 0.0;
#pragma omp parallel for private(k)
		for (k = 0; k < total_num; k++) {
			double u_maxAbsCoeff = Nodes[k].u.maxAbsCoeff();
			Vmax = Vmax > u_maxAbsCoeff ? Vmax : u_maxAbsCoeff;
		}
		//update m_dt:
		m_dt = m_lambda * m_h / max(Vmax, sqrt(m_stiffness));
		m_dt = m_lambda * m_h / sqrt(m_stiffness);
		std::cout << std::endl;

		return;
	}

	void Sph::dump_file(string file_name)
	{
		//本函数用于输出可以被软件OVITO可视化的文件 | This function can be used to dump a file which can be visualized by software OVITO
		ofstream fout(file_name);
		if (!fout)
		{
			std::cout << "(out_put function)Error! Can not write into this file." << std::endl;
			exit(-1);
		}

		fout << "Number of particles = " << total_num << std::endl;
		fout << "A = 1 Angstrom (basic length-scale)" << std::endl;

		fout << "H0(1,1) = " << Lx * 100 << " A" << std::endl;
		fout << "H0(1,2) = " << 0 << " A" << std::endl;
		fout << "H0(1,3) = " << 0 << " A" << std::endl;

		fout << "H0(2,1) = " << 0 << " A" << std::endl;
		fout << "H0(2,2) = " << Ly * 100 << " A" << std::endl;
		fout << "H0(2,3) = " << 0 << " A" << std::endl;

		fout << "H0(3,1) = " << 0 << " A" << std::endl;
		fout << "H0(3,2) = " << 0 << " A" << std::endl;
		fout << "H0(3,3) = " << Lz * 100 << " A" << std::endl;

		fout << ".NO_VELOCITY." << std::endl;
		fout << "entry_count = " << 3 << std::endl;

		for (int k = 0; k < total_num; k++){
			fout << m_mass << std::endl; // mass
			fout << Nodes[k].status << std::endl; // element type
			fout << Nodes[k].place.getx() * 100 << " " << Nodes[k].place.gety() * 100 << " " << Nodes[k].place.getz() * 100 << std::endl; // Info of place.
		}

		fout.close();

		return;
	}

	bool Sph::read_file(string file_name) {
		//本函数用于读取已经之前本程序输出过的文件
		ifstream fin(file_name);
		string get_str;
		if (!fin){
			std::cout << "(in_put function)Error! Can not read this file." << std::endl;
			return false;
		}
		for (int k = 0; k < 13; k++) {
			//line1: Number of particles = {$total_num}
			//line2: A = 1 Angstrom (basic length-scale)
			//line3: H0(1,1) = {$Lx * 100} A
			//line4: H0(1,2) = 0 A
			//line5: H0(1,3) = 0 A
			//line6: H0(2,1) = 0 A
			//line7: H0(2,2) = {$Ly * 100} A
			//line8: H0(2,3) = 0 A
			//line9: H0(3,1) = 0 A
			//line10: H0(3,2) = 0 A
			//line11: H0(3,3) = {$L * 100} A
			//line12: .NO_VELOCITY.
			//line13: entry_count = 3
			getline(fin, get_str);
		}
		int k = 0;
		double in_x, in_y, in_z;
		char status;
		while (k != total_num) {
			fin >> m_mass;
			fin >> status;
			Nodes[k].set_status(status);
			fin >> in_x >> in_y >> in_z;
			Nodes[k].setplace(in_x / 100.0, in_y / 100.0, in_z / 100.0);
			k++;
		}
		std::cout << "m_step: " << m_step << ", total_num: " << total_num << std::endl;
		return true;
	}

	void Sph::drawing_liquid(void)
	{
		GLint i = 0, j = 0, r = 0;

		glColor3f(0.0, 1.0, 0.0);
		//绘制流场外围
		glBegin(GL_LINE_STRIP);
		glVertex3i(0, 0, 0);
		glVertex3i(l, 0, 0);
		glVertex3i(l, m, 0);
		glVertex3i(0, m, 0);
		glVertex3i(0, 0, 0);
		glEnd();
		glBegin(GL_LINE_STRIP);
		glVertex3i(0, 0, h);
		glVertex3i(l, 0, h);
		glVertex3i(l, m, h);
		glVertex3i(0, m, h);
		glVertex3i(0, 0, h);
		glEnd();
		glBegin(GL_LINES);
		glVertex3i(0, 0, 0);
		glVertex3i(0, 0, h);
		glVertex3i(l, 0, 0);
		glVertex3i(l, 0, h);
		glVertex3i(0, m, 0);
		glVertex3i(0, m, h);
		glVertex3i(l, m, 0);
		glVertex3i(l, m, h);
		glEnd();
		//绘制流场外围结束

		GLfloat* vertices = nullptr;
		vertices = (GLfloat*)malloc(3 * total_num * (sizeof(GLfloat)));
		for (i = total_num - 1; i >= 0; i--)
		{
			glColor3f(0.5f, 0.5f, 1.0f);
			if (Nodes[i].isWall()) {
				continue;
			}
			if (Nodes[i].isRigidBody()) {
				glColor3f(1.0f, 1.0f, 0.0f);
			}
			//if (Nodes[i].place[0] < 0 || Nodes[i].place[0] > Lx) { continue; }
			//if (Nodes[i].place[1] < 0 || Nodes[i].place[1] > Ly) { continue; }
			//if (Nodes[i].place[2] < 0 || Nodes[i].place[2] > Lz) { continue; }
			glPushMatrix();
			glTranslated(Nodes[i].place.getx() / dx, Nodes[i].place.gety() / dy, Nodes[i].place.getz() / dz);
			glutSolidSphere(radius, 6, 6);
			glPopMatrix();

			//glVertex3f(Nodes[i].place[0] / dx, Nodes[i].place[1] / dy, Nodes[i].place[2] / dz);

			//*(vertices + 3 * i + 0) = Nodes[i].place[0] / dx;
			//*(vertices + 3 * i + 1) = Nodes[i].place[1] / dy;
			//*(vertices + 3 * i + 2) = Nodes[i].place[2] / dz;
		}

		//if(m_openglcontext != nullptr)
			//m_openglcontext->VBO_bind<GLfloat>(vertices);

		return;
	}

	void Sph::step() {
		double final_error = 1.0;
		record_velocity();
		drawing_liquid();
		compute();
		if (m_step++ % Nwri == 0)
		{
			final_error = comp_error();
			std::cout << "=====================" << std::endl;
			std::cout << "Step: " << m_step << ";  error: " << final_error << std::endl;
			char str_step[10];
			itoa(m_step, str_step, 10);
			if (if_dump) { dump_file(string("data") + str_step + ".cfg"); }
		}
	}

	void Sph::OnlyReadFileStep() {
		drawing_liquid();
		if (m_step++ % Nwri == 0)
		{
			char str_step[10];
			itoa(m_step, str_step, 10);
			if (!read_file(string("data") + str_step + ".cfg")) {
				m_step = 0;
			}
		}
	}
}