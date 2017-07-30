#include <iostream>
#include <iomanip>
#include "particle_simulator.hpp"
#include "io_util.hpp"
#include "parameter.hpp"
#include "f_calculator.hpp"
#include <ctime>

static_assert(Parameter::bond_leng != 0.0, "bond_leng should not be 0.0");

constexpr int sys = 16;

std::array<PS::F64vec, sys> x;
std::array<PS::F64vec, sys> f_ref, f;

constexpr int nanglelist = sys - 2;
int anglelist[nanglelist][3];

constexpr int bond_n = sys;
const double cf_b = 12.0;
const double cf_s = 10.0;

constexpr int nbondlist = sys - 1;
int bondlist[nbondlist][2];

void clear_force() {
  for(int i = 0; i < sys; i++)
    f_ref[i] = f[i] = 0.0;
}

void generate_chain() {
  PS::MT::init_genrand(static_cast<unsigned long>(time(NULL)));
  for(int i = 0; i < sys; i++) {
    const PS::F64vec dR(PS::MT::genrand_real1() * 0.5,
			PS::MT::genrand_real1() * 0.5,
			PS::MT::genrand_real1() * 0.5);
    x[i] = i * dR;
  }
}

void make_angle_list() {
  for(int i = 0; i < nanglelist; i++) {
    anglelist[i][0] = i;
    anglelist[i][1] = i + 1;
    anglelist[i][2] = i + 2;
  }
}

void make_bond_list() {
  for(int i = 0; i < nbondlist; i++) {
    bondlist[i][0] = i;
    bondlist[i][1] = i + 1;
  }
}

void StoreBondForceWithARLaw(const PS::F64vec&	__restrict dr,
			     const PS::F64&	__restrict inv_dr,
			     PS::F64vec&	__restrict d_vir,
			     PS::F64&		__restrict d_lap,
			     PS::F64vec*	__restrict F)
{
  const PS::F64 cf_bond = cf_s * (inv_dr - Parameter::ibond);
  const PS::F64vec Fbond(cf_bond * dr.x, cf_bond * dr.y, cf_bond * dr.z);

  //NOTE: The value of virial is twice.
  d_vir.x += 2.0 * dr.x * Fbond.x;
  d_vir.y += 2.0 * dr.y * Fbond.y;
  d_vir.z += 2.0 * dr.z * Fbond.z;

  //NOTE: The value of lap is twice.
  d_lap += 2.0 * cf_s * (6.0 * Parameter::ibond - 4.0 * inv_dr);
  F[0] -= Fbond;
  F[1] += Fbond;
}

void StoreBendForceWithARLaw(const PS::F64vec*	__restrict dr,
			     const PS::F64*	__restrict inv_dr,
			     const PS::F64*	__restrict dist,
			     PS::F64vec&	__restrict d_vir,
			     PS::F64&		__restrict d_lap,
			     PS::F64vec*	__restrict F)
{
  const PS::F64	inv_dr_prod     = inv_dr[0] * inv_dr[1];
  const PS::F64	inv_dist[2]     = { inv_dr[0] * inv_dr[0],
                                inv_dr[1] * inv_dr[1] };
  const PS::F64	in_prod         = dr[0] * dr[1];
  const PS::F64	cf_bd           = cf_b * inv_dr_prod;
  const PS::F64       cf_crs[2]	= { in_prod * inv_dist[0],
                                    in_prod * inv_dist[1] };

  const PS::F64vec Ftb0(cf_bd * (dr[1].x - cf_crs[0] * dr[0].x),
                        cf_bd * (dr[1].y - cf_crs[0] * dr[0].y),
                        cf_bd * (dr[1].z - cf_crs[0] * dr[0].z));
  const PS::F64vec Ftb1(cf_bd * (dr[0].x - cf_crs[1] * dr[1].x),
                        cf_bd * (dr[0].y - cf_crs[1] * dr[1].y),
                        cf_bd * (dr[0].z - cf_crs[1] * dr[1].z));

  //NOTE: The value of virial is twice.
  d_vir.x += 2.0 * (dr[0].x * Ftb0.x + dr[1].x * Ftb1.x);
  d_vir.y += 2.0 * (dr[0].y * Ftb0.y + dr[1].y * Ftb1.y);
  d_vir.z += 2.0 * (dr[0].z * Ftb0.z + dr[1].z * Ftb1.z);

  //NOTE: The value of lap is twice.
  d_lap += 2.0 * 2.0 * cf_bd * inv_dist[0] * inv_dist[1] * ( in_prod * ( in_prod + 2.0 * (dist[0] + dist[1]) ) + dist[0] * dist[1]);

  F[0] -= Ftb0;
  F[1] += Ftb0 - Ftb1;
  F[2] += Ftb1;
}

void bondsangles() {
  PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
  PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];
  PS::F64vec d_vir;
  PS::F64 d_lap = 0.0;

  pos_buf[0] = x[0];
  pos_buf[1] = x[1];

  dr[0] = pos_buf[1] - pos_buf[0];
  dist2[0] = dr[0] * dr[0];
  inv_dr[0] = 1.0 / std::sqrt(dist2[0]);

  StoreBondForceWithARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0]);

  for(PS::U32 unit = 2; unit < bond_n; unit++) {
    pos_buf[unit]    = x[unit];
    dr[unit - 1]     = pos_buf[unit] - pos_buf[unit - 1];
    dist2[unit - 1]  = dr[unit - 1] * dr[unit - 1];
    inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);

    StoreBondForceWithARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1]);
    StoreBendForceWithARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2]);
  }

  //Store the sum of force.
  for(PS::U32 unit = 0; unit < bond_n; unit++)
    f[unit] += Fbb[unit];
}

void bonds() {
  PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
  PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];
  PS::F64vec d_vir;
  PS::F64 d_lap = 0.0;

  pos_buf[0] = x[0];
  pos_buf[1] = x[1];

  dr[0] = pos_buf[1] - pos_buf[0];
  dist2[0] = dr[0] * dr[0];
  inv_dr[0] = 1.0 / std::sqrt(dist2[0]);

  StoreBondForceWithARLaw(dr[0], inv_dr[0], d_vir, d_lap, &Fbb[0]);

  for(PS::U32 unit = 2; unit < bond_n; unit++) {
    pos_buf[unit] = x[unit];
    dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
    dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
    inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);

    StoreBondForceWithARLaw(dr[unit - 1], inv_dr[unit - 1], d_vir, d_lap, &Fbb[unit - 1]);
  }

  //Store the sum of force.
  for(PS::U32 unit = 0; unit < bond_n; unit++)
    f[unit] = Fbb[unit];
}

void angles() {
  PS::F64vec Fbb[bond_n], pos_buf[bond_n], dr[bond_n - 1];
  PS::F64  dist2[bond_n - 1], inv_dr[bond_n - 1];
  PS::F64vec d_vir;
  PS::F64 d_lap = 0.0;

  pos_buf[0] = x[0];
  pos_buf[1] = x[1];

  dr[0] = pos_buf[1] - pos_buf[0];
  dist2[0] = dr[0] * dr[0];
  inv_dr[0] = 1.0 / std::sqrt(dist2[0]);

  for(PS::U32 unit = 2; unit < bond_n; unit++) {
    pos_buf[unit] = x[unit];
    dr[unit - 1] = pos_buf[unit] - pos_buf[unit - 1];
    dist2[unit - 1] = dr[unit - 1] * dr[unit - 1];
    inv_dr[unit - 1] = 1.0 / std::sqrt(dist2[unit - 1]);

    StoreBendForceWithARLaw(&dr[unit - 2], &inv_dr[unit - 2], dist2, d_vir, d_lap, &Fbb[unit - 2]);
  }

  //Store the sum of force.
  for(PS::U32 unit = 0; unit < bond_n; unit++)
    f[unit] = Fbb[unit];
}

void lammps_bonds() {
  for (int n = 0; n < nbondlist; n++) {
    const int i1 = bondlist[n][0];
    const int i2 = bondlist[n][1];

    const double delx = x[i1][0] - x[i2][0];
    const double dely = x[i1][1] - x[i2][1];
    const double delz = x[i1][2] - x[i2][2];

    const double rsq = delx*delx + dely*dely + delz*delz;
    const double r  = sqrt(rsq);
    const double dr = r - Parameter::bond_leng;
    const double rk = cf_s * dr;

    // force & energy
    const double fbond = -rk / r * Parameter::ibond;

    // apply force to each of 2 atoms
    f_ref[i1][0] += delx*fbond;
    f_ref[i1][1] += dely*fbond;
    f_ref[i1][2] += delz*fbond;

    f_ref[i2][0] -= delx*fbond;
    f_ref[i2][1] -= dely*fbond;
    f_ref[i2][2] -= delz*fbond;
  }
}

void lammps_angles() {
  double f1[3], f3[3];

  for (int n = 0; n < nanglelist; n++) {
    const int i1 = anglelist[n][0];
    const int i2 = anglelist[n][1];
    const int i3 = anglelist[n][2];

    // 1st bond

    const double delx1 = x[i1][0] - x[i2][0];
    const double dely1 = x[i1][1] - x[i2][1];
    const double delz1 = x[i1][2] - x[i2][2];

    const double rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    const double r1 = sqrt(rsq1);

    // 2nd bond

    const double delx2 = x[i3][0] - x[i2][0];
    const double dely2 = x[i3][1] - x[i2][1];
    const double delz2 = x[i3][2] - x[i2][2];

    const double rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    const double r2 = sqrt(rsq2);

    // c = cosine of angle

    double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // force & energy

    const double a11 = cf_b * c / rsq1;
    const double a12 = -cf_b / (r1*r2);
    const double a22 = cf_b * c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms
    f_ref[i1][0] += f1[0];
    f_ref[i1][1] += f1[1];
    f_ref[i1][2] += f1[2];

    f_ref[i2][0] -= f1[0] + f3[0];
    f_ref[i2][1] -= f1[1] + f3[1];
    f_ref[i2][2] -= f1[2] + f3[2];

    f_ref[i3][0] += f3[0];
    f_ref[i3][1] += f3[1];
    f_ref[i3][2] += f3[2];
  }
}

void check_error() {
  int err = 0;
  const double eps = 1.0e-13;
  for(int i = 0; i < sys; i++) {
    for(int j = 0; j < 3; j++) {
      if(fabs(f_ref[i][j] - f[i][j]) > eps) {
        std::cout << std::setprecision(15);
        std::cout << f_ref[i][j] << " " << f[i][j] << std::endl;
        err++;
      }
    }
  }
  if(err != 0)
    std::cerr << "error count is " << err << std::endl;
  else
    std::cerr << "Check PASS\n";
}

int main() {
  generate_chain();
  make_angle_list();
  make_bond_list();

  clear_force();
  angles();
  lammps_angles();
  check_error();

  clear_force();
  bonds();
  lammps_bonds();
  check_error();

  clear_force();
  lammps_bonds();
  lammps_angles();
  bondsangles();
  check_error();
}
