//
// main.cpp
//
// MD simulation code for water molecules using IPS/Tree method
//
// Kentaro NOMURA 2017/mm/dd
//
// Known problem
//

#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>

#include <string>
#include <vector>

#include "unit.h"
#include "pdb_manager.h"
#include "constraint.h"
#include "water_params.h"
#include "user_defined_class.h"

/* Cutoff functions */

const PS::F64 alpha[10] = {
  0.19578,
  0.0125143224110408,
  -0.603493863454666,
  11.7355819865242,
  -96.296895305654,
  216.649868508398,
  -197.409191110696,
  59.9544311773618,
  13.9564907382725,
  -8.66620089071555
  };

const int target = 0;

struct CalcForceEpEp{
  const PS::F64 rc_lj,rc_cl;

  CalcForceEpEp(const PS::F64 _rc_lj,const PS::F64 _rc_cl)
    : rc_lj(_rc_lj), rc_cl(_rc_cl){}

  template<class EPJ>
  void operator () (const EP * ep_i,
		    const PS::S32 n_ip,
		    const EPJ * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    const PS::F64 rc_lj_sq = rc_lj*rc_lj;
    const PS::F64 rc_cl_sq = rc_cl*rc_cl;
    const PS::F64 rc_inv = 1.0 / rc_cl;
    const PS::F64 rc2_inv = rc_inv * rc_inv;
    const PS::F64 rc3_inv = rc2_inv * rc_inv;

    const PS::F64 bound_ele_pot = rc_inv*1.296557;

    for(PS::S32 i=0; i<n_ip; i++){
      const PS::F64vec ri = ep_i[i].pos;
      //PS::F64vec ai = 0.0;
      PS::F64vec flj = 0.0;
      PS::F64vec fcl = 0.0;
      PS::F64 poti = 0.0;
      const PS::S64 idi = ep_i[i].id;
      const PS::F64 sigma_i   = ep_i[i].sigma;
      const PS::F64 epsilon_i = ep_i[i].epsilon;
      const PS::F64 charge_i  = ep_i[i].charge;
      for(PS::S32 j=0; j<n_jp; j++){
	PS::F64vec rij = ri - ep_j[j].pos;

	const PS::F64 r2 = rij*rij;
	if((idi/3) == (ep_j[j].id/3) && r2 < rc_lj_sq) continue;
	const PS::F64 r2_inv = 1.0/r2;
 	if(0.0 < r2 && r2 <= rc_cl_sq){
	  const PS::F64 r = sqrt(r2);
	  const PS::F64 r_inv = r2_inv * r;
	  const PS::F64 qq = charge_i * ep_j[j].getCharge();
	  PS::F64 rc = r * rc_inv;
	  PS::F64 r2c2 = rc*rc;
	  PS::F64 coef = (rc*rc - alpha[0]*alpha[0]);
	  PS::F64 coef2 = coef*coef;
	  const PS::F64 utmp = r2c2 * (alpha[1] + r2c2*
				       (alpha[2] + r2c2*
					(alpha[3] + r2c2*
					 (alpha[4] + r2c2*
					  (alpha[5] + r2c2*
					   (alpha[6] + r2c2*
					    (alpha[7] + r2c2*
					     (alpha[8] + (alpha[9] * r2c2)))))))));
	  const PS::F64 ftmp = (2.0*alpha[1] + r2c2*
				(4.0*alpha[2] + r2c2*
				 (6.0*alpha[3] + r2c2*
				  (8.0*alpha[4] + r2c2*
				   (10.0*alpha[5] + r2c2*
				    (12.0*alpha[6] + r2c2*
				     (14.0*alpha[7] + r2c2*
				      (16.0*alpha[8] + (18.0*alpha[9] * r2c2)))))))));
	  poti += qq * (r_inv - 0.5*rc_inv * coef * coef2 * utmp - bound_ele_pot);
	  fcl  += qq * (r2_inv*r_inv + 0.5*rc3_inv*coef2*(6.0*utmp + coef*ftmp))*rij;
#if 0
	  if(idi==target){
	    const PS::S32 type = ep_j[j].charge < 0.0 ? 0 : 1;
	    printf("%d %d %e %e %e %e %e\n",ep_j[j].id,type,rij.x,rij.y,rij.z,r,qq);
	  }
#endif
	}
	if (0.0 < r2 && r2 <= rc_lj_sq) {
	  const PS::F64 sigma    = 0.5*(sigma_i + ep_j[j].sigma);
	  const PS::F64 epsilon  = 4.0*sqrt(epsilon_i * ep_j[j].epsilon);
	  const PS::F64 rs2_inv  = sigma*sigma * r2_inv;
	  const PS::F64 rs6_inv  = rs2_inv * rs2_inv * rs2_inv;
	  const PS::F64 rs12_inv = rs6_inv * rs6_inv;
	  poti += epsilon * (rs12_inv - rs6_inv);
	  flj  += (epsilon * (12.0*rs12_inv - 6.0*rs6_inv) * r2_inv) * rij;
	}
      }
      force[i].f_lj += flj;
      force[i].f_coulomb += fcl;
      force[i].pot += 0.5*poti;
    }
  }
};

struct CalcForceEpSp{
  const PS::F64 rc_cl;
  CalcForceEpSp(const PS::F64 _rc_cl)
    : rc_cl(_rc_cl){}

  template <class SPJ>
  void operator () (const EP * ep_i,
		    const PS::S32 n_ip,
		    const SPJ * ep_j,
		    const PS::S32 n_jp,
		    Force * force){
    const PS::F64 rc_cl_sq = rc_cl*rc_cl;
    const PS::F64 rc_inv = 1.0 / rc_cl;
    const PS::F64 rc2_inv = rc_inv * rc_inv;
    const PS::F64 rc3_inv = rc2_inv * rc_inv;

    const PS::F64 bound_ele_pot = rc_inv*1.296557;

    for(PS::S32 i=0; i<n_ip; i++){
      const PS::F64vec ri = ep_i[i].pos;
      //PS::F64vec ai = 0.0;
      PS::F64vec flj = 0.0;
      PS::F64vec fcl = 0.0;
      PS::F64 poti = 0.0;
      const PS::S64 idi = ep_i[i].id;
      const PS::F64 charge_i  = ep_i[i].charge;
      for(PS::S32 j=0; j<n_jp; j++){
	PS::F64vec rij = ri - ep_j[j].getPos();
#if 0
	if(idi==target && rij*rij <= 900.){
	  int type = rij*rij <= rc_cl_sq ? 6 : 7;
	  printf("%d %d %e %e %e\n",j,type,rij.x,rij.y,rij.z);
	}
#endif
	// pseudo particles with positive charge
	for(int k=0;k<3;k++){
	  PS::F64vec rijk = rij - ep_j[j].ppp[k].pos;
	  const PS::F64 r2 = rijk*rijk;
	  if(r2 < 1.0 && ep_j[j].ppp[k].mass != 0.0){
	    printf("PPP is too close to epi: %d %d %lf %lf %lf : %lf\n",ep_i[i].id,j,rijk.x,rijk.y,rijk.z,ep_j[j].ppp[k].mass);
	    /*
	    printf("Distnce from GC: %lf\n",sqrt(ep_j[j].ppp[k].pos*ep_j[j].ppp[k].pos));
	    PS::F64vec tmp = ep_j[j].pos - ep_j[j].positive.pos;
 	    printf("Distnce between GC and center: %lf\n",sqrt(tmp*tmp));
	    //*/
	    ep_j[j].DebugPrint();
	  }
	  if(0.0 < r2 && r2 <= rc_cl_sq){
	    const PS::F64 r2_inv = 1.0/r2;
	    const PS::F64 r = sqrt(r2);
	    const PS::F64 r_inv = r2_inv * r;
	    const PS::F64 qq = charge_i * ep_j[j].ppp[k].mass;
	    PS::F64 rc = r * rc_inv;
	    PS::F64 r2c2 = rc*rc;
	    PS::F64 coef = (rc*rc - alpha[0]*alpha[0]);
	    PS::F64 coef2 = coef*coef;
	    const PS::F64 utmp = r2c2 * (alpha[1] + r2c2*
					 (alpha[2] + r2c2*
					  (alpha[3] + r2c2*
					   (alpha[4] + r2c2*
					    (alpha[5] + r2c2*
					     (alpha[6] + r2c2*
					      (alpha[7] + r2c2*
					       (alpha[8] + (alpha[9] * r2c2)))))))));
	    const PS::F64 ftmp = (2.0*alpha[1] + r2c2*
				  (4.0*alpha[2] + r2c2*
				   (6.0*alpha[3] + r2c2*
				    (8.0*alpha[4] + r2c2*
				     (10.0*alpha[5] + r2c2*
				      (12.0*alpha[6] + r2c2*
				       (14.0*alpha[7] + r2c2*
					(16.0*alpha[8] + (18.0*alpha[9] * r2c2)))))))));
	    poti += qq * (r_inv - 0.5*rc_inv * coef * coef2 * utmp - bound_ele_pot);
	    fcl  += qq * (r2_inv*r_inv + 0.5*rc3_inv*coef2*(6.0*utmp + coef*ftmp))*rijk;
	  }
	}
	for(int k=0;k<3;k++){
	  PS::F64vec rijk = rij - ep_j[j].ppn[k].pos;
	  const PS::F64 r2 = rijk*rijk;
	  const PS::F64 r2_inv = 1.0/r2;
	  if(0.0 < r2 && r2 <= rc_cl_sq){
	    const PS::F64 r = sqrt(r2);
	    const PS::F64 r_inv = r2_inv * r;
	    const PS::F64 qq = charge_i * ep_j[j].ppn[k].mass;
	    PS::F64 rc = r * rc_inv;
	    PS::F64 r2c2 = rc*rc;
	    PS::F64 coef = (rc*rc - alpha[0]*alpha[0]);
	    PS::F64 coef2 = coef*coef;
	    const PS::F64 utmp = r2c2 * (alpha[1] + r2c2*
					 (alpha[2] + r2c2*
					  (alpha[3] + r2c2*
					   (alpha[4] + r2c2*
					    (alpha[5] + r2c2*
					     (alpha[6] + r2c2*
					      (alpha[7] + r2c2*
					       (alpha[8] + (alpha[9] * r2c2)))))))));
	    const PS::F64 ftmp = (2.0*alpha[1] + r2c2*
				  (4.0*alpha[2] + r2c2*
				   (6.0*alpha[3] + r2c2*
				    (8.0*alpha[4] + r2c2*
				     (10.0*alpha[5] + r2c2*
				      (12.0*alpha[6] + r2c2*
				       (14.0*alpha[7] + r2c2*
					(16.0*alpha[8] + (18.0*alpha[9] * r2c2)))))))));
	    poti += qq * (r_inv - 0.5*rc_inv * coef * coef2 * utmp - bound_ele_pot);
	    fcl  += qq * (r2_inv*r_inv + 0.5*rc3_inv*coef2*(6.0*utmp + coef*ftmp))*rijk;
	  }
	}
      }
      force[i].f_coulomb += fcl;
      force[i].pot += 0.5*poti;
    }
  }
};

//#define IPS_TREE_FORCE_ERROR_CHECK

class ForceCalculator {
public:
#ifdef IPS_TREE_FORCE_ERROR_CHECK
  PS::TreeForForceShort<Force, EP, EP>::Scatter tree_short;
#endif
  PS::TreeForForceLong
  <Force, EP, EP,
   MomentQuadrupoleSixPseudoParticleCutoff,
   SPJQuadrupoleSixPseudoParticleCutoff>::WithCutoff tree_long;

  // Methods for water simulation
  template <class Tpsys>
  void initialize(const Tpsys &system,const PS::F64 theta) {
    PS::S32 numPtclLocal = system.getNumberOfParticleLocal();
    PS::U64 ntot = 3 * numPtclLocal;
#ifdef IPS_TREE_FORCE_ERROR_CHECK
    tree_short.initialize(ntot,0.0);
#endif
    tree_long.initialize(ntot,theta);
  };

  template <class Tpsys,class Tdinfo>
  void operator ()(Tpsys &system,Tdinfo &dinfo){
    // cut off radii settings
    const PS::F64 rcut_lj = 4.0*SIGMA_OXY;
    const PS::F64 rcut_cl = 28.0; // 28 angstrom

    // Reset potential and accelerations
    for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      system[i].pot = 0.0;
      system[i].acc = 0.0;
      // search radius must be larger than rcut_cl and shold have enough margin
      system[i].search_radius = 35.0;
    }
    tree_long.calcForceAll(CalcForceEpEp(rcut_lj,rcut_cl),
			   CalcForceEpSp(rcut_cl),
			   system, dinfo);
    for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      Force result = tree_long.getForce(i);
      system[i].pot += result.pot;
      system[i].acc += result.f_lj + result.f_coulomb;
    }
#ifdef IPS_TREE_FORCE_ERROR_CHECK
    tree_short.calcForceAll(CalcForceEpEp(rcut_lj,rcut_cl),
			    system, dinfo);
    PS::F64 error = 0.0;
    for (PS::S32 i=0; i<system.getNumberOfParticleLocal(); i++) {
      Force exact = tree_short.getForce(i);
      const PS::F64vec acc = exact.f_lj + exact.f_coulomb;
      const PS::F64vec diff = acc - system[i].acc;
      printf("%d x %lf %lf\n",i,acc.x,fabs(diff.x/acc.x));
      printf("%d y %lf %lf\n",i,acc.y,fabs(diff.y/acc.y));
      printf("%d z %lf %lf\n",i,acc.z,fabs(diff.z/acc.z));

      error += (diff.x/acc.x)*(diff.x/acc.x);
      error += (diff.y/acc.y)*(diff.y/acc.y);
      error += (diff.z/acc.z)*(diff.z/acc.z);
    }
    error /= 3.0 * system.getNumberOfParticleLocal();
    error = sqrt(error);
    fprintf(stderr,"error= %e\n",error);
#endif
    assert(system.getNumberOfParticleLocal()%3 == 0);
  }
};

template<class Tpsys>
void ReadPDBFile(Tpsys &psys,
		 const std::string filename){
  FileManager<PDBData> pdb;
  std::ifstream ifs_test(filename);
  assert(ifs_test.is_open());

  int natom = 0;
  PS::F64vec pos_min(1e32,1e32,1e32);
  for(std::string line;std::getline(ifs_test,line);){
    PDBData data = pdb.ReadLine(line);
    if(pos_min.x > data.coord[0]) pos_min.x = data.coord[0];
    if(pos_min.y > data.coord[1]) pos_min.y = data.coord[1];
    if(pos_min.z > data.coord[2]) pos_min.z = data.coord[2];
    if(data.record == "HETATM") natom++;
  }
  //std::cerr << "Number of atom in PDB file is " << natom << std::endl;
  psys.setNumberOfParticleLocal(natom);
  PS::MTTS mt;
  mt.init_genrand(PS::Comm::getRank()*PS::Comm::getNumberOfThread()+PS::Comm::getThreadNum());

  std::ifstream ifs(filename);
  assert(ifs.is_open());

  fprintf(stderr,"charge is OXY: %lf, HYD: %lf\n",CHARGE_OXY,CHARGE_HYD);
  int count = 0;
  for(std::string line;std::getline(ifs,line);){
    PDBData data = pdb.ReadLine(line);
    //std::cout << data << std::endl;
    if(data.record == "HETATM"){
      if(count%3==0) assert(data.resName == "O  ");
      else           assert(data.resName == "H  ");

      psys[count].pos.x = data.coord[0] - pos_min.x;
      psys[count].pos.y = data.coord[1] - pos_min.y;
      psys[count].pos.z = data.coord[2] - pos_min.z;

      const double v_max = 0.1;
      do {
	psys[count].vel[0] = (2.*mt.genrand_res53() - 1.) * v_max;
	psys[count].vel[1] = (2.*mt.genrand_res53() - 1.) * v_max;
	psys[count].vel[2] = (2.*mt.genrand_res53() - 1.) * v_max;
      }while(psys[count].vel * psys[count].vel >= v_max * v_max);

      psys[count].id      = count;
      psys[count].type    = count%3==0 ? 0 : 1;
      psys[count].mass    = (count%3==0 ? MASS_OXY    : MASS_HYD);
      psys[count].sigma   = (count%3==0 ? SIGMA_OXY   : SIGMA_HYD);
      psys[count].epsilon = (count%3==0 ? EPSILON_OXY : EPSILON_HYD);
      psys[count].charge  = (count%3==0 ? CHARGE_OXY  : CHARGE_HYD);
      count++;
    }
  }
  assert(count == natom);
}

template<class Tpsys>
void MakeIceLattice(Tpsys &psys,
		    PS::F64vec &cell_size){
  ReadPDBFile(psys,"ice_unit.pdb");

  //const PS::S32 nx=16,ny=16,nz=8;
  const PS::S32 nx=8,ny=8,nz=4;
  //const PS::S32 nx=2,ny=2,nz=2;
  const PS::F64vec unit(4.511,4.511,7.315);
  //const PS::F64vec unit(4.48,4.48,7.31);
  cell_size.x = unit.x * nx * sin(M_PI*60./180.);
  cell_size.y = unit.y * ny;
  cell_size.z = unit.z * nz;
  fprintf(stdout,"cellsize: %lf %lf %lf\n",cell_size.x,cell_size.y,cell_size.z);

  const int nunit = psys.getNumberOfParticleLocal();
  const int natom = nunit*nx*ny*nz;
  psys.setNumberOfParticleLocal(natom);

  fprintf(stdout,"density is %lf kg m^3\n", unit_density*(MASS_OXY+MASS_HYD*2)*(natom/3.0) / (cell_size.x*cell_size.y*cell_size.z) );

  int count = nunit;
  for(int x=0;x<nx;x++){
    for(int y=0;y<ny;y++){
      for(int z=0;z<nz;z++){
	if(x==0 && y==0 && z==0)continue;
	for(int i=0;i<12;i++){
	  psys[count] = psys[i];

	  psys[count].pos.x += unit.x*x*sin(M_PI*60./180.);
	  psys[count].pos.y += unit.x*x*cos(M_PI*60./180.) + unit.y*y;
	  psys[count].pos.z += unit.z*z;
	  psys[count].id = count;
	  count++;
	}
      }}}
  assert(count == natom);
}

template<class Tpsys,class Tdinfo>
void PeriodicBoundaryCondition(Tpsys & system,const Tdinfo &dinfo){
  PS::S32 n = system.getNumberOfParticleLocal();
  const PS::F64ort domain = dinfo.getPosRootDomain();
  const PS::F64vec length = domain.getFullLength();

  for(int i=0; i<n; i++){
    for(int k=0;k<3;k++){
      if(system[i].gpos.x <  domain.low_.x){
	system[i].pos.x  += length.x;
	system[i].gpos.x += length.x;
      }
      if(system[i].gpos.x >= domain.high_.x){
	system[i].pos.x  -= length.x;
	system[i].gpos.x -= length.x;
      }
      if(system[i].gpos.y <  domain.low_.y){
	system[i].pos.y  += length.y;
	system[i].gpos.y += length.y;
      }
      if(system[i].gpos.y >= domain.high_.y){
	system[i].pos.y  -= length.y;
	system[i].gpos.y -= length.y;
      }
      if(system[i].gpos.z <  domain.low_.z){
	system[i].pos.z  += length.z;
	system[i].gpos.z += length.z;
      }
      if(system[i].gpos.z >= domain.high_.z){
	system[i].pos.z  -= length.z;
	system[i].gpos.z -= length.z;
      }
    }
  }
}

template<class Tpsys>
void RemoveTotalMomentum(Tpsys &system){
  const PS::S32 n_loc = system.getNumberOfParticleLocal();
  PS::F64vec cm_vel_loc = 0.0;
  PS::F64    cm_mass_loc = 0.0;
  for(int i=0; i<n_loc; i++){
    cm_vel_loc  += system[i].mass * system[i].vel;
    cm_mass_loc += system[i].mass;
  }
  PS::F64vec cm_vel;
  PS::F64    cm_mass;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.x, &cm_vel.x, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.y, &cm_vel.y, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_vel_loc.z, &cm_vel.z, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&cm_mass_loc,  &cm_mass,  1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  cm_vel = cm_vel_loc;
  cm_mass = cm_mass_loc;
#endif
  cm_vel /= cm_mass;
  for(int i=0; i<n_loc; i++){
    system[i].vel -= cm_vel;
  }
}

template<class Tpsys>
void ScaleVelocity(Tpsys & system,const PS::F64 T){
  const PS::S32 natom_local = system.getNumberOfParticleLocal();
  PS::F64 ekin_loc = 0.0;
  for(PS::S32 i=0; i<natom_local; i++){
    ekin_loc += system[i].mass * system[i].vel * system[i].vel;
  }
  ekin_loc *= 0.5;
  PS::S32 natom;
  PS::F64 ekin;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&ekin_loc, &ekin, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&natom_local, &natom, 1, PS::GetDataType<PS::S32>(), MPI::SUM);
#else
  ekin = ekin_loc;
  natom = natom_local;
#endif
  const PS::F64 scaler = sqrt(1.5*natom*T / ekin);
  for(PS::S32 i=0;i<natom_local;i++) system[i].vel *= scaler;

  RemoveTotalMomentum(system);
}

template<class Tpsys>
void CalcEnergy(const Tpsys & system,
                PS::F64 & etot,
                PS::F64 & ekin,
                PS::F64 & epot,
                const bool clear=true){
  if(clear){
    etot = ekin = epot = 0.0;
  }
  PS::F64 ekin_loc = 0.0;
  PS::F64 epot_loc = 0.0;
  const PS::S32 nbody = system.getNumberOfParticleLocal();
  for(PS::S32 i=0; i<nbody; i++){
    ekin_loc += system[i].mass * (system[i].vel * system[i].vel);
    epot_loc += system[i].pot;
  }
  ekin_loc *= 0.5;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
  MPI::COMM_WORLD.Allreduce(&epot_loc, &epot, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
  MPI::COMM_WORLD.Allreduce(&ekin_loc, &ekin, 1, PS::GetDataType<PS::F64>(), MPI::SUM);
#else
  epot = epot_loc;
  ekin = ekin_loc;
#endif
  etot = epot + ekin;
}

template<class Tpsys>
void SetGravityCenterForAtoms(Tpsys & system){
  const PS::S32 n = system.getNumberOfParticleLocal();
  const unsigned int nmol = n/3;
  for(unsigned int i=0;i<nmol;i++){
    PS::F64vec3 gpos = 0.0;
    PS::F64 mass = 0.0;
    for(int j=0;j<3;j++){
      gpos += system[3*i+j].mass * system[3*i+j].pos;
      mass += system[3*i+j].mass;
    }
    gpos /= mass;
    for(int j=0;j<3;j++) system[3*i+j].gpos = gpos;
  }
}


int main(int argc, char *argv[]){
  PS::Initialize(argc, argv);

  PS::F64 Tbegin = PS::GetWtime();
  std::cout<<std::setprecision(15);
  std::cerr<<std::setprecision(15);
  PS::F64 theta = 0.2;
  const PS::S32 n_leaf_limit = 32;

  long long int nmol  = 500;
  PS::F64 density     = 0.997 * 1e3 / unit_density; // 1000 kg/m^3
  PS::F64 temperature = 80.0 / unit_temp;        // 300 K

  PS::F64 dt   = 2.0 * 1e-15 / unit_time; // 2 fs

  PS::S32 n_group_limit = 64;
  PS::S32 nstep      = 1000;
  PS::S32 nstep_eq   = 1000;
  PS::S32 nstep_snp  = 100;
  PS::S32 nstep_diag = 100;

  //char sinput[1024];
  char dir_name[1024];
  int c;
  sprintf(dir_name,"./result");
  while((c=getopt(argc,argv,"o:T:s:e:S:D:t:n:h")) != -1){
    switch(c){
    case 'o':
      sprintf(dir_name,optarg);
      break;
    case 'T':
      temperature = atof(optarg) / unit_temp;
      std::cerr<<"temperature="<< temperature << " = " << temperature*unit_temp << " K" << std::endl;
      break;
    case 's':
      nstep = atoi(optarg);
      std::cerr<<"nstep="<<nstep<<std::endl;
      break;
    case 'e':
      nstep_eq = atoi(optarg);
      std::cerr<<"nstep_eq="<<nstep_eq<<std::endl;
      break;
    case 'S':
      nstep_snp = atoi(optarg);
      std::cerr<<"nstep_snp="<<nstep_snp<<std::endl;
      break;
    case 'D':
      nstep_diag = atoi(optarg);
      std::cerr<<"nstep_diag="<<nstep_diag<<std::endl;
      break;
    case 't':
      theta = atof(optarg);
      std::cerr<<"theta="<<dt<<std::endl;
      break;
    case 'n':
      n_group_limit = atoi(optarg);
      std::cerr<<"n_group_limit="<<n_group_limit<<std::endl;
      break;
    case 'h':
      std::cerr<<"N: n_tot (default: 1000)"<<std::endl;
      std::cerr<<"d: number density (default: 1000 kg/m^3)"<<std::endl;
      std::cerr<<"T: temperature (default: 300 K)"<<std::endl;
      std::cerr<<"s: number of steps (default: 1000)"<<std::endl;
      std::cerr<<"S: time step for snapshot(default: 100)"<<std::endl;
      std::cerr<<"D: time step for diag(default: 100)"<<std::endl;
      std::cerr<<"e: number of steps for equilibration(default: 1000)"<<std::endl;
      std::cerr<<"o: dir name of output (default: ./result)"<<std::endl;
      std::cerr<<"t: parameter for accuracy of force calculation (default: 0.2)"<<std::endl;
      std::cerr<<"c: cutoff length (default: 15.0 angstrom)"<<std::endl;
      std::cerr<<"n: n_group_limit (default: 64.0)"<<std::endl;
      return 0;
    }
  }
  long long int n_tot = nmol;

  struct stat st;
  if(stat(dir_name, &st) != 0) {
    PS::S32 rank = PS::Comm::getRank();
    PS::S32 ret_loc, ret=0;
    if(rank == 0)
      ret_loc = mkdir(dir_name, 0777);
    PS::Comm::broadcast(&ret_loc, ret);
    if(ret == 0) {
      if(rank == 0)
	fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
    } else {
      fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
      PS::Abort();
      exit(0);
    }
  }

  std::ofstream fout_eng;
  std::ofstream fout_tcal;
  if(PS::Comm::getRank() == 0){
    char sout_de[1024];
    char sout_tcal[1024];
    sprintf(sout_de, "%s/t-de.dat", dir_name);
    sprintf(sout_tcal, "%s/t-tcal.dat", dir_name);
    std::cerr<<sout_de<<std::endl;
    std::cerr<<sout_tcal<<std::endl;
    fout_eng.open(sout_de);
    fout_tcal.open(sout_tcal);
  }

  PS::ParticleSystem<FP> system_water;
  system_water.initialize();

  PS::S32 n_grav_glb = n_tot;
  PS::F64vec box_size;
  MakeIceLattice(system_water,box_size);

  const PS::F64 coef_ema = 0.3;
  PS::DomainInfo dinfo;
  dinfo.initialize(coef_ema);
  dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
  dinfo.setPosRootDomain(PS::F64vec(0.0,0.0,0.0),
			 PS::F64vec(box_size.x,box_size.y,box_size.z));

  SetGravityCenterForAtoms(system_water);
  PeriodicBoundaryCondition(system_water,dinfo);

  Constraint constraint;
  for(int i=0;i<system_water.getNumberOfParticleLocal()/3;i++){
    for(int j=0;j<3;j++) system_water[3*i+j].KeepCurrentPos();
    constraint.Shake(&system_water[3*i],box_size,dt);
    constraint.Rattle(&system_water[3*i],box_size,dt);
  }

  dinfo.collectSampleParticle(system_water);
  dinfo.decomposeDomain();
  system_water.exchangeParticle(dinfo);

  ForceCalculator ips_tree;
  ips_tree.initialize(system_water,theta);
  ips_tree(system_water,dinfo);

  PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
  ScaleVelocity(system_water,temperature);
  CalcEnergy(system_water, Etot0, Ekin0, Epot0);
  if(PS::Comm::getRank() == 0) fprintf(stdout,"Etot = %lf, Epot = %lf, Ekin = %lf\n",Etot0,Epot0,Ekin0);
  PS::F64 Tloop = 0.0;

  PS::S32 snp_id = 0;
  PS::S32 time_snp = 0;
  PS::S32 time_diag = 0;
  PS::F64 time_sys = -dt * nstep_eq;
  bool isInitialized = false;

  PS::F64 Epot_ave = 0.0, Ekin_ave = 0.0;
  int n_loc = system_water.getNumberOfParticleLocal();
  const PS::S32 nstep_em = 0.5 * nstep_eq;

  // Main roop
  for(int s=-nstep_eq;s<nstep;s++){
    if(s < 0) ScaleVelocity(system_water,temperature);
    if(s%1000==0) RemoveTotalMomentum(system_water);

    if(s == time_snp){
      FileHeader header(box_size);
      char filename[256];
      sprintf(filename, "%s/%04d.cdv", dir_name, snp_id++);
      system_water.writeParticleAscii(filename, header);
      time_snp += nstep_snp;
    }

    if(!isInitialized && s == 0){
      CalcEnergy(system_water, Etot0, Ekin0, Epot0);
      if(PS::Comm::getRank() == 0) printf("Etot0 = %lf, Epot0 = %lf, Ekin0 = %lf\n",Etot0,Epot0,Ekin0);
      isInitialized = true;
    }

    for(int i=0;i<system_water.getNumberOfParticleLocal()/3;i++){
      for(int j=0;j<3;j++){
	system_water[3*i+j].KeepCurrentPos();
	system_water[3*i+j].IntegrateVel(0.5*dt);
	system_water[3*i+j].IntegratePos(dt,box_size);
      }

      constraint.Shake(&system_water[3*i],box_size,dt);
    }
    SetGravityCenterForAtoms(system_water);
    PeriodicBoundaryCondition(system_water,dinfo);
    if(s%10 == 0){
      dinfo.collectSampleParticle(system_water);
      dinfo.decomposeDomain();
    }
    system_water.exchangeParticle(dinfo);
    assert(system_water.getNumberOfParticleLocal()%3 == 0);

    ips_tree(system_water, dinfo);

    for(int i=0;i<system_water.getNumberOfParticleLocal();i++){
      system_water[i].IntegrateVel(0.5*dt);
    }
    for(int i=0;i<system_water.getNumberOfParticleLocal()/3;i++){
      constraint.Rattle(&system_water[3*i],box_size,dt);
    }
    CalcEnergy(system_water, Etot1, Ekin1, Epot1);
    if(s>=0){
      Epot_ave += Epot1;
      Ekin_ave += Ekin1;
    }
    if(s == time_diag) {
      if(PS::Comm::getRank() == 0){
	fout_eng<<time_sys<<"   "<< " " << Epot1 << " " << Ekin1 << " " <<(Etot1-Etot0)/Etot0<<std::endl;
	fprintf(stderr, "%10.7f %lf %lf %+e\n",
		time_sys, Epot1, Ekin1, (Etot1 - Etot0) / Etot0);
	time_diag += nstep_diag;
      }
    }
    time_sys += dt;
  }
  PS::Finalize();
  return 0;
}
