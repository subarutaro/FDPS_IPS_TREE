#ifndef H_USER_DEFINED_CLASS
#define H_USER_DEFINED_CLASS

#include "water_params.h"


class FileHeader{
public:
  const PS::F64vec cell_size;
  FileHeader(const PS::F64vec _cs) : cell_size(_cs){}
  PS::S32 readAscii(FILE * fp){
    return 0;
  }
  void writeAscii(FILE* fp) const{
    fprintf(fp,
	    "'box_sx=%lf,box_sy=%lf,box_sz=%lf,box_ex=%lf,box_ey=%lf,box_ez=%lf\n",
	    0.0,0.0,0.0,
	    cell_size.x,cell_size.y,cell_size.z);
  }
};

// =================
// Force class
// =================
class Force{
public:
  PS::F64vec3 f_lj;
  PS::F64vec3 f_coulomb;
  PS::F64 pot;
  void clear(){
    f_lj = 0.0;
    f_coulomb = 0.0;
    pot = 0.0;
  }
};

// =========================================
// Full particle and essential particle class
// =========================================
class FP{
public:
  PS::S64 id;
  PS::S64 type;

  PS::F64vec3 gpos;
  PS::F64vec3 pos;
  PS::F64vec3 vel;

  PS::F64vec3 acc;
  PS::F64vec3 force_long;

  PS::F64 mass;
  PS::F64 sigma;
  PS::F64 epsilon;
  PS::F64 charge;

  PS::F64 pot;

  PS::F64 search_radius;

  PS::F64vec prev_pos; // for RATTLE

  PS::F64 getRsearch() const {
    return search_radius;
  }
  PS::F64vec getPos() const {
    return gpos;
  }

  void copyFromForce(const Force& force){
    acc = force.f_lj + force.f_coulomb;
    pot = force.pot;
  }


  void writeAscii(FILE* fp) const{
    fprintf(fp, "%lld %lld %lf %lf %lf\n",
	    id, type, pos.x, pos.y, pos.z);
  }
  void readAscii(FILE* fp){
    fscanf(fp, "%lld %lld %lf %lf %lf\n",
	   &id, &type, &pos.x, &pos.y, &pos.z);
  }

  // for constraint
  void KeepCurrentPos(){
    prev_pos = pos;
  }
  void IntegratePos(const PS::F64 dt,const PS::F64vec cell_size){
    pos += vel*dt;
  }
  void IntegrateVel(const PS::F64 dth){
    vel += acc/mass*dth;
  }

  void EnergyMinimize(const PS::F64 dt){
    pos += acc/mass*dt;
  }
};

class EP{
public:
  PS::S64 id;
  PS::F64vec pos;

  PS::F64 sigma;
  PS::F64 epsilon;
  PS::F64 charge;

  PS::F64 search_radius;

  PS::F64    getRSearch() const {return search_radius;}
  PS::F64vec getPos()     const {return pos;}
  PS::F64    getCharge()  const {return charge;}
  void copyFromFP(const FP & fp){
    pos     = fp.pos;
    id      = fp.id;
    sigma   = fp.sigma;
    epsilon = fp.epsilon;
    charge  = fp.charge;
    search_radius = fp.search_radius;
  }

  void setPos(const PS::F64vec3 _pos){
    pos = _pos;
  }
};


// =================================
// Moment and Superparticle class for IPS/Tree
//==================================

#define N_PSEUDOPARTICLE 4

class PseudoParticle{
public:
  PS::F64    mass;
  PS::F64vec pos;

  void copy(const PseudoParticle &p){
    mass = p.mass;
    pos  = p.pos;
  }
  void init(){
    mass = 0.0;
    pos  = 0.0;
  }

  PS::F64vec getPos() const {
    return pos;
  }
  PS::F64 getCharge() const {
    return mass;
  }
  PS::F64 getRSearch() const {
    //return 35.0;
    return 40.0;
    //return 50.0;
  }
};

const PS::F64mat3asym mul(const PS::F64mat3asym& lhs,const PS::F64mat3asym& rhs){
  return PS::F64mat3asym(PS::F64vec(lhs.xx,lhs.xy,lhs.xz) * PS::F64vec(rhs.xx,rhs.yx,rhs.zx),
			 PS::F64vec(lhs.xx,lhs.xy,lhs.xz) * PS::F64vec(rhs.xy,rhs.yy,rhs.zy),
			 PS::F64vec(lhs.xx,lhs.xy,lhs.xz) * PS::F64vec(rhs.xz,rhs.yz,rhs.zz),

			 PS::F64vec(lhs.yx,lhs.yy,lhs.yz) * PS::F64vec(rhs.xx,rhs.yx,rhs.zx),
			 PS::F64vec(lhs.yx,lhs.yy,lhs.yz) * PS::F64vec(rhs.xy,rhs.yy,rhs.zy),
			 PS::F64vec(lhs.yx,lhs.yy,lhs.yz) * PS::F64vec(rhs.xz,rhs.yz,rhs.zz),

			 PS::F64vec(lhs.zx,lhs.zy,lhs.zz) * PS::F64vec(rhs.xx,rhs.yx,rhs.zx),
			 PS::F64vec(lhs.zx,lhs.zy,lhs.zz) * PS::F64vec(rhs.xy,rhs.yy,rhs.zy),
			 PS::F64vec(lhs.zx,lhs.zy,lhs.zz) * PS::F64vec(rhs.xz,rhs.yz,rhs.zz));
}

const PS::F64mat3asym trans(const PS::F64mat3asym& m){
  return PS::F64mat3asym(m.xx,m.yx,m.zx,
			 m.xy,m.yy,m.zy,
			 m.xz,m.yz,m.zz);
}

#define JACOBI_NO_TRIANGLE_FUNCTION
inline void Jacobi3x3(const PS::F64mat &A, PS::F64 lambda[3], PS::F64 U[9]){
  const PS::F64 eps = 1e-8;
  PS::F64mat3 A0 = A;
  PS::F64mat3 A1 = A;
  PS::F64mat3asym U0(1.0,0.0,0.0,
		     0.0,1.0,0.0,
		     0.0,0.0,1.0);
  PS::F64mat3asym U1 = U0;
  PS::F64 theta,cs,sn,cs2,sn2,cssn,a,b,c;

  PS::F64 max = 1.0;
  PS::S32 count = 0;
  //while(max > eps){
#pragma unroll 4
  for(int i=0;i<4;i++){
#ifdef JACOBI_NO_TRIANGLE_FUNCTION
      a = 0.5*(A0.xx - A0.yy);
      b = A0.xy;
      c = fabs(b)<eps ? 1.0 : fabs(a) / sqrt(a*a + b*b);
      sn = sqrt(0.5 - 0.5*c) * (a*b<0.0 ? -1.0 : 1.0);
      cs = sqrt(0.5 + 0.5*c);
#else
      if(fabs(A0.xx-A0.yy)<eps) theta = 0.25*M_PI;
      else                      theta = 0.5*atan(2.0*A0.xy / (A0.xx - A0.yy));
      cs = cos(theta); sn = sin(theta);
#endif
      cs2 = cs*cs; sn2 = sn*sn; cssn = cs*sn;
      A1.xx = A0.xx*cs2 + 2.0*A0.xy*cssn + A0.yy*sn2;
      A1.yy = A0.yy*cs2 - 2.0*A0.xy*cssn + A0.xx*sn2;
      //A1.zz = A0.zz;
      A1.xy = A1.xy = A0.xy*(cs2-sn2) + (A0.yy-A0.xx)*cssn;
      A1.xz = A0.xz*cs + A0.yz*sn;
      A1.yz = A0.yz*cs - A0.xz*sn;
      A0 = A1;

      U1.xx = U0.xx*cs + U0.xy*sn;
      U1.xy = U0.xy*cs - U0.xx*sn;
      //U1.xz = U0.xz;
      U1.yx = U0.yx*cs + U0.yy*sn;
      U1.yy = U0.yy*cs - U0.yx*sn;
      //U1.yz = U0.yz;
      U1.zx = U0.zx*cs + U0.zy*sn;
      U1.zy = U0.zy*cs - U0.zx*sn;
      //U1.zz = U0.zz;
      U0 = U1;

#ifdef JACOBI_NO_TRIANGLE_FUNCTION
      a = 0.5*(A0.yy - A0.zz);
      b = A0.yz;
      c = fabs(b)<eps ? 1.0 : fabs(a) / sqrt(a*a + b*b);
      sn = sqrt(0.5 - 0.5*c) * (a*b<0.0 ? -1.0 : 1.0);
      cs = sqrt(0.5 + 0.5*c);
#else
      if(fabs(A0.yy-A0.zz)<eps) theta = 0.25*M_PI;
      else                      theta = 0.5*atan(2.0*A0.yz / (A0.yy - A0.zz));
      cs = cos(theta); sn = sin(theta);
#endif
      cs2 = cs*cs; sn2 = sn*sn; cssn = cs*sn;
      //A1.xx = A0.xx
      A1.yy = A0.yy*cs2 + 2.0*A0.yz*cssn + A0.zz*sn2;
      A1.zz = A0.zz*cs2 - 2.0*A0.yz*cssn + A0.yy*sn2;
      A1.xy = A0.xy*cs + A0.xz*sn;
      A1.xz = A0.xz*cs - A0.xy*sn;
      A1.yz = A0.yz*(cs2 - sn2) + (A0.zz - A0.yy)*cssn;
      A0 = A1;

      //U1.xx = U0.xx;
      U1.xy = U0.xy*cs + U0.xz*sn;
      U1.xz = U0.xz*cs - U0.xy*sn;
      //U1.yx = U0.yx;
      U1.yy = U0.yy*cs + U0.yz*sn;
      U1.yz = U0.yz*cs - U0.yy*sn;
      //U1.zx = U0.zx;
      U1.zy = U0.zy*cs + U0.zz*sn;
      U1.zz = U0.zz*cs - U0.zy*sn;
      U0 = U1;

#ifdef JACOBI_NO_TRIANGLE_FUNCTION
      a = 0.5*(A0.xx - A0.zz);
      b = A0.xz;
      c = fabs(b)<eps ? 1.0 : fabs(a) / sqrt(a*a + b*b);
      sn = sqrt(0.5 - 0.5*c) * (a*b<0.0 ? -1.0 : 1.0);
      cs = sqrt(0.5 + 0.5*c);
#else
      if(fabs(A0.xx-A0.zz)<eps) theta = 0.25*M_PI;
      else                      theta = 0.5*atan(2.0*A0.xz / (A0.xx - A0.zz));
      cs = cos(theta); sn = sin(theta);
#endif
      cs2 = cs*cs; sn2 = sn*sn; cssn = cs*sn;
      A1.xx = A0.xx*cs2 + 2.0*A0.xz*cssn + A0.zz*sn2;
      //A1.yy = A0.yy;
      A1.zz = A0.zz*cs2 - 2.0*A0.xz*cssn + A0.xx*sn2;
      A1.xy = A0.xy*cs + A0.yz*sn;
      A1.xz = A0.xz*(cs2-sn2) + (A0.zz - A0.xx)*cssn;
      A1.yz = A0.yz*cs - A0.xy*sn;
      A0 = A1;

      U1.xx = U0.xx*cs + U0.xz*sn;
      //U1.xy =  U0.xy
      U1.xz = U0.xz*cs - U0.xx*sn;
      U1.yx = U0.yx*cs + U0.yz*sn;
      //U1.yy =  U0.yy*cs + U0.yz*sn;
      U1.yz = U0.yz*cs - U0.yx*sn;
      U1.zx = U0.zx*cs + U0.zz*sn;
      //U1.zy =  U0.zy*cs + U0.zz*sn;
      U1.zz = U0.zz*cs - U0.zx*sn;
      U0 = U1;

    max = std::max(fabs(A0.xy),std::max(fabs(A0.yz),fabs(A0.xz)));
  }
  // normalize eigenvector
#if 0
  printf("A0:(count = %d)\n",count);
  printf("%e %e %e\n",A0.xx,A0.xy,A0.xz);
  printf("%e %e %e\n",A0.yx,A0.yy,A0.yz);
  printf("%e %e %e\n",A0.zx,A0.zy,A0.zz);
#endif
  PS::F64 evl[3] = {A0.xx,A0.yy,A0.zz};
  PS::F64vec evc[3];
  evc[0] = PS::F64vec(U0.xx,U0.yx,U0.zx);
  evc[1] = PS::F64vec(U0.xy,U0.yy,U0.zy);
  evc[2] = PS::F64vec(U0.xz,U0.yz,U0.zz);
#if 0
  printf("eingenvalue:\n");
  printf("%e %e %e\n",evl[0],evl[1],evl[2]);
  printf("eingenvectors:\n");
  printf("%e %e %e\n",evc[0][0],evc[0][1],evc[0][2]);
  printf("%e %e %e\n",evc[1][0],evc[1][1],evc[1][2]);
  printf("%e %e %e\n",evc[2][0],evc[2][1],evc[2][2]);
  printf("check: %e %e %e\n",evc[0]*evc[1],evc[0]*evc[2],evc[1]*evc[2]);
#endif
  evc[0] = evc[0] / sqrt(evc[0]*evc[0]);
  evc[1] = evc[1] / sqrt(evc[1]*evc[1]);
  evc[2] = evc[2] / sqrt(evc[2]*evc[2]);
  // rearrange eigenvalues in ascending order
  if(evl[0]>evl[1]){
    std::swap(evl[0],evl[1]);
    std::swap(evc[0],evc[1]);
  }
  if(evl[0]>evl[2]){
    std::swap(evl[0],evl[2]);
    std::swap(evc[0],evc[2]);
  }
  if(evl[1]>evl[2]){
    std::swap(evl[1],evl[2]);
    std::swap(evc[1],evc[2]);
  }
  for(int i=0;i<3;i++){
    lambda[i] = evl[i];
    for(int j=0;j<3;j++)
      U[3*i+j] = evc[i][j];
  }
}

PS::F64mat CalcQuadrupoleMoment(const PS::F64mat &quad){
  PS::F64mat qmom;

  PS::F64 tmp = quad.xx + quad.yy + quad.zz;
  qmom.xx = 1.5*quad.xx - 0.5*tmp;
  qmom.yy = 1.5*quad.yy - 0.5*tmp;
  qmom.zz = 1.5*quad.zz - 0.5*tmp;
  qmom.xy = 1.5*quad.xy;
  qmom.xz = 1.5*quad.xz;
  qmom.yz = 1.5*quad.yz;

  return qmom;
}

PS::F64mat CalcQuadrupole(const PseudoParticle pp[4]){
  PS::F64mat quad = 0.0;
  for(int i=0;i<4;i++){
    quad.xx += pp[i].mass * pp[i].pos.x * pp[i].pos.x;
    quad.yy += pp[i].mass * pp[i].pos.y * pp[i].pos.y;
    quad.zz += pp[i].mass * pp[i].pos.z * pp[i].pos.z;
    quad.xy += pp[i].mass * pp[i].pos.x * pp[i].pos.y;
    quad.xz += pp[i].mass * pp[i].pos.x * pp[i].pos.z;
    quad.yz += pp[i].mass * pp[i].pos.y * pp[i].pos.z;
  }
  return quad;
}

class MomentQuadrupole{
 private:
void GeneratePseudoParticles(){
  MomentQuadrupole mom(*this);
  //generate a pseudoparticle which has same dipole of mom
  if(mom.di*mom.di == 0.0){
    pp[3].pos.x = 0.0;
    pp[3].pos.y = 0.0;
    pp[3].pos.z = 0.0;
    pp[3].mass = 0.0;
  }else{
    PS::F64 m;
    if(mom.mass > 0.0) m = -5.0*sqrtf(mom.di*mom.di);
    else               m =  5.0*sqrtf(mom.di*mom.di);
    pp[3].pos.x = mom.di.x / m;
    pp[3].pos.y = mom.di.y / m;
    pp[3].pos.z = mom.di.z / m;
    pp[3].mass  = m;
    //remove multipoles of pp[3] from moment  
    const PS::F64vec r = pp[3].pos;
    mom.mass    -= pp[3].mass;
    mom.di      -= pp[3].mass * r;
    mom.quad.xx -= pp[3].mass * r.x * r.x;
    mom.quad.yy -= pp[3].mass * r.y * r.y;
    mom.quad.zz -= pp[3].mass * r.z * r.z;
    mom.quad.xy -= pp[3].mass * r.x * r.y;
    mom.quad.xz -= pp[3].mass * r.x * r.z;
    mom.quad.yz -= pp[3].mass * r.y * r.z;
    assert(fabs(mom.di.x)<1e-5 && fabs(mom.di.y)<1e-5 && fabs(mom.di.z)<1e-5);
  }
  //generate other 3 pseudopoarticles for mono and quad
  PS::F64mat q = CalcQuadrupoleMoment(mom.quad);
  if(mom.mass < 0){
    q = PS::F64mat(-q.xx,-q.yy,-q.zz,-q.xy,-q.xz,-q.yz);
  }
  double w[3],A[9];
  Jacobi3x3(q,w,A); // eigenvalues w is in 昇順
  const PS::F64 mass = fabs(mom.mass);
  assert(w[2]>=w[1] && w[1]>=w[0]);
  //assert(fabs(w[2]+w[1]+w[0]) < 1e-5);
  const PS::F64 alpha = sqrt((2.0*w[2] + w[1]) / mass);
  const PS::F64 beta  =
    (w[2]+2.0*w[1] <= 0.0) ? 0.0 : sqrt((w[2] + 2.0*w[1]) / (3.0*mass));

  // center of PPs must be (0,0,0)
  pp[0].pos = PS::F64vec(   0.0, 2.0*beta, 0.0);
  pp[1].pos = PS::F64vec( alpha,    -beta, 0.0);
  pp[2].pos = PS::F64vec(-alpha,    -beta, 0.0);
  pp[0].mass = pp[1].mass = pp[2].mass = mom.mass/3.0;
  // convert to original coordinate
  const PS::F64vec e2(A[0],A[1],A[2]);
  const PS::F64vec e1(A[3],A[4],A[5]);
  const PS::F64vec e0(A[6],A[7],A[8]);
  for(int k=0;k<3;k++){
    const PS::F64vec tmp = pp[k].pos;
    pp[k].pos = tmp.x*e0 + tmp.y*e1 + tmp.z*e2;
  }
#if 0
  if((pp[0].pos*pp[0].pos > 100.) || (pp[1].pos*pp[1].pos > 100.) || (pp[2].pos*pp[2].pos > 100.) || (pp[3].pos*pp[3].pos > 100.)   )
  const PS::F64vec tmp = PS::F64vec(mom.quad.xx,mom.quad.yy,mom.quad.zz);
  printf("%16.14lf %e %e\n",sqrt(pp[0].pos*pp[0].pos),pp[0].mass,sqrt(tmp*tmp));
  printf("%16.14lf %e %e\n",sqrt(pp[1].pos*pp[1].pos),pp[1].mass,sqrt(tmp*tmp));
  printf("%16.14lf %e %e\n",sqrt(pp[2].pos*pp[2].pos),pp[2].mass,sqrt(tmp*tmp));
  printf("%16.14lf %e %e\n",sqrt(pp[3].pos*pp[3].pos),pp[3].mass,sqrt(tmp*tmp));
#endif
#if 0
  static int count = 0;
  printf("%d %d %e %e %e\n",count++,0,pp[0].pos.x,pp[0].pos.y,pp[0].pos.z);
  printf("%d %d %e %e %e\n",count++,1,pp[1].pos.x,pp[1].pos.y,pp[1].pos.z);
  printf("%d %d %e %e %e\n",count++,2,pp[2].pos.x,pp[2].pos.y,pp[2].pos.z);
  printf("%d %d %e %e %e\n",count++,3,pp[3].pos.x,pp[3].pos.y,pp[3].pos.z);
#endif
}
  
 public:
  PS::S32    nptcl;
  PS::F64vec pos;
  PS::F64    mass;
  PS::F64vec di;
  PS::F64mat quad;
  PseudoParticle pp[N_PSEUDOPARTICLE];
  PS::F64ort vertex_out_;
  
  MomentQuadrupole(){ init(); }
  MomentQuadrupole(const MomentQuadrupole& m){
    nptcl = m.nptcl;
    pos   = m.pos;
    mass  = m.mass;
    di    = m.di;
    quad  = m.quad;
    for(int i=0;i<N_PSEUDOPARTICLE;i++) pp[i] = m.pp[i];
    vertex_out_ = m.vertex_out_;
  }

  PS::F64vec getPos() const {return pos;}
  PS::F64ort getVertexOut() const {return vertex_out_;}

  void init(){
    nptcl = 0;
    pos   = 0.0;
    mass  = 0.0;
    di    = 0.0;
    quad  = 0.0;
    for(int i=0;i<N_PSEUDOPARTICLE;i++) pp[i].init();
  }
  template<class Tepj>
  void accumulateAtLeaf(const Tepj& epj){
    nptcl++;
    pos  += epj.getPos();
    mass += epj.getCharge();
    vertex_out_.merge(epj.getPos(),epj.getRSearch());
    //vertex_out_.merge(epj.getPos());
  }
  void accumulate(const MomentQuadrupole& m){
    nptcl += m.nptcl;
    pos   += m.nptcl * m.pos;
    mass  += m.mass;
    vertex_out_.merge(m.vertex_out_);
  }
  void set(){
    pos = pos / (PS::F64)nptcl;
  }
  template <class Tepj>
  void accumulateAtLeaf2(const Tepj& epj){
    const PS::F64vec r = epj.getPos() - pos;
    di += epj.getCharge() * r;
    quad.xx += epj.getCharge() * r.x * r.x;
    quad.yy += epj.getCharge() * r.y * r.y;
    quad.zz += epj.getCharge() * r.z * r.z;
    quad.xy += epj.getCharge() * r.x * r.y;
    quad.xz += epj.getCharge() * r.x * r.z;
    quad.yz += epj.getCharge() * r.y * r.z;
  }
  void accumulate2(const MomentQuadrupole& m){
    const PS::F64vec r = m.pos - pos;
    di += m.di - m.mass * r;
    quad.xx += m.quad.xx - 2.0 * m.di.x * r.x + m.mass * r.x * r.x;
    quad.yy += m.quad.yy - 2.0 * m.di.y * r.y + m.mass * r.y * r.y;
    quad.zz += m.quad.zz - 2.0 * m.di.z * r.z + m.mass * r.z * r.z;
    quad.xy += m.quad.xy - (m.di.x*r.y + m.di.y*r.x) + m.mass * r.x * r.y;
    quad.xz += m.quad.xz - (m.di.x*r.z + m.di.z*r.x) + m.mass * r.x * r.z;
    quad.yz += m.quad.yz - (m.di.y*r.z + m.di.z*r.y) + m.mass * r.y * r.z;
  }
  void set2(){
    GeneratePseudoParticles();
#ifdef IPS_TREE_PSEUDOPARTICLE_QUADRUPOLE_MOMENT_CHECK
    // check generated pseudo particle reproduce physical quadrupole moment
    {
      const PS::F64mat apprx = CalcQuadrupoleMoment(CalcQuadrupole(pp));
      const PS::F64mat exact = CalcQuadrupoleMoment(quad);
      /*
      printf("pseudo:   %e %e %e %e %e %e\n",
	     apprx.xx, apprx.yy, apprx.zz,
	     apprx.xy, apprx.xz, apprx.yz);
      printf("original: %e %e %e %e %e %e\n",
  	     exact.xx,exact.yy,exact.zz,
	     exact.xy,exact.xz,exact.yz);
      //*/
      PS::F64 tmp = std::max(exact.xx,std::max(exact.yy,std::max(exact.zz,std::max(exact.xy,std::max(exact.xz,exact.yz)))));
      assert(fabs((exact.xx-apprx.xx)/tmp) < 1e-6);
      assert(fabs((exact.yy-apprx.yy)/tmp) < 1e-6);
      assert(fabs((exact.zz-apprx.zz)/tmp) < 1e-6);
      assert(fabs((exact.xy-apprx.xy)/tmp) < 1e-6);
      assert(fabs((exact.xz-apprx.xz)/tmp) < 1e-6);
      assert(fabs((exact.yz-apprx.yz)/tmp) < 1e-6);
    }
#endif
    for(int i=0;i<N_PSEUDOPARTICLE;i++) pp[i].pos += pos;
  }
};

class SPJQuadrupole{
public:
  PseudoParticle pp[N_PSEUDOPARTICLE];
  void copyFromMoment(const MomentQuadrupole & mom){
    for(int i=0;i<N_PSEUDOPARTICLE;i++) pp[i] = mom.pp[i];
  }
  void clear(){
    for(int i=0;i<N_PSEUDOPARTICLE;i++) pp[i].init();
  }

  PS::F32vec getPos() const {
    PS::S32    nptcl = 0;
    PS::F64vec pos   = 0.0;
    for(int i=0;i<N_PSEUDOPARTICLE;i++){
      nptcl++;
      pos += pp[i].pos;
    }
    return pos/(PS::F64)nptcl;
  }
  void setPos(const PS::F64vec & pos_new) {
    const PS::F64vec diff = pos_new - getPos();
    for(int i=0;i<N_PSEUDOPARTICLE;i++) pp[i].pos += diff;
  }
  MomentQuadrupole convertToMoment() const {
    MomentQuadrupole mom;
    for(int i=0;i<N_PSEUDOPARTICLE;i++)
      mom.accumulateAtLeaf(pp[i]);
    mom.set();
    for(int i=0;i<N_PSEUDOPARTICLE;i++)
      mom.accumulateAtLeaf2(pp[i]);
    mom.vertex_out_.init();
    return MomentQuadrupole(mom);
  }

  void DebugPrint() const {
    for(int k=0;k<N_PSEUDOPARTICLE;k++){
      printf("pp[%d]: %lf %lf %lf %lf\n",k,pp[k].mass,pp[k].pos.x,pp[k].pos.y,pp[k].pos.z);
    }
  }
};
#endif
