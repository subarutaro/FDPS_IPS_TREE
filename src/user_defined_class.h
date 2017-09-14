#ifndef H_USER_DEFINED_CLASS
#define H_USER_DEFINED_CLASS

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
    return this->search_radius;
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
  PS::F64 getRSearch() const {
    return search_radius;
  }
  PS::F64vec getPos() const { return pos;}
  PS::F64    getCharge() const {return charge;}
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

class MomentUnit{
 public:
  PS::S32    nptcl;
  PS::F64vec pos;
  PS::F64    mass;
  PS::F64vec dipole;
  PS::F64mat quadrupole;

  MomentUnit(){
    init();
  }
  MomentUnit(const MomentUnit& m){
    nptcl      = m.nptcl;
    pos        = m.pos;
    mass       = m.mass;
    dipole     = m.dipole;
    quadrupole = m.quadrupole;
  }

  void init(){
    nptcl      = 0;
    pos        = 0.0;
    mass       = 0.0;
    dipole     = 0.0;
    quadrupole = 0.0;
  }
  void set(){
    if(mass == 0.0) pos = 0.0;
    else            pos = pos / mass; // gravity center as a center of multipole expansion
  }
  template<class Tepj>
  void accumulateAtLeaf(const Tepj& epj){
    nptcl++;
    pos  += epj.getCharge()*epj.getPos();
    mass += epj.getCharge();
  }
  template <class Tepj>
  void accumulateAtLeaf2(const Tepj& epj){
    const PS::F64vec r = epj.getPos() - pos;
    dipole += epj.getCharge() * r;
    quadrupole.xx += epj.getCharge() * r.x * r.x;
    quadrupole.yy += epj.getCharge() * r.y * r.y;
    quadrupole.zz += epj.getCharge() * r.z * r.z;
    quadrupole.xy += epj.getCharge() * r.x * r.y;
    quadrupole.xz += epj.getCharge() * r.x * r.z;
    quadrupole.yz += epj.getCharge() * r.y * r.z;
  }
  void accumulate(const MomentUnit& m){
    nptcl += m.nptcl;
    pos   += m.mass * m.pos;
    mass  += m.mass;
  }
  void accumulate2(const MomentUnit& m){
    const PS::F64vec r = m.pos - pos;
    dipole += m.dipole - m.mass * r;
    quadrupole.xx += m.quadrupole.xx - 2.0 * m.dipole.x * r.x + m.mass * r.x * r.x;
    quadrupole.yy += m.quadrupole.yy - 2.0 * m.dipole.y * r.y + m.mass * r.y * r.y;
    quadrupole.zz += m.quadrupole.zz - 2.0 * m.dipole.z * r.z + m.mass * r.z * r.z;
    quadrupole.xy += m.quadrupole.xy - (m.dipole.x*r.y + m.dipole.y*r.x) + m.mass * r.x * r.y;
    quadrupole.xz += m.quadrupole.xz - (m.dipole.x*r.z + m.dipole.z*r.x) + m.mass * r.x * r.z;
    quadrupole.yz += m.quadrupole.yz - (m.dipole.y*r.z + m.dipole.z*r.y) + m.mass * r.y * r.z;
  }
};

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

  PS::F64 getCharge() const {
    return mass;
  }
  PS::F64vec getPos() const {
    return pos;
  }
};

extern "C"
int dsyev_
(const char *jobz,const char *uplo,
 int *n, double *a, int *lda,
 double *w, double *work, int *lwork, int *info);

void calc_eigen_value(const PS::F64mat &quadrupole,double w[3],double A[9]){
  A[0]=quadrupole.xx; A[1]=quadrupole.xy; A[2]=quadrupole.xz;
  A[3]=quadrupole.xy; A[4]=quadrupole.yy; A[5]=quadrupole.yz;
  A[6]=quadrupole.xz; A[7]=quadrupole.yz; A[8]=quadrupole.zz;
  int n = 3;
  double work;
  int lwork=-1,info;
  dsyev_("V","U",&n,A,&n,w,&work,&lwork,&info);
  assert(info == 0);
  lwork = (int) work;
  double work2[lwork>1 ? lwork:1];
  dsyev_("V","L",&n,A,&n,w,work2,&lwork,&info);
  assert(info == 0);
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

PS::F64mat CalcQuadrupole(const PseudoParticle pp[3]){
  PS::F64mat quad = 0.0;
  for(int i=0;i<3;i++){
    quad.xx += pp[i].mass * pp[i].pos.x * pp[i].pos.x;
    quad.yy += pp[i].mass * pp[i].pos.y * pp[i].pos.y;
    quad.zz += pp[i].mass * pp[i].pos.z * pp[i].pos.z;
    quad.xy += pp[i].mass * pp[i].pos.x * pp[i].pos.y;
    quad.xz += pp[i].mass * pp[i].pos.x * pp[i].pos.z;
    quad.yz += pp[i].mass * pp[i].pos.y * pp[i].pos.z;
  }
  return quad;
}

void GeneratePseudoParticles(PseudoParticle pp[3],
			     const MomentUnit& mom){
  if(mom.mass == 0.0){
    pp[0].mass = pp[1].mass = pp[2].mass = 0.0;
    pp[0].pos = pp[1].pos = pp[2].pos = 0.0;
    return;
  }
#if 0
  if(mom.nptcl < 3){
    // place a particle with mass M at center of mass (= 0.0)
    pp[0].pos  = 0.0;
    pp[0].mass = mom.mass;
    // other particle does not exist (masses must be 0.0)
    pp[1].pos  = pp[2].pos  = 0.0;
    pp[1].mass = pp[2].mass = 0.0;
    return;
  }
#endif
  PS::F64mat quad = CalcQuadrupoleMoment(mom.quadrupole);
  if(mom.mass < 0.0){
    quad = PS::F64mat(-quad.xx,-quad.yy,-quad.zz,-quad.xy,-quad.xz,-quad.yz);
  }
  double w[3],A[9];
  calc_eigen_value(quad,w,A); // eigenvalues w is in 昇順

  const PS::F64 mass = fabs(mom.mass);
  assert(w[2]>=w[1] && w[1]>=w[0]);
  assert(fabs(w[2]+w[1]+w[0]) < 1e-5);
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
}

class MomentQuadrupoleSixPseudoParticleCutoff{
public:
  PS::F64    mass;
  PS::F64vec pos;
  MomentUnit positive;
  MomentUnit negative;
  PS::F64ort vertex_out_;

  MomentQuadrupoleSixPseudoParticleCutoff(){
    mass = 0.0;
    pos  = 0.0;
    positive.init();
    negative.init();
    vertex_out_.init();
  }
  MomentQuadrupoleSixPseudoParticleCutoff(const PS::F64 & _m,
					  const PS::F64vec & _p,
					  const MomentUnit & _pos,
					  const MomentUnit & _neg,
					  const PS::F64ort & v_out){
    mass = _m;
    pos  = _p;
    positive = _pos;
    negative = _neg;
    vertex_out_ = v_out;
  }
  void init(){
    mass = 0.0;
    pos  = 0.0;
    positive.init();
    negative.init();
    vertex_out_.init();
  }

  PS::F64vec getPos() const {
    return pos;
  }
  PS::F64ort getVertexOut() const { return vertex_out_; }

  template<class Tepj>
  void accumulateAtLeaf(Tepj & epj){
    if(epj.getCharge() > 0.0){
      mass += epj.getCharge();
      pos  += epj.getCharge()*epj.getPos();
      positive.accumulateAtLeaf(epj);
    }else if(epj.getCharge() < 0.0){
      mass -= epj.getCharge();
      pos  -= epj.getCharge()*epj.getPos();
      negative.accumulateAtLeaf(epj);
    }else{
      assert(false && "error: no charge particle exist!");
    }
    vertex_out_.merge(epj.getPos(), epj.getRSearch());
  }
  template<class Tepj>
  void accumulateAtLeaf2(Tepj & epj){
    if(epj.getCharge() > 0.0){
      positive.accumulateAtLeaf2(epj);
    }else if(epj.getCharge() < 0.0){
      negative.accumulateAtLeaf2(epj);
    }
  }
  void set(){
    pos = pos / mass;
    positive.set();
    negative.set();
  }
  void accumulate(const MomentQuadrupoleSixPseudoParticleCutoff & mom){
    const PS::F64 m = fabs(mom.mass);
    mass += m;
    pos  += m * mom.pos;
    positive.accumulate(mom.positive);
    negative.accumulate(mom.negative);
    vertex_out_.merge(mom.vertex_out_);
  }
  void accumulate2(const MomentQuadrupoleSixPseudoParticleCutoff & mom){
    positive.accumulate2(mom.positive);
    negative.accumulate2(mom.negative);
  }
};


class SPJQuadrupoleSixPseudoParticleCutoff{
public:
  PseudoParticle pp[6]; // pseudo particle of positive(3) and negative(3) charge
  void copyFromMoment(const MomentQuadrupoleSixPseudoParticleCutoff & mom){
    GeneratePseudoParticles(pp  ,mom.positive);
    GeneratePseudoParticles(pp+3,mom.negative);
#ifdef IPS_TREE_PSEUDOPARTICLE_QUADRUPOLE_MOMENT_CHECK
    // check generated pseudo particle reproduce physical quadrupole moment
    {
      const PS::F64mat qm_ppp = CalcQuadrupoleMoment(CalcQuadrupole(pp));
      const PS::F64mat qm_mom = CalcQuadrupoleMoment(mom.positive.quadrupole);
      printf("pseudo: %e %e %e %e %e %e\n",
	     qm_ppp.xx, qm_ppp.yy, qm_ppp.zz,
	     qm_ppp.xy, qm_ppp.xz, qm_ppp.yz);
      printf("original: %e %e %e %e %e %e\n",
  	     qm_mom.xx,qm_mom.yy,qm_mom.zz,
	     qm_mom.xy,qm_mom.xz,qm_mom.yz);
      assert(fabs(qm_mom.xx) < 1e-6 || fabs((qm_mom.xx-qm_ppp.xx)/qm_mom.xx) < 1e-6);
      assert(fabs(qm_mom.yy) < 1e-6 || fabs((qm_mom.yy-qm_ppp.yy)/qm_mom.yy) < 1e-6);
      assert(fabs(qm_mom.zz) < 1e-6 || fabs((qm_mom.zz-qm_ppp.zz)/qm_mom.zz) < 1e-6);
      assert(fabs(qm_mom.xy) < 1e-6 || fabs((qm_mom.xy-qm_ppp.xy)/qm_mom.xy) < 1e-6);
      assert(fabs(qm_mom.xz) < 1e-6 || fabs((qm_mom.xz-qm_ppp.xz)/qm_mom.xz) < 1e-6);
      assert(fabs(qm_mom.yz) < 1e-6 || fabs((qm_mom.yz-qm_ppp.yz)/qm_mom.yz) < 1e-6);
    }
    {
      const PS::F64mat qm_ppn = CalcQuadrupoleMoment(CalcQuadrupole(pp+3));
      const PS::F64mat qm_mom = CalcQuadrupoleMoment(mom.negative.quadrupole);
      fprintf(stderr,
	      "pseudo: %e %e %e %e %e %e\n",
	      qm_ppn.xx, qm_ppn.yy, qm_ppn.zz,
	      qm_ppn.xy, qm_ppn.xz, qm_ppn.yz);
      fprintf(stderr,
	      "original: %e %e %e %e %e %e\n",
	      qm_mom.xx,qm_mom.yy,qm_mom.zz,
	      qm_mom.xy,qm_mom.xz,qm_mom.yz);
      assert(fabs(qm_mom.xx) < 1e-6 || fabs((qm_mom.xx-qm_ppn.xx)/qm_mom.xx) < 1e-6);
      assert(fabs(qm_mom.yy) < 1e-6 || fabs((qm_mom.yy-qm_ppn.yy)/qm_mom.yy) < 1e-6);
      assert(fabs(qm_mom.zz) < 1e-6 || fabs((qm_mom.zz-qm_ppn.zz)/qm_mom.zz) < 1e-6);
      assert(fabs(qm_mom.xy) < 1e-6 || fabs((qm_mom.xy-qm_ppn.xy)/qm_mom.xy) < 1e-6);
      assert(fabs(qm_mom.xz) < 1e-6 || fabs((qm_mom.xz-qm_ppn.xz)/qm_mom.xz) < 1e-6);
      assert(fabs(qm_mom.yz) < 1e-6 || fabs((qm_mom.yz-qm_ppn.yz)/qm_mom.yz) < 1e-6);
    }
#endif
    for(int i=0;i<3;i++){
      pp[i  ].pos += mom.positive.pos;
      pp[i+3].pos += mom.negative.pos;
    }
  }
  void clear(){
    pp[0].init();
    pp[1].init();
    pp[2].init();
    pp[3].init();
    pp[4].init();
    pp[5].init();
  }

  PS::F32vec getPos() const {
    PS::F64    mass = 0.0;
    PS::F64vec pos  = 0.0;
    for(int i=0;i<6;i++){
      mass += fabs(pp[i].mass);
      pos  += fabs(pp[i].mass)*pp[i].pos;
    }    
    return pos/mass;
  }
  void setPos(const PS::F64vec & pos_new) {
    PS::F64    mass = 0.0;
    PS::F64vec pos  = 0.0;
    for(int i=0;i<6;i++){
      mass += fabs(pp[i].mass);
      pos  += fabs(pp[i].mass)*pp[i].pos;
    }
    pos = pos/mass;

    const PS::F64vec diff = pos_new - pos;
    for(int i=0;i<6;i++){
      if(pp[i].mass != 0.0) pp[i].pos += diff;
    }
  }
  MomentQuadrupoleSixPseudoParticleCutoff convertToMoment() const {
    PS::F64    mass = 0.0;
    PS::F64vec pos  = 0.0;
    for(int i=0;i<6;i++){
      mass += fabs(pp[i].mass);
      pos  += fabs(pp[i].mass)*pp[i].pos;
    }
    pos = pos / mass;
    MomentUnit positive,negative;
    for(int i=0;i<3;i++){
      positive.accumulateAtLeaf(pp[i]);
      negative.accumulateAtLeaf(pp[i+3]);
    }
    positive.set();
    negative.set();

    PS::F64ort vtx; // dummy vertex
    vtx.init();
    /*
    for(int i=0;i<6;i++){
      vtx.merge(pp[i].pos);
    }
    //*/
    return MomentQuadrupoleSixPseudoParticleCutoff(mass,pos,positive,negative,vtx);
  }

  void DebugPrint() const {
    for(int k=0;k<6;k++){
      printf("pp[%d]: %lf %lf %lf %lf\n",k,pp[k].mass,pp[k].pos.x,pp[k].pos.y,pp[k].pos.z);
    }
  }
};
#endif
