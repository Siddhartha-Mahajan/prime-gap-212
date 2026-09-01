/* Independent directed-rounding verifier for a fixed, dyadic sieve witness.
   No eigenvalue routines, FFTs, transcendental functions, or external libraries.
   Compile: g++ -O3 -std=c++17 -fopenmp -frounding-math -ffp-contract=off
                 -fno-fast-math certify.cpp -o certify
   Run: OMP_NUM_THREADS=4 ./certify witness.hex i j
   Output endpoints are exact hexadecimal representations of binary64 numbers.
*/
#pragma STDC FENV_ACCESS ON
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits==53, "IEEE binary64 required");
using V=std::vector<double>;
// A polynomial stores coefficients at indices off, off+1, ... .
struct Poly {int off=0; V a;};
struct Data {int k,R,n,D,den,T,L,maxr;std::vector<int>B;std::vector<V>g,h;};
void mode(int m){if(std::fesetround(m))throw std::runtime_error("fesetround failed");}
void selftest(){
 volatile double a=1.,b=0x1p-54;
 mode(FE_DOWNWARD);volatile double lo=a+b;
 mode(FE_UPWARD);volatile double hi=a+b;
 if(lo!=1. || !(hi>1.))throw std::runtime_error("rounding-mode self-test failed");
 volatile double tiny=std::numeric_limits<double>::denorm_min(),half=0.5;
 mode(FE_DOWNWARD);volatile double tl=tiny*half;
 mode(FE_UPWARD);volatile double tu=tiny*half;
 if(tl!=0. || tu!=tiny)throw std::runtime_error("subnormal rounding self-test failed");
 mode(FE_TONEAREST);
}
double readhex(std::istream &in){std::string s;in>>s;if(!in)throw std::runtime_error("short witness");
 char*end=nullptr;double x=std::strtod(s.c_str(),&end);
 if(!end||*end||!std::isfinite(x))throw std::runtime_error("bad hexadecimal input");return x;}
Data load(const char* path){
 mode(FE_TONEAREST);std::ifstream f(path);if(!f)throw std::runtime_error("cannot open witness");
 Data d;f>>d.k>>d.R>>d.n>>d.D>>d.den>>d.T>>d.L>>d.maxr;
 if(!f||d.k<2||d.k>50||d.R<1||d.R>16||d.n<1||d.n>1000000||
    d.D<1||d.D>=d.n||d.den<1||(d.den&(d.den-1))||d.T!=d.n+d.k-1||
    d.L<d.k-1||d.L>d.T||d.maxr<0||d.maxr>d.k-2)
  throw std::runtime_error("invalid witness dimensions");
 d.B.resize(d.maxr+2);for(auto &b:d.B)f>>b;
 if(!f||d.B[0]!=0||d.B[1]<=d.D||d.B[1]>d.n||d.B.back()/d.D!=d.maxr)
  throw std::runtime_error("invalid support bounds");
 for(int r=2;r<int(d.B.size());r++)
  if(d.B[r]<d.B[r-1]||d.B[r]>d.B[r-1]+d.D)
   throw std::runtime_error("support increment condition failed");
 d.g.assign(d.R,V(d.n));d.h.assign(d.R,V(d.n));
 for(int a=0;a<d.n;a++){for(int l=0;l<d.R;l++){d.g[l][a]=readhex(f);if(d.g[l][a]<=0)throw std::runtime_error("nonpositive g");}
  for(int l=0;l<d.R;l++)d.h[l][a]=readhex(f);}
 std::string extra;if(f>>extra)throw std::runtime_error("unexpected trailing data");
 return d;
}
// Direct, truncated convolution of nonnegative coefficients.
// Any summation order preserves the one-sided rounding invariant.
Poly convolve(const Poly&a,const Poly&b,int cap,int rounding){
 Poly c;c.off=a.off+b.off;
 int len=std::min(cap-c.off+1,int(a.a.size()+b.a.size())-1);
 if(len<=0){c.a={};return c;}c.a.resize(len);
 const int an=a.a.size(),bn=b.a.size();
 #pragma omp parallel
 {
  mode(rounding);
  #pragma omp for schedule(static)
  for(int z=0;z<len;z++){
   int first=std::max(0,z-bn+1),last=std::min(an-1,z);double sum=0.;
   #pragma omp simd reduction(+:sum)
   for(int t=first;t<=last;t++)sum+=a.a[t]*b.a[z-t];
   c.a[z]=sum;
  }
 }
 for(double x:c.a)if(x<0||!std::isfinite(x))throw std::runtime_error("invalid positive convolution");
 return c;
}
double at(const Poly&p,int z){z-=p.off;return z>=0&&z<int(p.a.size())?p.a[z]:0.;}
// The exact integer recurrence is safe for n <= 50; check overflow anyway.
double choose(int n,int r){
 uint64_t x=1;
 for(int j=1;j<=r;j++){
  uint64_t factor=uint64_t(n-j+1);
  if(factor && x>std::numeric_limits<uint64_t>::max()/factor)
   throw std::runtime_error("binomial overflow");
  x=(x*factor)/uint64_t(j);
 }
 if(x>9007199254740992ULL)throw std::runtime_error("inexact binomial");
 return double(x);
}
struct Powers {std::vector<Poly>sp,bp;V norm;};
// P.sp[j] and P.bp[r] are the small- and large-coordinate powers.
// P.norm[s] includes the binomial factors in the norm formula.
Powers powers(const Data&d,int i,int j,int rounding){
 mode(rounding);double mesh=1./d.den;
 Poly qs;qs.off=0;qs.a.resize(d.D);
 for(int a=0;a<d.D;a++)qs.a[a]=mesh*(d.g[i][a]*d.g[j][a]);
 Poly qb;qb.off=d.D;qb.a.resize(d.B[1]-d.D);
 for(int a=d.D;a<d.B[1];a++)qb.a[a-d.D]=mesh*(d.g[i][a]*d.g[j][a]);
 Powers P;P.sp.resize(d.k+1);P.bp.resize(d.maxr+1);P.norm.assign(d.n,0.);
 Poly cur{0,V{1.}};
 for(int r=1;r<=d.k;r++){
  cur=convolve(cur,qs,d.n-1,rounding);
  if(r>=d.k-d.maxr-1)P.sp[r]=cur;
 }
 P.bp[0]=Poly{0,V{1.}};
 for(int r=1;r<=d.maxr;r++)P.bp[r]=convolve(P.bp[r-1],qb,d.B[r]-r,rounding);
 for(int r=0;r<=d.maxr;r++){
  Poly t=convolve(P.sp[d.k-r],P.bp[r],d.n-1,rounding);
  mode(rounding);double mult=choose(d.k,r);
  for(int a=0;a<int(t.a.size());a++)P.norm[t.off+a]+=mult*t.a[a];
 }
 return P;
}
// The radial product can have either sign: a negative product reverses
// the choice of a positive coefficient endpoint.
std::pair<double,double> norm_integral(const Data&d,int i,int j,const Powers&lo,const Powers&hi){
 double out[2];
 for(int phase=0;phase<2;phase++){
  mode(phase==0?FE_DOWNWARD:FE_UPWARD);double sum=0.;
  for(int s=0;s<d.n;s++){
   bool neg=(d.h[i][s]<0)!=(d.h[j][s]<0);
   const V &v=((phase==0)^neg)?lo.norm:hi.norm;
   double prod=d.h[i][s]*d.h[j][s];sum+=v[s]*prod;
  }out[phase]=sum;
 }if(!std::isfinite(out[0])||!std::isfinite(out[1])||out[0]>out[1])throw std::runtime_error("invalid I enclosure");
 return {out[0],out[1]};
}
// Group common-coordinate coefficients by their final fiber cutoff u.
void weights(V&W,const Data&d,const Powers&P,int s,int last,int rounding){
 mode(rounding);std::fill(W.begin(),W.begin()+last+1,0.);
 for(int r=0;r<=d.maxr;r++){
  int lastb=std::min(s,int(P.bp[r].off+P.bp[r].a.size())-1);
  double mult=choose(d.k-1,r);
  for(int b=P.bp[r].off;b<=lastb;b++){
   int u=std::min(d.n-1-s,std::max(d.D-1,d.B[r+1]-r-1-b));
   if(u<0||u>last)throw std::runtime_error("weight subscript");
   W[u]+=mult*(at(P.bp[r],b)*at(P.sp[d.k-1-r],s-b));
  }
 }
}
// Enclose the signed fiber prefixes H_i(s,u) and H_j(s,u).
void prefix(V&pi,V&pj,const Data&d,int i,int j,int s,int last,int rounding){
 mode(rounding);double si=0.,sj=0.;
 for(int u=0;u<=last;u++){si+=d.g[i][u]*d.h[i][s+u];sj+=d.g[j][u]*d.h[j][s+u];pi[u]=si;pj[u]=sj;}
}
// The four-endpoint product encloses all signed cross terms in J_box.
// Each thread has its own rounding mode and one-sided accumulators.
std::pair<double,double> marginal_integral(const Data&d,int i,int j,const Powers&lo,const Powers&hi){
 const int maxs=std::min(d.n-1,d.L-(d.k-1));
 V lower(omp_get_max_threads(),0.),upper(omp_get_max_threads(),0.);
 #pragma omp parallel
 {
  int th=omp_get_thread_num();int nu=d.B[1];V wl(nu),wh(nu),il(nu),ih(nu),jl(nu),jh(nu);
  double sl=0.,sh=0.;
  #pragma omp for schedule(static)
  for(int s=0;s<=maxs;s++){
   int last=std::min(d.n-1-s,nu-1);
   weights(wl,d,lo,s,last,FE_DOWNWARD);prefix(il,jl,d,i,j,s,last,FE_DOWNWARD);
   weights(wh,d,hi,s,last,FE_UPWARD);prefix(ih,jh,d,i,j,s,last,FE_UPWARD);
   mode(FE_DOWNWARD);double row=0.;
   for(int u=0;u<=last;u++){
    double z=std::min(std::min(il[u]*jl[u],il[u]*jh[u]),std::min(ih[u]*jl[u],ih[u]*jh[u]));
    row+=(z<0?wh[u]:wl[u])*z;
   }sl+=row;
   mode(FE_UPWARD);row=0.;
   for(int u=0;u<=last;u++){
    double z=std::max(std::max(il[u]*jl[u],il[u]*jh[u]),std::max(ih[u]*jl[u],ih[u]*jh[u]));
    row+=(z<0?wl[u]:wh[u])*z;
   }sh+=row;
  }
  lower[th]=sl;upper[th]=sh;
 }
 mode(FE_DOWNWARD);double a=0.;for(double v:lower)a+=v;a=(a/d.den)/d.den;
 mode(FE_UPWARD);double b=0.;for(double v:upper)b+=v;b=(b/d.den)/d.den;
 if(!std::isfinite(a)||!std::isfinite(b)||a>b)throw std::runtime_error("invalid J enclosure");
 return {a,b};
}
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("usage: certify witness.hex i j");selftest();Data d=load(argv[1]);
 int i=std::atoi(argv[2]),j=std::atoi(argv[3]);if(i<0||j<i||j>=d.R)throw std::runtime_error("bad pair");
 double start=omp_get_wtime();
 Powers lo=powers(d,i,j,FE_DOWNWARD);Powers hi=powers(d,i,j,FE_UPWARD);
 auto I=norm_integral(d,i,j,lo,hi);auto J=marginal_integral(d,i,j,lo,hi);
 mode(FE_TONEAREST);
 std::cout<<"{\"i\":"<<i<<",\"j\":"<<j<<",\"Ilo\":\""<<std::hexfloat<<I.first<<"\",\"Ihi\":\""<<I.second
 <<"\",\"Jlo\":\""<<J.first<<"\",\"Jhi\":\""<<J.second<<"\",\"seconds\":"<<std::defaultfloat<<std::setprecision(10)<<(omp_get_wtime()-start)<<"}\n";
 }catch(std::exception&e){std::cerr<<e.what()<<"\n";return 1;}return 0;}
