/* Integer-interval verifier for the finite box witness.
   Every mathematical operation uses exact checked 256-bit integers.
   An interval [lo,hi] represents [lo/2^96,hi/2^96]. Products and
   divisions are rounded outwards by explicit integer floor/ceiling.
   Overflow raises an exception and cannot yield a certificate.
   Compile: g++ -O3 -std=c++17 -fopenmp certify_integer.cpp -o certify_integer
   Run: OMP_NUM_THREADS=4 ./certify_integer witness.fix i j
   Boost.Multiprecision is header-only. Elapsed time is not part of the proof.
*/
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>
namespace mp=boost::multiprecision;
using Z=mp::number<mp::cpp_int_backend<256,256,mp::signed_magnitude,mp::checked,void>>;
using V=std::vector<Z>;
constexpr int PREC=96;
const Z ONE=Z(1)<<PREC;
Z floor_scale(const Z&x){return x>=0?(x>>PREC):-(((-x)+ONE-1)>>PREC);}
Z ceil_scale(const Z&x){return x>=0?((x+ONE-1)>>PREC):-((-x)>>PREC);}
Z floor_div(const Z&x,const Z&d){return x>=0?x/d:-(((-x)+d-1)/d);}
Z ceil_div(const Z&x,const Z&d){return x>=0?(x+d-1)/d:-((-x)/d);}
struct IV{Z lo=0,hi=0;};
IV product(const IV&a,const IV&b){
 std::array<Z,4> p{a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
 return {floor_scale(*std::min_element(p.begin(),p.end())),ceil_scale(*std::max_element(p.begin(),p.end()))};
}
IV pos_product(const IV&a,const IV&b){ // a is nonnegative
 return {floor_scale((b.lo<0?a.hi:a.lo)*b.lo),ceil_scale((b.hi<0?a.lo:a.hi)*b.hi)};
}
struct Poly{int off=0;V lo,hi;};
struct Data{int k,R,n,D,den,T,L,maxr;std::vector<int>B;std::vector<V>g,h;};
Z readint(std::istream&in){Z z;in>>z;if(!in)throw std::runtime_error("invalid integer input");return z;}
Data load(const char*path){
 std::ifstream f(path);if(!f)throw std::runtime_error("cannot open witness");int p;f>>p;
 if(!f||p!=PREC)throw std::runtime_error("precision mismatch");
 Data d;f>>d.k>>d.R>>d.n>>d.D>>d.den>>d.T>>d.L>>d.maxr;
 if(!f||d.k<2||d.k>50||d.R<1||d.R>16||d.n<1||d.n>1000000||
    d.D<1||d.D>=d.n||d.den<1||(d.den&(d.den-1))||d.T!=d.n+d.k-1||
    d.L<d.k-1||d.L>d.T||d.maxr<0||d.maxr>d.k-2)
  throw std::runtime_error("invalid witness dimensions");
 d.B.resize(d.maxr+2);for(auto&b:d.B)f>>b;
 if(!f||d.B[0]!=0||d.B[1]<=d.D||d.B[1]>d.n||d.B.back()/d.D!=d.maxr)
  throw std::runtime_error("invalid support bounds");
 for(int r=2;r<int(d.B.size());r++)
  if(d.B[r]<d.B[r-1]||d.B[r]>d.B[r-1]+d.D)
   throw std::runtime_error("support increment condition failed");
 d.g.assign(d.R,V(d.n));d.h.assign(d.R,V(d.n));
 for(int s=0;s<d.n;s++){
  for(int l=0;l<d.R;l++){d.g[l][s]=readint(f);if(d.g[l][s]<=0)throw std::runtime_error("nonpositive g");}
  for(int l=0;l<d.R;l++)d.h[l][s]=readint(f);
 }
 std::string extra;if(f>>extra)throw std::runtime_error("trailing witness data");return d;
}
// All polynomial coefficients are nonnegative. Accumulation is exact;
// one outward division by 2^96 is performed per output coefficient.
Poly convolve(const Poly&a,const Poly&b,int cap){
 Poly c;c.off=a.off+b.off;if(a.lo.empty()||b.lo.empty())return c;
 int len=std::min(cap-c.off+1,int(a.lo.size()+b.lo.size())-1);
 if(len<=0)return c;c.lo.resize(len);c.hi.resize(len);
 int an=a.lo.size(),bn=b.lo.size();
 #pragma omp parallel for schedule(static)
 for(int z=0;z<len;z++){
  int first=std::max(0,z-bn+1),last=std::min(an-1,z);Z sl=0,sh=0;
  for(int t=first;t<=last;t++){sl+=a.lo[t]*b.lo[z-t];sh+=a.hi[t]*b.hi[z-t];}
  c.lo[z]=floor_scale(sl);c.hi[z]=ceil_scale(sh);
 }
 return c;
}
IV at(const Poly&p,int z){z-=p.off;return z>=0&&z<int(p.lo.size())?IV{p.lo[z],p.hi[z]}:IV{0,0};}
uint64_t choose(int n,int r){
 uint64_t x=1;for(int j=1;j<=r;j++){
  uint64_t factor=uint64_t(n-j+1);
  if(factor && x>std::numeric_limits<uint64_t>::max()/factor)throw std::runtime_error("binomial overflow");
  x=(x*factor)/uint64_t(j);
 }return x;
}
struct Powers{std::vector<Poly>sp,bp;std::vector<IV>norm;};
Powers powers(const Data&d,int i,int j){
 Poly qs{0,V(d.D),V(d.D)},qb{d.D,V(d.B[1]-d.D),V(d.B[1]-d.D)};
 for(int a=0;a<d.B[1];a++){
  Z z=d.g[i][a]*d.g[j][a],den=Z(d.den)*ONE;IV q{floor_div(z,den),ceil_div(z,den)};
  if(a<d.D){qs.lo[a]=q.lo;qs.hi[a]=q.hi;}else{qb.lo[a-d.D]=q.lo;qb.hi[a-d.D]=q.hi;}
 }
 Powers P;P.sp.resize(d.k+1);P.bp.resize(d.maxr+1);P.norm.resize(d.n);
 Poly cur{0,V{ONE},V{ONE}};
 for(int r=1;r<=d.k;r++){cur=convolve(cur,qs,d.n-1);if(r>=d.k-d.maxr-1)P.sp[r]=cur;}
 P.bp[0]=Poly{0,V{ONE},V{ONE}};
 for(int r=1;r<=d.maxr;r++)P.bp[r]=convolve(P.bp[r-1],qb,d.B[r]-r);
 for(int r=0;r<=d.maxr;r++){
  Poly t=convolve(P.sp[d.k-r],P.bp[r],d.n-1);uint64_t mult=choose(d.k,r);
  for(int s=0;s<int(t.lo.size());s++){P.norm[t.off+s].lo+=mult*t.lo[s];P.norm[t.off+s].hi+=mult*t.hi[s];}
 }
 return P;
}
IV norm_integral(const Data&d,const Powers&P,int i,int j){
 IV out;for(int s=0;s<d.n;s++){
  Z z=d.h[i][s]*d.h[j][s];IV p{floor_scale(z),ceil_scale(z)};
  IV q=pos_product(P.norm[s],p);out.lo+=q.lo;out.hi+=q.hi;
 }return out;
}
// Group the common-coordinate weights by their fiber cutoff.
void weights(std::vector<IV>&W,const Data&d,const Powers&P,int s,int last){
 std::fill(W.begin(),W.begin()+last+1,IV{0,0});
 for(int r=0;r<=d.maxr;r++){
  int lastb=std::min(s,int(P.bp[r].off+P.bp[r].lo.size())-1);uint64_t mult=choose(d.k-1,r);
  for(int b=P.bp[r].off;b<=lastb;b++){
   int u=std::min(d.n-1-s,std::max(d.D-1,d.B[r+1]-r-1-b));
   if(u<0||u>last)throw std::runtime_error("weight subscript");
   IV a=at(P.bp[r],b),v=at(P.sp[d.k-1-r],s-b);
   W[u].lo+=mult*(a.lo*v.lo);W[u].hi+=mult*(a.hi*v.hi);
  }
 }
 for(int u=0;u<=last;u++)W[u]={floor_scale(W[u].lo),ceil_scale(W[u].hi)};
}
// Products g*h are accumulated exactly at scale 2^192 before rounding.
void prefixes(std::vector<IV>&pi,std::vector<IV>&pj,const Data&d,int i,int j,int s,int last){
 Z si=0,sj=0;for(int u=0;u<=last;u++){
  si+=d.g[i][u]*d.h[i][s+u];sj+=d.g[j][u]*d.h[j][s+u];
  pi[u]={floor_scale(si),ceil_scale(si)};pj[u]={floor_scale(sj),ceil_scale(sj)};
 }
}
IV marginal_integral(const Data&d,const Powers&P,int i,int j){
 int maxs=std::min(d.n-1,d.L-(d.k-1));std::vector<IV>totals(omp_get_max_threads());
 #pragma omp parallel
 {
  int th=omp_get_thread_num(),nu=d.B[1];std::vector<IV>W(nu),pi(nu),pj(nu);IV sum;
  #pragma omp for schedule(static)
  for(int s=0;s<=maxs;s++){
   int last=std::min(d.n-1-s,nu-1);weights(W,d,P,s,last);prefixes(pi,pj,d,i,j,s,last);
   for(int u=0;u<=last;u++){
    IV v=pos_product(W[u],product(pi[u],pj[u]));sum.lo+=v.lo;sum.hi+=v.hi;
   }
  }
  totals[th]=sum;
 }
 IV out;for(const auto&v:totals){out.lo+=v.lo;out.hi+=v.hi;}
 Z den=Z(d.den)*d.den;out={floor_div(out.lo,den),ceil_div(out.hi,den)};
 if(out.lo>out.hi)throw std::runtime_error("invalid marginal enclosure");return out;
}
int main(int argc,char**argv){try{
 if(argc!=4)throw std::runtime_error("usage: certify_integer witness.fix i j");
 Data d=load(argv[1]);int i=std::stoi(argv[2]),j=std::stoi(argv[3]);
 if(i<0||j<i||j>=d.R)throw std::runtime_error("invalid pair");
 double start=omp_get_wtime();Powers P=powers(d,i,j);
 IV I=norm_integral(d,P,i,j),J=marginal_integral(d,P,i,j);
 if(I.lo>I.hi||J.lo>J.hi)throw std::runtime_error("invalid interval");
 std::cout<<"{\"i\":"<<i<<",\"j\":"<<j<<",\"bits\":"<<PREC
 <<",\"Ilo\":\""<<I.lo<<"\",\"Ihi\":\""<<I.hi
 <<"\",\"Jlo\":\""<<J.lo<<"\",\"Jhi\":\""<<J.hi
 <<"\",\"seconds\":"<<(omp_get_wtime()-start)<<"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}return 0;
}
