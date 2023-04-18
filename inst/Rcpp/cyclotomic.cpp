// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/multiprecision/gmp.hpp>
typedef boost::multiprecision::mpq_rational                 gmpq;

struct cyclotomic {
  int order;
  std::map<int, gmpq> terms;
};

void insertWithPlus(std::map<int, gmpq> mp, int k, gmpq v) {
  if(mp.contains(k)) {
    v += mp.at(k);
  }
  mp[k] = v;
}

std::vector<int> replacements(int n, int p, int r) {
  int s = n / p;
  std::vector<int> rpl(0);
  int x = r - s;
  while(x >= 0) {
    rpl.push_back(x);
    x -= s;
  }
  x = r + s;
  while(x < n) {
    rpl.push_back(x);
    x += s;
  }
  return rpl;
}

void replace(int n, int p, int r, std::map<int, gmpq> trms) {
  if(trms.contains(r)) {
    gmpq minusrat = -trms.at(r);
    trms.erase(r);
    std::vector<int> rpl = replacements(n, p, r);
    int l = rpl.size();
    for(int k : rpl) {
      insertWithPlus(trms, k, minusrat);
    }
  }
}

std::vector<int> includeMods(n, q, start) {
  std::vector<int> out = {start};
  int x = start - q;
  while(x >= 0) {
    out.push_back(x);
    x -= q;
  }
  x = start + q;
  while(x < n) {
    out.push_back(x);
    x += q;
  }
  return out;
}

std::vector<int> removeExps(int n, int p, int q) {
  int ndivq = n / q;
  std::vector<int> out(0);
  if(p == 2) {
    int x = q / 2;
    out.push_back(ndivq * x);
    while(x < q - 1) {
      x++;
      out.push_back(ndivq * x);
    }
  } else {
    int m = (q / p - 1) / 2;
    int x = -m;
    out.push_back(ndivq * x);
    while(x < m) {
      x++;
      out.push_back(ndivq * x);
    }
  }
  return out;
}
