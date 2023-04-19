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

gmpq getOrZero(std::map<int, gmpq> mp, int k) {
  if(mp.contains(k)) {
    return mp.at(k);
  } else {
    return 0;
  }
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
    for(int k : rpl) {
      insertWithPlus(trms, k, minusrat);
    }
  }
}

std::vector<int> includeMods(int n, int q, int start) {
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

// [[Rcpp::export]]
std::vector<int> Rcpp_removeExps(int n, int p, int q) {
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

std::optional<gmpq> equalReplacements(int p, int r, cyclotomic cyc) {
  int ord = cyc.order;
  std::map<int, gmpq> trms = cyc.terms;
  std::vector<int> rpl = replacements(ord, p, r);
  gmpq x1 = getOrZero(trms, rpl[0]);
  int l = rpl.size();
  gmpq x;
  for(int i = 1; i < l; i++) {
    x = getOrZero(trms, rpl[i]);
    if(x != x1) {
      return std::nullopt;
    }
  }
  return x1;
}

cyclotomic reduceByPrime(int p, cyclotomic cyc) {
  int n = cyc.order;
  std::vector<gmpq> cfs(0);
  std::optional<gmpq> x = equalReplacements(p, 0, cyc);
  int r = p;
  while(r < n - p && x) {
    cfs.push_back(-(*x));
    x = equalReplacements(p, r, cyc);
    r += p;
  }
  if(x) {
    int ndivp = n / p;
    std::map<int, gmpq> trms;
    int i = 0;
    int l = cfs.size();
    while(i < ndivp && i < l) {
      gmpq coef = cfs[i];
      if(coef != 0) {
        trms[i] = coef;
      }
      i++;
    }
    cyclotomic out;
    out.order = ndivp;
    out.terms = trms;
    return out;
  } else {
    return cyc;
  }
}

void removeZeros(std::map<int, gmpq> mp) {
  std::map<int, gmpq>::iterator it = mp.begin();
  while(it != mp.end()) {
    if(it->second == 0) {
      it = mp.erase(it);
    } else {
      ++it;
    }
  }
}

int gcdCyc(cyclotomic cyc) {
  Rcpp::Function f("R_gcdList");
  std::map<int, gmpq> trms = cyc.terms;
  int l = trms.size();
  Rcpp::IntegerVector x(l+1);
  x(0) = cyc.order;
  int i = 1;
  for(const auto& item : trms) {
    x(i++) = item.first;
  }
  return Rcpp::as<int>(f(x));
}

cyclotomic gcdReduce(cyclotomic cyc) {
  int d = gcdCyc(cyc);
  if(d == 1) {
    return cyc;
  }
  std::map<int, gmpq> trms = cyc.terms;
  std::map<int, gmpq> newtrms;
  for(const auto& item : trms) {
    int key = item.first;
    newtrms[key / d] = item.second;
  }
  cyclotomic out;
  out.order = cyc.order / d;
  out.terms = newtrms;
  return out;
}

std::optional<gmpq> equalCoefficients(cyclotomic cyc) {
  std::map<int, gmpq> trms = cyc.terms;
  if(trms.size() == 0) {
    return std::nullopt;
  }
  std::map<int, gmpq>::iterator it = trms.begin();
  gmpq firstcoef = it->second;
  while(it != trms.end()) {
    it++;
    if(it->second != firstcoef) {
      return std::nullopt;
    }
  }
  return firstcoef;
}

int intpow(int base, unsigned exp){
  int result = 1;
  while(exp){
    if(exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}

void convertToBase(int n, std::map<int, gmpq> trms) {
  if(n > 1) {
    Rcpp::Function f("R_extraneousPowers");
    Rcpp::IntegerMatrix epows = Rcpp::as<Rcpp::IntegerMatrix>(f(n));
    int l = epows.nrow();
    for(int i = l-1; i = 0; i--) {
      Rcpp::IntegerVector pr = epows(i, Rcpp::_);
      replace(n, pr(1), pr(2), trms);
    }
  }
}
