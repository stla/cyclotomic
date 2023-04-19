#include <Rcpp.h>
#include <boost/multiprecision/gmp.hpp>
typedef boost::multiprecision::mpq_rational                 gmpq;

struct cyclotomic {
  int order;
  std::map<int, gmpq> terms;
};

void insertWithPlus(std::map<int, gmpq>& mp, int k, gmpq v) {
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

void replace(int n, int p, int r, std::map<int, gmpq>& trms) {
  if(trms.contains(r)) {
    gmpq minusrat = -trms.at(r);
    trms.erase(r);
    std::vector<int> rpl = replacements(n, p, r);
    for(auto& k : rpl) {
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

void append(std::vector<int>& v1, std::vector<int> v2) {
  v1.reserve(v1.size() + v2.size());
  for(const auto& k : v2) {
    v1.emplace_back(k);
  }
}

// [[Rcpp::export]]
std::vector<int> Rcpp_removeExps(int n, int p, int q) {
  int ndivq = n / q;
  std::vector<int> out(0);
  if(p == 2) {
    int x = q / 2;
    append(out, includeMods(n, q, ndivq * x));
    while(x < q - 1) {
      x++;
      append(out, includeMods(n, q, ndivq * x));
    }
  } else {
    int m = (q / p - 1) / 2;
    int x = -m;
    append(out, includeMods(n, q, ndivq * x));
    while(x < m) {
      x++;
      append(out, includeMods(n, q, ndivq * x));
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

void removeZeros(std::map<int, gmpq>& mp) {
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
  ++it;
  while(it != trms.end()) {
    if(it->second != firstcoef) {
      return std::nullopt;
    }
    ++it;
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

void convertToBase(int n, std::map<int, gmpq>& trms) {
  if(n > 1) {
    Rcpp::Function f("R_extraneousPowers");
    Rcpp::IntegerMatrix epows = Rcpp::as<Rcpp::IntegerMatrix>(f(n));
    int l = epows.nrow();
    for(int i = l-1; i >= 0; i--) {
      Rcpp::IntegerVector pr = epows(i, Rcpp::_);
      replace(n, pr(0), pr(1), trms);
    }
  }
}

cyclotomic fromRational(gmpq rat) {
  std::map<int, gmpq> trms;
  if(rat != 0) {
    trms[0] = rat;
  }
  cyclotomic cyc;
  cyc.order = 1;
  cyc.terms = trms;
  return cyc;
}

cyclotomic tryRational(cyclotomic cyc) {
  Rcpp::Function f("R_phiNrpSqfree");
  Rcpp::List pns = Rcpp::as<Rcpp::List>(f(cyc.order));
  int phi     = Rcpp::as<int>(pns["phi"]);
  int nrp     = Rcpp::as<int>(pns["nrp"]);
  bool sqfree = Rcpp::as<bool>(pns["sqfree"]);
  if(sqfree && cyc.terms.size() == phi) {
    std::optional<gmpq> mayberat = equalCoefficients(cyc);
    if(mayberat) {
      gmpq rat = *mayberat;
      if(nrp % 2 == 0) {
        return fromRational(rat);
      } else {
        return fromRational(-rat);
      }
    } else {
      return cyc;
    }
  } else {
    return cyc;
  }
}

cyclotomic tryReduce(cyclotomic cyc) {
  Rcpp::Function f("R_squareFreeOddFactors");
  Rcpp::IntegerVector factors = Rcpp::as<Rcpp::IntegerVector>(f(cyc.order));
  int l = factors.size();
  if(l == 0) {
    return cyc;
  }
  for(int i = 0; i < l; i++) {
    cyc = reduceByPrime(factors(i), cyc);
  }
  return cyc;
}

cyclotomic cyclotomic0(int ord, std::map<int, gmpq> trms) {
  removeZeros(trms);
  cyclotomic cyc;
  cyc.order = ord;
  cyc.terms = trms;
  return tryReduce(tryRational(gcdReduce(cyc)));
}

cyclotomic mkCyclotomic(int ord, std::map<int, gmpq> trms) {
  convertToBase(ord, trms);
  return cyclotomic0(ord, trms);
}

cyclotomic zeta(int n) {
  std::map<int, gmpq> trms;
  gmpq one(1);
  if(n == 1) {
    trms[0] = one;
    cyclotomic cyc;
    cyc.order = 1;
    cyc.terms = trms;
    return cyc;
  }
  trms[1] = one;
  convertToBase(n, trms);
  return cyclotomic0(n, trms);
}

std::map<int, gmpq> unionWithPlus(
    std::map<int, gmpq> mp1, std::map<int, gmpq> mp2
) {
  for(auto& item : mp1) {
    int key = item.first;
    if(mp2.contains(item.first)) {
      mp1[key] += mp2[key];
      mp2.erase(key);
    }
  }
  for(auto& item : mp2) {
    int key = item.first;
    mp1[key] = mp2[key];
  }
  return mp1;
}

bool isZero(cyclotomic cyc) {
  return cyc.terms.size() == 0;
}

cyclotomic sumCyc(cyclotomic cyc1, cyclotomic cyc2) {
  if(isZero(cyc1)) {
    return cyc2;
  }
  if(isZero(cyc2)) {
    return cyc1;
  }
  int o1 = cyc1.order;
  std::map<int, gmpq> trms1 = cyc1.terms;
  int o2 = cyc2.order;
  std::map<int, gmpq> trms2 = cyc2.terms;
  Rcpp::Function scm("R_scm");
  int ord = Rcpp::as<int>(scm(o1, o2));
  int m1 = ord / o1;
  int m2 = ord / o2;
  std::map<int, gmpq> mp1;
  for(const auto& item : trms1) {
    int key = item.first;
    mp1[m1 * key] = item.second;
  }
  std::map<int, gmpq> mp2;
  for(const auto& item : trms2) {
    int key = item.first;
    mp2[m2 * key] = item.second;
  }
  std::map<int, gmpq> trms = unionWithPlus(mp1, mp2);
  return mkCyclotomic(ord, trms);
}

cyclotomic zeroCyc() {
  std::map<int, gmpq> trms;
  cyclotomic cyc;
  cyc.order = 1;
  cyc.terms = trms;
  return cyc;
}

cyclotomic prodCyc(cyclotomic cyc1, cyclotomic cyc2) {
  if(isZero(cyc1)) {
    return zeroCyc();
  }
  if(isZero(cyc2)) {
    return zeroCyc();
  }
  int o1 = cyc1.order;
  std::map<int, gmpq> trms1 = cyc1.terms;
  int o2 = cyc2.order;
  std::map<int, gmpq> trms2 = cyc2.terms;
  Rcpp::Function scm("R_scm");
  int ord = Rcpp::as<int>(scm(o1, o2));
  int m1 = ord / o1;
  int m2 = ord / o2;
  std::map<int, gmpq> trms;
  for(const auto& item1 : trms1) {
    int k1 = item1.first;
    gmpq c1 = item1.second;
    for(const auto& item2 : trms2) {
      int k2 = item2.first;
      gmpq c2 = item2.second;
      int k = (m1 * k1 + m2 * k2) % ord;
      insertWithPlus(trms, k, c1 * c2);
    }
  }
  return mkCyclotomic(ord, trms);
}

cyclotomic fromInteger(int n) {
  gmpq nrat(n);
  return fromRational(nrat);
}

cyclotomic prodRatCyc(gmpq rat, cyclotomic cyc) {
  if(rat == 0) {
    return zeroCyc();
  }
  std::map<int, gmpq> trms = cyc.terms;
  std::map<int, gmpq> newtrms;
  for(const auto& item : trms) {
    newtrms[item.first] = rat * item.second;
  }
  cyclotomic result;
  result.order = cyc.order;
  result.terms = newtrms;
  return result;
}

cyclotomic minusCyc(cyclotomic cyc) {
  gmpq o(-1);
  return prodRatCyc(o, cyc);
}

cyclotomic powerCyc(cyclotomic cyc, int p){ // TODO: p < 0 - requires invCyc
  cyclotomic result = fromInteger(1);
  if(p >= 0) {
    while(p) {
      if(p & 1) {
        result = prodCyc(result, cyc);
      }
      p >>= 1;
      cyc = prodCyc(cyc, cyc);
    }
    return result;
  }
  // if p < 0 TODO
  return zeroCyc();
}


//// | sqrt stuff ////////
//const cyclotomic zeta8 = zeta(8); impossible because of R_extraneousPowers not found!
//const cyclotomic sqrt2 = sumCyc(zeta8, minusCyc(powerCyc(zeta8, 3)));
//const cyclotomic im = zeta(4);

cyclotomic eb(int n) {
  if(n == 1) {
    return zeroCyc();
  }
  cyclotomic zetan = zeta(n);
  cyclotomic result = zetan;
  cyclotomic addme;
  int l = (n - 1) / 2; // =0 si n=2
  for(int k = 2; k <= l; k++) {
    addme = powerCyc(zetan, (k * k) % n);
    result = sumCyc(result, addme);
  }
  return result;
}















void display(std::map<int, gmpq>& mp) {
  for(auto& it : mp) {
    Rcpp::Rcout << "Key: " << it.first << ", value: " << it.second << "\n";
  }
}

// [[Rcpp::export]]
void test() {
  cyclotomic e4 = zeta(4);
  cyclotomic e9 = zeta(9);
  cyclotomic cyc = eb(8); //powerCyc(prodCyc(e4, e9), 3);
  Rcpp::Rcout << "Order: " << cyc.order << "\n";
  display(cyc.terms);
}
