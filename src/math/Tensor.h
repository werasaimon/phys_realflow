#pragma once
// Tensors of any rank in any dimension, with the Einstein summation convention - written the way a
// textbook writes them, without template magic.
//
// A tensor is a machine that takes directions (vectors) and rulers (covectors) and returns a
// number, linear in each of them. Its components sit in a table with one index per slot; each
// slot is either UPPER (^, contravariant: the slot holds an arrow, a vector component) or LOWER
// (_, covariant: the slot holds a ruler, a gradient component). The variance string names them
// slot by slot: "^__" is Gamma^l_mn, "__" the metric g_mn, "^^" its inverse g^mn.
//
//   Tensor g(4, "__");                         // 4 x 4 zeros, both indices down
//   g(0, 0) = -1; g(1, 1) = g(2, 2) = g(3, 3) = 1;
//   Tensor v(4, "^");  v(0) = 1;               // a vector
//   double n = einstein("_ab ^a ^b ->", g, v, v).scalar();   // g_ab v^a v^b = -1
//
// The Einstein convention (Einstein 1916, "Die Grundlage der allgemeinen Relativitaetstheorie",
// par. 5): a letter written twice in a product, once up and once down, is summed over all values.
// einstein(spec, ...) takes the formula as text: "^a_bcd ^b ^c ^d -> ^a" is R^a_bcd u^b x^c u^d.
// A marker ^ or _ applies to the letters after it until the next marker; spaces separate the
// tensors; after "->" come the free indices of the result (nothing: a scalar). A letter twice up
// (or twice down) is refused with a message saying why - that is the whole point of the rule:
// only an arrow can be measured with a ruler.
//
// Storage: the components in one std::vector, the first index slowest (row-major), dim^rank of
// them. This is an ANALYSIS layer (curvature of a metric, checks, teaching), not a simulation hot
// path: it allocates freely and reports a misuse by throwing std::invalid_argument with a message
// in Russian and English. The solvers never use it.
#include <stdexcept>
#include <string>
#include <vector>

namespace rf {

class Tensor {
public:
    Tensor() = default;
    // dim^rank zeros; rank = variance.size(), each character '^' or '_'.
    Tensor(int dim, const std::string& variance);

    int dim() const { return dim_; }
    int rank() const { return int(variance_.size()); }
    const std::string& variance() const { return variance_; }
    bool isUpper(int slot) const { return variance_[size_t(slot)] == '^'; }
    size_t size() const { return data_.size(); }
    std::vector<double>& data() { return data_; }
    const std::vector<double>& data() const { return data_; }

    // Components by index: as many ints as the rank (0 for a scalar).
    double& operator()() { return at(nullptr); }
    double& operator()(int i) { const int a[1] = {i}; return at(a); }
    double& operator()(int i, int j) { const int a[2] = {i, j}; return at(a); }
    double& operator()(int i, int j, int k) { const int a[3] = {i, j, k}; return at(a); }
    double& operator()(int i, int j, int k, int l) { const int a[4] = {i, j, k, l}; return at(a); }
    double operator()() const { return at(nullptr); }
    double operator()(int i) const { const int a[1] = {i}; return at(a); }
    double operator()(int i, int j) const { const int a[2] = {i, j}; return at(a); }
    double operator()(int i, int j, int k) const { const int a[3] = {i, j, k}; return at(a); }
    double operator()(int i, int j, int k, int l) const { const int a[4] = {i, j, k, l}; return at(a); }
    // Any rank: idx holds rank() indices.
    double& at(const int* idx) { return data_[offset(idx)]; }
    double at(const int* idx) const { return data_[offset(idx)]; }
    double scalar() const; // the value of a rank-0 tensor

    // Sums need the same dimension and the same variance: an arrow plus a ruler means nothing.
    Tensor operator+(const Tensor& b) const;
    Tensor operator-(const Tensor& b) const;
    Tensor operator*(double s) const;
    double maxAbs() const; // the largest |component|

private:
    size_t offset(const int* idx) const;
    int dim_ = 0;
    std::string variance_;
    std::vector<double> data_;
};

// The tensor product (a outer b)_{...} = a_{...} b_{...}: the slots of a, then those of b.
Tensor outer(const Tensor& a, const Tensor& b);
// Raise slot `slot` with the inverse metric ginv ("^^"): T^a = g^ab T_b. Lower it with g ("__").
Tensor raise(const Tensor& t, int slot, const Tensor& ginv);
Tensor lower(const Tensor& t, int slot, const Tensor& g);
// The slots renamed and reordered: permute(t, "abc", "acb") swaps the last two slots (with their
// variance). `from` names t's slots, `to` is the order of the result.
Tensor permute(const Tensor& t, const std::string& from, const std::string& to);
// (T_ab + T_ba) / 2 and (T_ab - T_ba) / 2 over slots s1, s2 (same variance), other slots kept.
Tensor symmetrize(const Tensor& t, int s1, int s2);
Tensor antisymmetrize(const Tensor& t, int s1, int s2);
// The trace over one upper and one lower slot: T^a_a (both slots disappear).
Tensor trace(const Tensor& t, int s1, int s2);

// The Einstein convention: see the file head. The inputs in the order of the spec's tokens.
Tensor einstein(const std::string& spec, const std::vector<const Tensor*>& inputs);
Tensor einstein(const std::string& spec, const Tensor& a);
Tensor einstein(const std::string& spec, const Tensor& a, const Tensor& b);
Tensor einstein(const std::string& spec, const Tensor& a, const Tensor& b, const Tensor& c);
Tensor einstein(const std::string& spec, const Tensor& a, const Tensor& b, const Tensor& c, const Tensor& d);

} // namespace rf
