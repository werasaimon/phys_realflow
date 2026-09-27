// Tensors and the Einstein convention: the storage, the index operations (outer product, raising
// and lowering, permutation, symmetrization, trace) and einstein(spec, ...), which reads the
// formula as text, checks it against the rule "a repeated letter is once up and once down", and
// sums by walking every value of every letter. See Tensor.h.
#include "math/Tensor.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

namespace rf {

namespace {

[[noreturn]] void fail(const std::string& ru, const std::string& en) {
    throw std::invalid_argument("einstein/tensor: " + ru + " / " + en);
}

size_t power(int base, int exponent) {
    size_t p = 1;
    for (int i = 0; i < exponent; ++i) p *= size_t(base);
    return p;
}

// Advances the multi-index idx[0..n) through 0 .. dim-1 like an odometer (the last digit fastest).
// False when it wrapped round: every combination has been visited.
bool nextIndex(std::vector<int>& idx, int dim) {
    for (int k = int(idx.size()) - 1; k >= 0; --k) {
        if (++idx[size_t(k)] < dim) return true;
        idx[size_t(k)] = 0;
    }
    return false;
}

} // namespace

// ---------------------------------------------------------------------------------------------
// Storage and arithmetic
// ---------------------------------------------------------------------------------------------
Tensor::Tensor(int dim, const std::string& variance) : dim_(dim), variance_(variance) {
    for (char c : variance)
        if (c != '^' && c != '_') fail("вариантность из символов ^ и _, а не '" + variance + "'",
                                       "variance is made of ^ and _, not '" + variance + "'");
    data_.assign(power(dim, int(variance.size())), 0.0);
}

size_t Tensor::offset(const int* idx) const {
    size_t o = 0;
    for (int s = 0; s < rank(); ++s) o = o * size_t(dim_) + size_t(idx[s]);
    return o;
}

double Tensor::scalar() const {
    if (rank() != 0) fail("scalar() у тензора ранга " + std::to_string(rank()), "scalar() of a rank " + std::to_string(rank()) + " tensor");
    return data_[0];
}

static void requireSameShape(const Tensor& a, const Tensor& b) {
    if (a.dim() != b.dim() || a.variance() != b.variance())
        fail("складывать можно тензоры одной размерности и вариантности: " + a.variance() + " и " + b.variance(),
             "only tensors of the same dimension and variance add: " + a.variance() + " and " + b.variance());
}

Tensor Tensor::operator+(const Tensor& b) const {
    requireSameShape(*this, b);
    Tensor r = *this;
    for (size_t k = 0; k < data_.size(); ++k) r.data_[k] += b.data_[k];
    return r;
}

Tensor Tensor::operator-(const Tensor& b) const {
    requireSameShape(*this, b);
    Tensor r = *this;
    for (size_t k = 0; k < data_.size(); ++k) r.data_[k] -= b.data_[k];
    return r;
}

Tensor Tensor::operator*(double s) const {
    Tensor r = *this;
    for (double& v : r.data_) v *= s;
    return r;
}

double Tensor::maxAbs() const {
    double m = 0;
    for (double v : data_) m = std::max(m, std::fabs(v));
    return m;
}

// ---------------------------------------------------------------------------------------------
// Index operations, each written with einstein() where that reads like the formula
// ---------------------------------------------------------------------------------------------
Tensor outer(const Tensor& a, const Tensor& b) {
    if (a.dim() != b.dim()) fail("разные размерности", "different dimensions");
    Tensor r(a.dim(), a.variance() + b.variance());
    size_t k = 0;
    for (double x : a.data())
        for (double y : b.data()) r.data()[k++] = x * y;
    return r;
}

// The letters a, b, c, ... for the slots of t, with its own variance markers: "^a_b_c".
static std::string slotLetters(const Tensor& t, int renamedSlot, char renamedLetter) {
    std::string s;
    for (int i = 0; i < t.rank(); ++i) {
        s += t.variance()[size_t(i)];
        s += i == renamedSlot ? renamedLetter : char('a' + i);
    }
    return s;
}

static std::string resultLetters(const Tensor& t, int slot, char newVariance) {
    std::string s;
    for (int i = 0; i < t.rank(); ++i) {
        s += i == slot ? newVariance : t.variance()[size_t(i)];
        s += char('a' + i);
    }
    return s;
}

// T^..a.. = g^{a z} T_..z..: the slot's letter becomes z, contracted with the inverse metric.
Tensor raise(const Tensor& t, int slot, const Tensor& ginv) {
    if (t.isUpper(slot)) fail("индекс уже сверху", "the index is already up");
    const std::string spec = slotLetters(t, slot, 'z') + " ^" + std::string(1, char('a' + slot)) + "z -> " + resultLetters(t, slot, '^');
    return einstein(spec, t, ginv);
}

Tensor lower(const Tensor& t, int slot, const Tensor& g) {
    if (!t.isUpper(slot)) fail("индекс уже снизу", "the index is already down");
    const std::string spec = slotLetters(t, slot, 'z') + " _" + std::string(1, char('a' + slot)) + "z -> " + resultLetters(t, slot, '_');
    return einstein(spec, t, g);
}

Tensor permute(const Tensor& t, const std::string& from, const std::string& to) {
    if (int(from.size()) != t.rank() || to.size() != from.size())
        fail("permute: имена слотов не совпадают с рангом", "permute: slot names do not match the rank");
    std::string in, out;
    for (int i = 0; i < t.rank(); ++i) {
        in += t.variance()[size_t(i)];
        in += from[size_t(i)];
    }
    for (char c : to) {
        const size_t i = from.find(c);
        if (i == std::string::npos) fail(std::string("permute: нет слота ") + c, std::string("permute: no slot ") + c);
        out += t.variance()[i];
        out += c;
    }
    return einstein(in + " -> " + out, t);
}

// T with slots s1 and s2 swapped (same variance required).
static Tensor swapSlots(const Tensor& t, int s1, int s2) {
    if (t.variance()[size_t(s1)] != t.variance()[size_t(s2)])
        fail("симметризовать можно только слоты одной вариантности", "only slots of the same variance can be symmetrized");
    std::string from, to;
    for (int i = 0; i < t.rank(); ++i) from += char('a' + i);
    to = from;
    std::swap(to[size_t(s1)], to[size_t(s2)]);
    return permute(t, from, to);
}

Tensor symmetrize(const Tensor& t, int s1, int s2) { return (t + swapSlots(t, s1, s2)) * 0.5; }
Tensor antisymmetrize(const Tensor& t, int s1, int s2) { return (t - swapSlots(t, s1, s2)) * 0.5; }

Tensor trace(const Tensor& t, int s1, int s2) {
    std::string in, out;
    for (int i = 0; i < t.rank(); ++i) {
        const char letter = (i == s1 || i == s2) ? 'z' : char('a' + i);
        in += t.variance()[size_t(i)];
        in += letter;
        if (letter != 'z') {
            out += t.variance()[size_t(i)];
            out += letter;
        }
    }
    return einstein(in + " -> " + out, t); // einstein() refuses two uppers or two lowers
}

// ---------------------------------------------------------------------------------------------
// The Einstein convention
// ---------------------------------------------------------------------------------------------
namespace {

struct Slot {
    char letter;
    bool up;
};
using Token = std::vector<Slot>;

// "^a_bc" -> a up, b down, c down. Every letter must follow a marker.
Token parseToken(const std::string& text) {
    Token t;
    char marker = 0;
    for (char c : text) {
        if (c == '^' || c == '_') { marker = c; continue; }
        if (!std::isalpha(static_cast<unsigned char>(c)))
            fail(std::string("в формуле допустимы буквы, ^ и _, а не '") + c + "'", std::string("only letters, ^ and _ are allowed, not '") + c + "'");
        if (!marker) fail(std::string("перед индексом ") + c + " нет ^ или _", std::string("index ") + c + " has no ^ or _ before it");
        t.push_back({c, marker == '^'});
    }
    return t;
}

void parseSpec(const std::string& spec, std::vector<Token>& inputs, Token& output) {
    const size_t arrow = spec.find("->");
    if (arrow == std::string::npos) fail("в формуле нет '->'", "the formula has no '->'");
    std::istringstream left(spec.substr(0, arrow));
    std::string word;
    while (left >> word) inputs.push_back(parseToken(word));
    std::istringstream right(spec.substr(arrow + 2));
    std::string rest, all;
    while (right >> rest) all += rest;
    output = parseToken(all);
}

// The rule itself: a letter once is free (and must be in the result with the same variance),
// twice is summed (once up, once down), more often is an error.
void checkRule(const std::vector<Token>& inputs, const Token& output) {
    std::string seen;
    for (const Token& t : inputs)
        for (const Slot& s : t) seen += s.letter;
    for (char c : std::string(seen)) {
        int ups = 0, downs = 0;
        for (const Token& t : inputs)
            for (const Slot& s : t)
                if (s.letter == c) (s.up ? ups : downs)++;
        const std::string L(1, c);
        if (ups + downs > 2) fail("индекс " + L + " встречается больше двух раз", "index " + L + " appears more than twice");
        if (ups == 2) fail("индекс " + L + " дважды сверху: свёртка по Эйнштейну требует одного верхнего и одного нижнего; опустите индекс метрикой",
                           "index " + L + " is up twice: Einstein summation needs one upper and one lower index; lower one with the metric");
        if (downs == 2) fail("индекс " + L + " дважды снизу: свёртка по Эйнштейну требует одного верхнего и одного нижнего; поднимите индекс обратной метрикой",
                             "index " + L + " is down twice: Einstein summation needs one upper and one lower index; raise one with the inverse metric");
        const bool inOutput = std::any_of(output.begin(), output.end(), [&](const Slot& s) { return s.letter == c; });
        if (ups + downs == 2 && inOutput) fail("индекс " + L + " суммируется и не может стоять в результате", "index " + L + " is summed and cannot be in the result");
        if (ups + downs == 1 && !inOutput) fail("свободный индекс " + L + " пропущен в результате", "free index " + L + " is missing from the result");
    }
    for (const Slot& s : output) {
        int count = 0;
        bool up = false;
        for (const Token& t : inputs)
            for (const Slot& q : t)
                if (q.letter == s.letter) { ++count; up = q.up; }
        const std::string L(1, s.letter);
        if (count != 1) fail("индекс результата " + L + " должен быть свободным индексом входа", "result index " + L + " must be a free index of the inputs");
        if (up != s.up) fail("индекс " + L + " меняет положение (верх/низ) в результате: сначала поднимите или опустите его метрикой",
                             "index " + L + " changes position (up/down) in the result: raise or lower it with the metric first");
    }
}

// For one tensor: which letter each slot reads, and the slot's stride in the flat storage.
struct Access {
    std::vector<int> letter;
    std::vector<size_t> stride;
};

Access makeAccess(const Token& t, const std::string& letters, int dim) {
    Access a;
    size_t stride = 1;
    a.letter.resize(t.size());
    a.stride.resize(t.size());
    for (int s = int(t.size()) - 1; s >= 0; --s) {
        a.letter[size_t(s)] = int(letters.find(t[size_t(s)].letter));
        a.stride[size_t(s)] = stride;
        stride *= size_t(dim);
    }
    return a;
}

size_t offsetOf(const Access& a, const std::vector<int>& value) {
    size_t o = 0;
    for (size_t s = 0; s < a.letter.size(); ++s) o += size_t(value[size_t(a.letter[s])]) * a.stride[s];
    return o;
}

std::string varianceOf(const Token& t) {
    std::string v;
    for (const Slot& s : t) v += s.up ? '^' : '_';
    return v;
}

} // namespace

Tensor einstein(const std::string& spec, const std::vector<const Tensor*>& inputs) {
    std::vector<Token> tokens;
    Token output;
    parseSpec(spec, tokens, output);
    if (tokens.size() != inputs.size())
        fail("в формуле " + std::to_string(tokens.size()) + " тензоров, передано " + std::to_string(inputs.size()),
             "the formula names " + std::to_string(tokens.size()) + " tensors, " + std::to_string(inputs.size()) + " given");
    const int dim = inputs.empty() ? 0 : inputs[0]->dim();
    for (size_t i = 0; i < tokens.size(); ++i) {
        if (inputs[i]->variance() != varianceOf(tokens[i]))
            fail("тензор " + std::to_string(i + 1) + " имеет вариантность " + inputs[i]->variance() + ", а в формуле " + varianceOf(tokens[i]),
                 "tensor " + std::to_string(i + 1) + " has variance " + inputs[i]->variance() + ", the formula says " + varianceOf(tokens[i]));
        if (inputs[i]->dim() != dim) fail("тензоры разной размерности", "tensors of different dimensions");
    }
    checkRule(tokens, output);
    std::string letters; // the result's letters first, then the summed ones
    for (const Slot& s : output) letters += s.letter;
    for (const Token& t : tokens)
        for (const Slot& s : t)
            if (letters.find(s.letter) == std::string::npos) letters += s.letter;
    std::vector<Access> access;
    for (const Token& t : tokens) access.push_back(makeAccess(t, letters, dim));
    const Access outAccess = makeAccess(output, letters, dim);
    Tensor result(dim, varianceOf(output));
    std::vector<int> value(letters.size(), 0);
    do { // every value of every letter: the sum over the repeated ones happens by adding up
        double product = 1.0;
        for (size_t i = 0; i < inputs.size(); ++i) product *= inputs[i]->data()[offsetOf(access[i], value)];
        result.data()[offsetOf(outAccess, value)] += product;
    } while (nextIndex(value, dim));
    return result;
}

Tensor einstein(const std::string& spec, const Tensor& a) { return einstein(spec, std::vector<const Tensor*>{&a}); }
Tensor einstein(const std::string& spec, const Tensor& a, const Tensor& b) { return einstein(spec, std::vector<const Tensor*>{&a, &b}); }
Tensor einstein(const std::string& spec, const Tensor& a, const Tensor& b, const Tensor& c) {
    return einstein(spec, std::vector<const Tensor*>{&a, &b, &c});
}
Tensor einstein(const std::string& spec, const Tensor& a, const Tensor& b, const Tensor& c, const Tensor& d) {
    return einstein(spec, std::vector<const Tensor*>{&a, &b, &c, &d});
}

} // namespace rf
