#pragma once
// The small LCP of one contact patch: at most four points pushing between two bodies (or one body
// and a support that does not give). Find the normal impulses x with
//
//     K x + b >= 0,     x >= 0,     x . (K x + b) = 0
//
// where K is the effective-mass matrix of the points (how an impulse at point j changes the normal
// velocity at point i) and b the normal velocities the points would have with no impulse at all.
// In words: every point either pushes (x > 0) and ends at exactly the velocity it aims for, or does
// not push (x = 0) and is not being pulled into contact.
//
// Solved exactly by total enumeration of the active sets - the guesses of which points push - as
// Box2D's two-point block solver (Catto, "b2ContactSolver", 2009), generalised to four points: each
// guess is a small linear system, and the right guess is the one whose answer is consistent. With a
// small CFM on the diagonal K is positive definite, the solution unique, and the search order only
// changes how fast it is found. Used by the two-sided block solve (ContactSolver.cpp) and by the
// one-sided solve of shock propagation (ShockPropagation.cpp).
#include "math/MatrixNxN.h"

#include <bitset>

namespace rf {

struct ContactLcp {
    int n = 0;
    float bb[4] = {0, 0, 0, 0};   // the velocities the points would have with no impulse at all
    float Kc[4][4] = {};          // the effective-mass matrix with the CFM on its diagonal
    float x[4] = {0, 0, 0, 0};    // the solution found by the last successful tryActiveSet

    // Solves the points of `mask` as equalities and checks the complementarity of the rest.
    bool tryActiveSet(int mask) {
        int idx[4], k = 0;
        for (int i = 0; i < n; ++i)
            if (mask & (1 << i)) idx[k++] = i;
        float xs4[4] = {0, 0, 0, 0};
        if (k > 0) {
            float Ms[4][4], rs[4], xs[4];
            for (int i = 0; i < k; ++i) {
                rs[i] = -bb[idx[i]];
                for (int j = 0; j < k; ++j) Ms[i][j] = Kc[idx[i]][idx[j]];
            }
            if (!solveSmall(k, Ms, rs, xs)) return false;
            for (int i = 0; i < k; ++i)
                if (xs[i] < 0.0f) return false;
            for (int i = 0; i < k; ++i) xs4[idx[i]] = xs[i];
        }
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) continue;
            float w = bb[i];
            for (int j = 0; j < n; ++j) w += Kc[i][j] * xs4[j];
            if (w < -1e-5f) return false;
        }
        for (int i = 0; i < n; ++i) x[i] = xs4[i];
        return true;
    }

    // The search: first nobody pushes (the whole patch separates - common for speculative points),
    // then the guess `hint` (last time's active set: it rarely changes), then every set, largest
    // first. Returns the active set found, or -1 (a numerical corner case: then x = 0).
    int solve(int hint) {
        const int full = (1 << n) - 1;
        if (tryActiveSet(0)) return 0;
        if (hint > 0 && hint <= full && tryActiveSet(hint)) return hint;
        for (int size = n; size >= 1; --size)
            for (int mask = full; mask >= 1; --mask) {
                if (int(std::bitset<32>(unsigned(mask)).count()) != size || mask == hint) continue;
                if (tryActiveSet(mask)) return mask;
            }
        for (int i = 0; i < n; ++i) x[i] = 0;
        return -1;
    }
};

} // namespace rf
