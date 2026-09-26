// Math and spatial structures against exact answers: vector/matrix/quaternion identities,
// the N x N solvers, mesh primitives (closed manifolds, volumes), BVH and AABB-tree queries
// against brute force.
#include "TestRunner.h"
#include "Tests.h"

void testMath() {
    auto near = [](float a, float b, float eps = 1e-5f) { return std::fabs(a - b) <= eps; };
    // Vector2 / Vector3 / Vector4
    CHECK(near(cross(Vector2(1, 0), Vector2(0, 1)), 1.0f) && near(dot(perp(Vector2(2, 3)), Vector2(2, 3)), 0.0f), "Vector2");
    Vector3 c = cross(Vector3(1, 0, 0), Vector3(0, 1, 0));
    CHECK(near(c.x, 0) && near(c.y, 0) && near(c.z, 1), "Vector3 cross");
    Vector3 p = anyPerpendicular(Vector3(0.3f, -2.0f, 0.7f));
    CHECK(near(dot(p, Vector3(0.3f, -2.0f, 0.7f)), 0.0f) && near(length(p), 1.0f), "anyPerpendicular");
    CHECK(near(length(Vector4(1, 2, 2, 4)), 5.0f) && near(Vector4(Vector3(1, 2, 3), 1).xyz().z, 3.0f), "Vector4");

    // Matrix3x3: inverse, determinant, skew, eigen
    Matrix3x3 A = Matrix3x3::fromRows({4, 1, 2}, {0, 3, 1}, {1, 0, 2});
    Matrix3x3 I = A * A.inverse();
    float err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) err = std::max(err, std::fabs(I.m[i][j] - (i == j ? 1.0f : 0.0f)));
    CHECK(err < 1e-5f && near(A.determinant(), 4 * 6 - 1 * (-1) + 2 * (-3)), "Matrix3x3 inverse/det (err %e, det %f)", err,
          A.determinant());
    Vector3 s = Matrix3x3::skew({1, 2, 3}) * Vector3(4, 5, 6), sc = cross(Vector3(1, 2, 3), Vector3(4, 5, 6));
    CHECK(near(s.x, sc.x) && near(s.y, sc.y) && near(s.z, sc.z), "skew * v == cross");
    CHECK(Matrix3x3::zero().inverse().m[0][0] == 0.0f, "singular inverse -> zero");
    Matrix3x3 S = Matrix3x3::fromRows({2, 1, 0}, {1, 3, 1}, {0, 1, 4}), V;
    Vector3 ev;
    symmetricEigen(S, ev, V);
    Matrix3x3 R = V * Matrix3x3::diag(ev) * V.transposed();
    err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) err = std::max(err, std::fabs(R.m[i][j] - S.m[i][j]));
    CHECK(err < 1e-5f, "3x3 eigen reconstruction error %e", err);

    // Quaternion: matrix round trip, log/exp, slerp, shortest arc
    Quaternion q = Quaternion::fromAxisAngle({1, 2, 3}, 1.1f);
    Quaternion q2 = Quaternion::fromMatrix3x3(q.toMatrix3x3());
    CHECK(near(std::fabs(dot(q, q2)), 1.0f), "quaternion <-> matrix");
    Quaternion q3 = Quaternion::fromRotationVector(q.log());
    CHECK(near(std::fabs(dot(q, q3)), 1.0f), "exp(log(q)) == q");
    CHECK(near(q.angle(), 1.1f), "angle %f", q.angle());
    Quaternion h = slerp(Quaternion(), Quaternion::fromAxisAngle({0, 0, 1}, 1.0f), 0.5f);
    CHECK(near(h.angle(), 0.5f), "slerp halfway angle %f", h.angle());
    Vector3 to = Quaternion::fromTwoVectors({1, 0, 0}, {0, 0, -1}).rotate({1, 0, 0});
    CHECK(near(to.z, -1.0f) && near(Quaternion::fromTwoVectors({1, 0, 0}, {-1, 0, 0}).rotate({1, 0, 0}).x, -1.0f),
          "fromTwoVectors");
    Quaternion qi = q * q.inverse();
    CHECK(near(qi.w, 1.0f) && near(qi.x, 0.0f), "q * q^-1 == 1");

    // Matrix4x4: TRS, inverse, projection
    Matrix4x4 M = Matrix4x4::trs({1, 2, 3}, q, {2, 3, 4});
    Vector3 x(0.3f, -0.7f, 1.1f);
    Vector3 back = M.inverse().transformPoint(M.transformPoint(x));
    CHECK(length(back - x) < 1e-5f, "Matrix4x4 inverse (err %e)", length(back - x));
    Matrix4x4 P = Matrix4x4::perspective(1.0f, 1.5f, 0.1f, 100.0f);
    CHECK(near(P.transformPoint({0, 0, -0.1f}).z, -1.0f, 1e-4f) && near(P.transformPoint({0, 0, -100}).z, 1.0f, 1e-3f),
          "perspective depth range");
    Vector3 eye(1, 2, 3);
    CHECK(length(Matrix4x4::lookAt(eye, {0, 0, 0}, {0, 1, 0}).transformPoint(eye)) < 1e-5f, "lookAt puts the eye at the origin");

    // MatrixNxN: LU, determinant, inverse, Cholesky, Jacobi (Hilbert-like SPD system, n = 6)
    const int n = 6;
    MatrixNxN H(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) H(i, j) = 1.0 / (i + j + 1) + (i == j ? 1.0 : 0.0);
    std::vector<double> xs(n), b, xl, xc;
    for (int i = 0; i < n; ++i) xs[i] = i - 2.5;
    b = H * xs;
    bool okL = H.solveLU(b, xl), okC = H.solveCholesky(b, xc);
    double eL = 0, eC = 0;
    for (int i = 0; i < n; ++i) eL = std::max(eL, std::fabs(xl[i] - xs[i])), eC = std::max(eC, std::fabs(xc[i] - xs[i]));
    CHECK(okL && okC && eL < 1e-12 && eC < 1e-12, "MatrixNxN solve: LU err %e, Cholesky err %e", eL, eC);
    MatrixNxN Hi;
    CHECK(H.inverse(Hi), "inverse");
    MatrixNxN E = H * Hi;
    double eI = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) eI = std::max(eI, std::fabs(E(i, j) - (i == j ? 1.0 : 0.0)));
    CHECK(eI < 1e-12, "H * H^-1 error %e", eI);
    MatrixNxN D = MatrixNxN::identity(3) * 2.0;
    D(0, 1) = 5;
    CHECK(std::fabs(D.determinant() - 8.0) < 1e-12, "determinant of triangular %f", D.determinant());
    MatrixNxN N(3, 3);
    N(0, 0) = 1; N(1, 1) = -1; // indefinite
    std::vector<double> dummy;
    CHECK(!N.solveCholesky({1, 1, 1}, dummy) && !MatrixNxN(3, 3).solveLU({1, 1, 1}, dummy), "Cholesky/LU must reject bad matrices");
    std::vector<double> vals;
    MatrixNxN Vn;
    H.symmetricEigen(vals, Vn);
    double eE = 0;
    for (int k = 0; k < n; ++k) { // H v_k = lambda_k v_k
        std::vector<double> v(n);
        for (int i = 0; i < n; ++i) v[i] = Vn(i, k);
        std::vector<double> hv = H * v;
        for (int i = 0; i < n; ++i) eE = std::max(eE, std::fabs(hv[i] - vals[k] * v[i]));
    }
    CHECK(eE < 1e-10, "N x N eigen residual %e", eE);

    // solveSmall (contact block solver)
    float Ms[4][4] = {{4, 1, 0, 0}, {1, 4, 1, 0}, {0, 1, 4, 1}, {0, 0, 1, 4}}, r4[4] = {5, 6, 6, 5}, x4[4];
    CHECK(solveSmall(4, Ms, r4, x4) && near(x4[0], 1) && near(x4[1], 1) && near(x4[2], 1) && near(x4[3], 1), "solveSmall");
}

void testPrimitives() {
    TriMesh c = primitives::cylinder(0.5f, 2.0f, 96);
    float v = c.signedVolume(), ve = kPi * 0.25f * 2.0f;
    CHECK(std::fabs(v - ve) / ve < 0.01f, "cylinder volume %f vs %f", v, ve);
    TriMesh s = primitives::sphere(1.0f, 96, 48);
    float vs = s.signedVolume(), vse = 4.0f / 3.0f * kPi;
    CHECK(std::fabs(vs - vse) / vse < 0.01f, "sphere volume %f vs %f", vs, vse);
    TriMesh b = primitives::box({0.5f, 1.0f, 1.5f});
    CHECK(std::fabs(b.signedVolume() - 6.0f) < 1e-4f, "box volume %f", b.signedVolume());
    TriMesh w = primitives::nacaWing("2412", 1.0f, 2.0f, 80);
    // NACA xx12 section area ~ 0.0822 c^2
    float vw = w.signedVolume();
    CHECK(std::fabs(vw - 0.0822f * 2.0f) / (0.0822f * 2.0f) < 0.03f, "wing volume %f", vw);
    TriMesh tb = primitives::streamlinedBody(1.0f, 0.3f);
    CHECK(tb.signedVolume() > 0, "streamlined volume %f", tb.signedVolume());
    // Welding a triangle soup restores shared vertices.
    TriMesh soup;
    for (auto t : b.triangles) {
        uint32_t base = uint32_t(soup.positions.size());
        for (int k = 0; k < 3; ++k) soup.positions.push_back(b.positions[t[k]]);
        soup.triangles.push_back({base, base + 1, base + 2});
    }
    soup.weld(1e-5f);
    CHECK(soup.positions.size() == 8, "welded vertex count %zu", soup.positions.size());
}

void testBVH() {
    TriMesh m = primitives::sphere(0.7f, 64, 32);
    Matrix3x3 R = Quaternion::fromEuler(0.3f, 0.2f, 0.1f).toMatrix3x3();
    m.transform(R, {1.0f, 0.6f, 0.8f}, {0.1f, -0.2f, 0.3f});
    MeshBVH bvh;
    bvh.build(m);
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> U(-1.2f, 1.2f);
    int mismatch = 0, signErr = 0;
    for (int n = 0; n < 2000; ++n) {
        Vector3 p(U(rng), U(rng), U(rng));
        ClosestHit h;
        bool ok = bvh.closestPoint(p, kInf, h);
        // brute force
        float best = kInf;
        for (size_t t = 0; t < m.triangles.size(); ++t) {
            const auto& tri = m.triangles[t];
            // sample-free exact distance through a tiny BVH-free routine: reuse MeshBVH on single tri is
            // overkill; use a dense barycentric search instead.
            Vector3 a = m.positions[tri[0]], b = m.positions[tri[1]], c = m.positions[tri[2]];
            for (int i = 0; i <= 8; ++i)
                for (int j = 0; i + j <= 8; ++j) {
                    Vector3 q = a + (b - a) * (i / 8.0f) + (c - a) * (j / 8.0f);
                    best = std::min(best, length(p - q));
                }
        }
        if (!ok || h.distance > best + 1e-5f || h.distance < best - 0.02f) ++mismatch;
        // inside test vs ellipsoid equation (approximate near the surface)
        Vector3 l = R.transposed() * (p - Vector3(0.1f, -0.2f, 0.3f));
        float e = sqr(l.x / 0.7f) + sqr(l.y / (0.7f * 0.6f)) + sqr(l.z / (0.7f * 0.8f));
        if (std::fabs(e - 1.0f) > 0.05f && ((e < 1.0f) != (h.signedDistance < 0))) ++signErr;
    }
    CHECK(mismatch == 0, "closest point mismatches: %d", mismatch);
    CHECK(signErr == 0, "inside/outside errors: %d", signErr);

    RayHit hit;
    TriMesh sp = primitives::sphere(1.0f, 128, 64);
    MeshBVH sb;
    sb.build(sp);
    bool ok = sb.raycast({-5, 0.01f, 0.02f}, {1, 0, 0}, 100, hit);
    CHECK(ok && std::fabs(hit.t - 4.0f) < 0.01f, "ray t=%f", hit.t);
}

void testAABBTree() {
    // Random boxes inserted, moved and removed: after every operation the tree must be valid, its
    // queries must return exactly the objects a brute-force scan finds, and it must stay balanced.
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> U(-10, 10), S(0.1f, 1.0f), M(-0.3f, 0.3f), J(-0.01f, 0.01f);
    AABBTree tree(0.1f);
    std::vector<AABB> boxes;
    std::vector<int> proxy;
    auto randomBox = [&]() {
        Vector3 c(U(rng), U(rng), U(rng)), h(S(rng), S(rng), S(rng));
        return AABB(c - h, c + h);
    };
    bool valid = true, queriesOk = true;
    for (int i = 0; i < 1000; ++i) {
        boxes.push_back(randomBox());
        proxy.push_back(tree.insert(boxes.back(), i));
    }
    valid &= tree.validate();
    const int heightFull = tree.height();
    std::vector<char> alive(boxes.size(), 1);
    int reinserted = 0;
    for (int step = 0; step < 50; ++step) {
        for (size_t i = 0; i < boxes.size(); ++i) { // everything drifts a little
            if (!alive[i]) continue;
            Vector3 d(M(rng), M(rng), M(rng));
            boxes[i].lo += d;
            boxes[i].hi += d;
            reinserted += tree.update(proxy[i], boxes[i]);
        }
        for (int k = 0; k < 10; ++k) { // some leave, some come back
            size_t i = size_t(rng() % boxes.size());
            if (alive[i]) tree.remove(proxy[i]), alive[i] = 0;
            else boxes[i] = randomBox(), proxy[i] = tree.insert(boxes[i], int(i)), alive[i] = 1;
        }
        valid &= tree.validate();
        // Query: every truly overlapping object is reported; every reported one has an
        // overlapping fat box.
        AABB q = randomBox();
        q.hi += Vector3(2.0f);
        std::vector<int> found;
        tree.query(q, [&](int o) { found.push_back(o); });
        std::sort(found.begin(), found.end());
        for (size_t i = 0; i < boxes.size(); ++i) {
            if (!alive[i]) continue;
            bool reported = std::binary_search(found.begin(), found.end(), int(i));
            if (boxes[i].overlaps(q) && !reported) queriesOk = false;
            if (reported && !tree.fatBox(proxy[i]).overlaps(q)) queriesOk = false;
        }
    }
    // Ray: the leaf of a box straight ahead must be reported.
    AABBTree small(0.0f);
    small.insert(AABB({4, -1, -1}, {6, 1, 1}), 7);
    small.insert(AABB({4, 5, -1}, {6, 7, 1}), 8);
    std::vector<int> hits;
    small.raycast({0, 0, 0}, {1, 0, 0}, 100.0f, [&](int o) { hits.push_back(o); });
    const int n = tree.objectCount();
    const float logN = std::log2(float(n));
    // Resting objects jiggle in place (+-1 cm): their fat boxes must absorb it, no re-insertion.
    int jiggleReinserts = 0;
    for (int step = 0; step < 50; ++step)
        for (size_t i = 0; i < boxes.size(); ++i) {
            if (!alive[i]) continue;
            Vector3 d(J(rng), J(rng), J(rng));
            AABB b = boxes[i];
            b.lo += d;
            b.hi += d;
            jiggleReinserts += tree.update(proxy[i], b);
        }
    std::printf("  AABB tree: %d objects, height %d (full tree %d, log2 n = %.1f), re-insertions: %d drifting, %d jiggling in place\n", n,
                tree.height(), heightFull, logN, reinserted, jiggleReinserts);
    CHECK(valid, "tree invariants broken");
    CHECK(queriesOk, "tree query differs from brute force");
    CHECK(tree.height() <= int(2.0f * logN) + 2, "tree not balanced: height %d for %d objects", tree.height(), n);
    CHECK(hits.size() == 1 && hits[0] == 7, "raycast");
    CHECK(jiggleReinserts < n / 10, "jiggling objects must not touch the tree (%d re-insertions)", jiggleReinserts);
}

void testMassProperties() {
    // Box as a polyhedron: exact volume and inertia; the eigen solver must diagonalise a rotated tensor.
    TriMesh b = primitives::box({0.5f, 1.0f, 1.5f});
    MassProperties mp = computeMassProperties(b);
    CHECK(std::fabs(mp.volume - 6.0f) < 1e-4f, "volume %f", mp.volume);
    float m = 6.0f; // unit density
    float Ixx = m / 12 * (4 + 9), Iyy = m / 12 * (1 + 9), Izz = m / 12 * (1 + 4);
    CHECK(std::fabs(mp.inertia.m[0][0] - Ixx) < 1e-3f && std::fabs(mp.inertia.m[1][1] - Iyy) < 1e-3f &&
              std::fabs(mp.inertia.m[2][2] - Izz) < 1e-3f,
          "inertia %f %f %f", mp.inertia.m[0][0], mp.inertia.m[1][1], mp.inertia.m[2][2]);
    Matrix3x3 R = Quaternion::fromEuler(0.4f, 0.7f, -0.3f).toMatrix3x3();
    TriMesh rb = b;
    rb.transform(R, Vector3(1.0f), {0.3f, -0.2f, 0.1f});
    MassProperties mr = computeMassProperties(rb);
    Vector3 eig;
    Matrix3x3 V;
    symmetricEigen(mr.inertia, eig, V);
    float sorted[3] = {eig.x, eig.y, eig.z};
    std::sort(sorted, sorted + 3);
    CHECK(std::fabs(sorted[0] - Izz) < 1e-2f && std::fabs(sorted[2] - Ixx) < 1e-2f, "eigen %f %f %f", sorted[0], sorted[1], sorted[2]);
    CHECK(length(mr.com - Vector3(0.3f, -0.2f, 0.1f)) < 1e-4f, "com");
    Quaternion q = Quaternion::fromMatrix3x3(R);
    Matrix3x3 R2 = q.toMatrix3x3();
    float err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) err = std::max(err, std::fabs(R.m[i][j] - R2.m[i][j]));
    CHECK(err < 1e-5f, "quat <-> mat3 roundtrip %e", err);
}
