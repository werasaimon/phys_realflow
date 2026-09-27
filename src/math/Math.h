#pragma once
// All math types of the project. Include this, or a single header for just one type.
//
//   Scalar.h      kPi, kInf, clampv, sqr, lerp, degToRad
//   ElementaryFunctions.h  rf::sin, cos, asin, acos, atan, atan2, exp, pow (one place for them)
//   Vector2.h     Vector2
//   Vector3.h     Vector3
//   Vector4.h     Vector4
//   Matrix3x3.h   Matrix3x3, symmetricEigen (3x3)
//   Quaternion.h  Quaternion, slerp
//   Matrix4x4.h   Matrix4x4 (affine transforms, projections)
//   MatrixNxN.h   MatrixNxN (LU, Cholesky, Jacobi eigen), solveSmall (n <= 4)
//   AABB.h        AABB

#include "math/AABB.h"
#include "math/ElementaryFunctions.h"
#include "math/Matrix3x3.h"
#include "math/Matrix4x4.h"
#include "math/MatrixNxN.h"
#include "math/Quaternion.h"
#include "math/Scalar.h"
#include "math/Vector2.h"
#include "math/Vector3.h"
#include "math/Vector4.h"
