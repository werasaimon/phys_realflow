#pragma once
// Non-convex models of the samples (teapot, bunny), decomposed once into convex parts and shared
// by every body that uses them.
#include "rigid/Shapes.h"

#include <memory>

namespace rf {

std::shared_ptr<const CompoundShape> teapotShape();
std::shared_ptr<const CompoundShape> bunnyShape();

} // namespace rf
