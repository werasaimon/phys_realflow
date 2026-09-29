#pragma once
// Non-convex models of the samples (teapot, bunny, ring, bowl, cup - and a wheel and a table built of
// boxes), decomposed once into convex parts and shared
// by every body that uses them.
#include "rigid/Shapes.h"

#include <memory>

namespace rf {

std::shared_ptr<const CompoundShape> teapotShape();
std::shared_ptr<const CompoundShape> bunnyShape();
// A steel ring of a chain: 10 cm across, its tube 24 mm thick.
std::shared_ptr<const CompoundShape> ringShape();
// A wheel of boxes: a hub, six spokes and a rim of 48 flat segments; 0.6 m across, axis z.
std::shared_ptr<const CompoundShape> wheelShape();
// A table of boxes: a top 0.6 x 0.04 x 0.4 m on four legs, 0.44 m tall.
std::shared_ptr<const CompoundShape> tableShape();
// A bowl: a half shell 1.2 m across, 4 cm thick, open at the top.
std::shared_ptr<const CompoundShape> bowlShape();
// A tapered cup 0.1 m tall, 8 to 12 cm across, its wall 4 mm: identical cups nest 20 mm apart.
std::shared_ptr<const CompoundShape> cupShape();

} // namespace rf
