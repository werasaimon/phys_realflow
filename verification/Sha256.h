#pragma once
// SHA-256 (FIPS 180-4, "Secure Hash Standard", NIST 2015), for the commitments of blinded
// references: the registry carries only the hash of "value|uncertainty|salt", so the reference
// cannot be read from the sources, yet anyone can check later that the sealed file was not
// changed. Plain C++: 64 rounds over 512-bit blocks, as in the standard's section 6.2.

#include <string>

namespace rf::verify {

// The hash of the bytes of `message` as 64 lowercase hex digits.
std::string sha256Hex(const std::string& message);

} // namespace rf::verify
