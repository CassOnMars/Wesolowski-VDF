# Wesolowski VDF .NET Library

An implementation of the Wesolowski VDF. Heavily borrows formulae from Chia
project:

https://raw.githubusercontent.com/Chia-Network/chiavdf/main/classgroups.pdf

Note: this does not produce proofs compatible with Chia Network, as the hash
function chosen differs (SHA 3 vs SHA-256) as well as the serialization of the
classgroups when fed into the hash function also differs (see HashPrime).

This is not fully hardened for production use (i.e. does not check all places
for division by zero, etc.) and may break, be susceptible to timing attacks,
may cause your printer to halt and catch fire, etc. Use at your own risk.

## Building

Due to [license issues](https://github.com/NethermindEth/Math.Gmp.Native/pull/6)
with the latest GMP .NET port, I have forked the repository and built this
with a relative path reference to maintain license purity. You could however
swap this out with the normal Nethermind fork from NuGet, however be aware if
you use license compatibility analysis tools the buck will stop at this library
for providing the correct type of license for conflict purposes.

## License

Given its use of the GNU Multiple Precision library, it must be licensed in a
way compatible with the GMP licensing, and is thus offered under the same
dual licensing, as either GPLv2 or LGPLv3:

The Wesolowski VDF .NET library is free software; you can redistribute it and/or
modify it under the terms of either:

  * the GNU Lesser General Public License as published by the Free Software
    Foundation; either version 3 of the License, or (at your option) any later
    version.

or

  * the GNU General Public License as published by the Free Software Foundation;
    either version 2 of the License, or (at your option) any later version.

or both in parallel, as here.

The Wesolowski VDF .NET library is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details.

You should have received copies of the GNU General Public License and the GNU
Lesser General Public License along with the GNU MP Library.  If not, see
https://www.gnu.org/licenses/.