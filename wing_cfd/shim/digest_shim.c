/* Work-around for the Ubuntu 24.04 OpenFOAM v1912 package: Foam::OSHA1stream
 * goes bad on the first write, so Foam::dictionary::digest() aborts every run
 * that has function objects. The digest is only used to detect changed
 * function-object dictionaries on re-read, so a constant value is harmless
 * (function objects are simply kept as they are when controlDict is edited).
 * Build: gcc -O2 -shared -fPIC -o libdigestshim.so digest_shim.c
 * Use:   LD_PRELOAD=$PWD/libdigestshim.so simpleFoam ...                  */
#include <string.h>

struct SHA1Digest { unsigned char v[20]; };

struct SHA1Digest _ZNK4Foam10dictionary6digestEv(const void *self)
{
    struct SHA1Digest d;
    (void)self;
    memset(d.v, 0x5a, sizeof d.v);
    return d;
}
