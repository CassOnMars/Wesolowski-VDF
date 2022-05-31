namespace VDF.Wesolowski
{
    using System.Linq;
    using System.Runtime.InteropServices;
    using System.Text;
    using MathGmp.Native;
    using Org.BouncyCastle.Crypto.Digests;

    // Example implementation of the Wesolowski VDF, implementing the IQF
    // classgroups math described in
    // https://raw.githubusercontent.com/Chia-Network/chiavdf/main/classgroups.pdf
    public class Program
    {
        public static void Main(string[] args)
        {
            var pt = ProofOfTime(
                Encoding.UTF8.GetBytes("hi bob"),
                10000000,
                40);
            var verify = VerifyProof(
                pt.Item1,
                pt.Item2,
                Encoding.UTF8.GetBytes("hi bob"),
                10000000,
                pt.Item3);
            Console.WriteLine(verify);
        }

        public static byte[] RandomBytesFromSeed(byte[] seed, uint count)
        {
            if (count > (32 * 65535))
            {
                throw new ArgumentException("count is too large", "count");
            }

            var blob = new List<byte>((int)count);
            ushort extra = 0;

            while (blob.Count < count)
            {
                var digest = new Sha3Digest();
                digest.BlockUpdate(seed, 0, seed.Length);

                var extraBits = new byte[]{
                    (byte)((extra & 0xFF00) >> 8),
                    (byte)(extra & 0xFF)
                };
                digest.BlockUpdate(seed, 0, extraBits.Length);

                var hash = new byte[digest.GetDigestSize()];
                digest.DoFinal(hash, 0);
                blob.AddRange(hash);
                extra += 1;
            }

            return blob.Take((int)count).ToArray();
        }

        public static mpz_t CreateDiscriminant(byte[] seed, short length)
        {
            mpz_t n = new mpz_t();
            ushort residue;

            byte extra = (byte)((byte)length & (byte)7);
            uint randomBytesLength = (uint)(((length + 7) >> 3) + 2);
            var randomBytes = RandomBytesFromSeed(seed, randomBytesLength);

            var nBytes = randomBytes.Take((int)randomBytesLength - 2).ToArray();
            var numeratorBytes = randomBytes.Skip((int)randomBytesLength - 2)
                .ToArray();
            var numerator = ((uint)numeratorBytes[0] << 8) +
                (uint)numeratorBytes[1];
            
            gmp_lib.mpz_init(n);
            size_t nBytesLength = new size_t((ulong)nBytes.Count());
            IntPtr nBytesPtr = Marshal.AllocHGlobal(nBytes.Count());
            Marshal.Copy(nBytes, 0, nBytesPtr, nBytes.Count());
            gmp_lib.mpz_import(
                n,
                nBytesLength,
                1,
                new size_t(1),
                0,
                new size_t(0),
                new void_ptr(nBytesPtr));
            gmp_lib.mpz_fdiv_q_2exp(
                n,
                n,
                new mp_bitcnt_t((uint)(8 - extra) & 7));
            Marshal.FreeHGlobal(nBytesPtr);
            
            residue = (ushort)VDF.Wesolowski.PrimesResidues.Residues[
                numerator % VDF.Wesolowski.PrimesResidues.Residues.Count()
            ];
            gmp_lib.mpz_setbit(n, new mp_bitcnt_t((uint)length - 1));

            var remainder = gmp_lib.mpz_fdiv_ui(
                n,
                VDF.Wesolowski.PrimesResidues.M);
            
            var value = new mpz_t();
            gmp_lib.mpz_init(value);

            long rawValue = 0;

            if (residue > remainder)
            {
                rawValue = (long)residue - (long)remainder;
            }
            else
            {
                rawValue = -((long)remainder - (long)residue);
            }

            char_ptr cp = new char_ptr(rawValue.ToString());
            gmp_lib.mpz_set_str(value, cp, 10);
            gmp_lib.mpz_add(n, n, value);

            while (true)
            {
                var sieve = new List<bool>(65536);
                for (var i = 0; i < 65536; i++) { sieve.Add(true); }

                foreach (var si in VDF.Wesolowski.PrimesResidues.SieveInfo)
                {
                    uint i = (gmp_lib.mpz_cdiv_ui(n, si.Item1) * si.Item2)
                        % si.Item1;

                    while (i < sieve.Count())
                    {
                        sieve[(int)i] = false;
                        i += si.Item1;
                    }
                }

                for (var i = 0; i < sieve.Count(); i++)
                {
                    if (sieve[i])
                    {
                        var q = VDF.Wesolowski.PrimesResidues.M * i;
                        var qValue = new mpz_t();
                        gmp_lib.mpz_init(qValue);
                        cp = new char_ptr(q.ToString());
                        gmp_lib.mpz_set_str(qValue, cp, 10);
                        gmp_lib.mpz_add(n, n, qValue);

                        if (gmp_lib.mpz_probab_prime_p(n, 2) > 0)
                        {
                            gmp_lib.mpz_neg(n, n);
                            return n;
                        }

                        gmp_lib.mpz_sub(n, n, qValue);
                    }
                }

                var mValue = new mpz_t();
                gmp_lib.mpz_init(mValue);
                cp = new char_ptr(
                    (VDF.Wesolowski.PrimesResidues.M * 65536).ToString());
                gmp_lib.mpz_set_str(mValue, cp, 10);
                gmp_lib.mpz_add(n, n, mValue);
            }
        }

        public static (ClassGroup, ClassGroup, ushort) ProofOfTime(
            byte[] challenge,
            ulong iterations,
            ushort bits)
        {
            if (iterations >= ((ulong)1) << 53)
            {
                throw new ArgumentException(
                    "iterations is too large",
                    "iterations");
            }
            
            var d = CreateDiscriminant(challenge, (short)bits);
            
            // d = b^2 - 4ac
            // b = 1
            // a = 2
            // d = 1 - 8c
            // c = (1 - d)/8
            
            char_ptr cpA = new char_ptr("2");
            char_ptr cpFourA = new char_ptr("8");
            char_ptr cpB = new char_ptr("1");

            var a = new mpz_t();
            var fourA = new mpz_t();
            var b = new mpz_t();
            var bSquared = new mpz_t();
            var c =  new mpz_t();
            gmp_lib.mpz_init(a);
            gmp_lib.mpz_init(fourA);
            gmp_lib.mpz_init(b);
            gmp_lib.mpz_init(bSquared);
            gmp_lib.mpz_init(c);
            gmp_lib.mpz_set_str(a, cpA, 10);
            gmp_lib.mpz_set_str(fourA, cpFourA, 10);
            gmp_lib.mpz_set_str(b, cpB, 10);
            gmp_lib.mpz_set_str(bSquared, cpB, 10);
            
            gmp_lib.mpz_sub(c, bSquared, d);
            gmp_lib.mpz_fdiv_q(c, c, fourA);

            var x = new ClassGroup();
            gmp_lib.mpz_set(x.a, a);
            gmp_lib.mpz_set(x.b, b);
            gmp_lib.mpz_set(x.c, c);
            gmp_lib.mpz_set(x.discriminant, d);

            var h = new ClassGroup();
            gmp_lib.mpz_set(h.a, a);
            gmp_lib.mpz_set(h.b, b);
            gmp_lib.mpz_set(h.c, c);
            gmp_lib.mpz_set(h.discriminant, d);

            var memory = Math.Log(10000000f, 2);
            var t = Math.Log((double)iterations, 2);
            var l = (ulong)(t - memory > 0f ?
                Math.Ceiling(Math.Pow(2, memory - 20)) :
                1);
            var intermediate = t * Math.Log(2) / (2f * l);
            var k = (byte)(Math.Max(
                Math.Round(
                    (Math.Log(intermediate) - Math.Log(Math.Log(intermediate)))
                    + 0.25f),
                1f));

            var q = l * k;
            var powers = new Dictionary<ulong, ClassGroup>();
            var powersChain = new List<ulong>();

            for (ulong i = 0; i < (ulong)(iterations / q); i++)
            {
                powersChain.Add((ulong)(i * q));
            }

            powersChain.Add(iterations);
            
            ulong previous = 0;
            foreach (var current in powersChain)
            {
                x.RepeatedSquare(current - previous);
                var state = new ClassGroup();
                gmp_lib.mpz_set(state.a, x.a);
                gmp_lib.mpz_set(state.b, x.b);
                gmp_lib.mpz_set(state.c, x.c);
                gmp_lib.mpz_set(state.discriminant, x.discriminant);
                powers.Add(current, state);
                previous = current;
            }


            return (x, GenerateProof(h, iterations, k, l, powers, bits), bits);
        }

        // Produces a hash of the two classgroups which is highly probable of
        // being a prime.
        public static mpz_t HashPrime(ClassGroup x, ClassGroup y)
        {
            var xa = new char_ptr("");
            var xb = new char_ptr("");
            var ya = new char_ptr("");
            var yb = new char_ptr("");
            
            gmp_lib.mpz_get_str(xa, 10, x.a);
            gmp_lib.mpz_get_str(xb, 10, x.b);
            gmp_lib.mpz_get_str(ya, 10, y.a);
            gmp_lib.mpz_get_str(yb, 10, y.b);

            var b = new mpz_t();
            gmp_lib.mpz_init(b);

            byte[] seed = Encoding.UTF8.GetBytes(
                xa.ToString() + ":" + xb.ToString() + ":" + ya.ToString() +
                ":" + yb.ToString());
            ulong i = 0;

            while (true)
            {
                var digest = new Sha3Digest();
                digest.BlockUpdate(Encoding.UTF8.GetBytes("prime"), 0, 5);
                var iBytes = BitConverter.GetBytes(i);
                digest.BlockUpdate(iBytes, 0, iBytes.Length);
                digest.BlockUpdate(seed, 0, seed.Length);

                var hash = new byte[digest.GetDigestSize()];
                digest.DoFinal(hash, 0);

                IntPtr hashPtr = Marshal.AllocHGlobal(hash.Count());
                Marshal.Copy(hash, 0, hashPtr, hash.Count());
                gmp_lib.mpz_import(
                    b,
                    new size_t((ulong)hash.Length),
                    1,
                    new size_t(1),
                    0,
                    new size_t(0),
                    new void_ptr(hashPtr));

                if (gmp_lib.mpz_probab_prime_p(b, 256) > 0) {
                    break;
                }

                i++;
            }

            return b;
        }

        public static ClassGroup GenerateProof(
            ClassGroup x,
            ulong iterations,
            byte k,
            ulong l,
            Dictionary<ulong, ClassGroup> powers,
            ushort bits)
        {
            var b = HashPrime(x, powers[iterations]);

            var k1 = k >> 1;
            var k0 = k - k1;
            var h = x.Identity();
            var ident = x.Identity();
            ulong kExponent = (ulong)1 << k;
            ulong k0Exponent = (ulong)1 << k0;
            ulong k1Exponent = (ulong)1 << k1;

            for (long j = (long)l - 1; j >= 0; j--)
            {
                h.Pow(kExponent);

                var yDict = new Dictionary<string, ClassGroup>();

                for (ulong jj = 0; jj < (ulong)(1 << k); jj++)
                {
                    var bMPZ = new mpz_t();
                    gmp_lib.mpz_init(bMPZ);
                    var bCP = new char_ptr(jj.ToString());
                    gmp_lib.mpz_set_str(bMPZ, bCP, 10);

                    var id = new ClassGroup();
                    gmp_lib.mpz_set(id.a, ident.a);
                    gmp_lib.mpz_set(id.b, ident.b);
                    gmp_lib.mpz_set(id.c, ident.c);
                    gmp_lib.mpz_set(id.discriminant, ident.discriminant);

                    var cp = new char_ptr("");
                    gmp_lib.mpz_get_str(cp, 10, bMPZ);
                    yDict.Add(cp.ToString(), id);
                }

                ulong end = (ulong)Math.Ceiling((double)iterations / (k * l));

                for (ulong jj = 0; jj < end; jj++)
                {
                    if (iterations < (k * (jj * l + (ulong)j + 1)))
                    {
                        continue;
                    }

                    var res = new mpz_t();
                    gmp_lib.mpz_init(res);
                    var two = new mpz_t();
                    gmp_lib.mpz_init(two);
                    var cp = new char_ptr("2");
                    gmp_lib.mpz_set_str(two, cp, 10);

                    var exp = new mpz_t();
                    gmp_lib.mpz_init(exp);

                    cp = new char_ptr(iterations.ToString());
                    var iterationsMPZ = new mpz_t();
                    gmp_lib.mpz_init(iterationsMPZ);
                    gmp_lib.mpz_set_str(iterationsMPZ, cp, 10);

                    cp = new char_ptr(k.ToString());
                    var kMPZ = new mpz_t();
                    gmp_lib.mpz_init(kMPZ);
                    gmp_lib.mpz_set_str(kMPZ, cp, 10);

                    cp = new char_ptr((jj * l + 1).ToString());
                    var iMPZ = new mpz_t();
                    gmp_lib.mpz_init(iMPZ);
                    gmp_lib.mpz_set_str(iMPZ, cp, 10);

                    gmp_lib.mpz_mul(exp, kMPZ, iMPZ);
                    gmp_lib.mpz_sub(exp, iterationsMPZ, exp);

                    cp = new char_ptr(((ulong)(1 << k)).ToString());
                    var bMod = new mpz_t();
                    gmp_lib.mpz_init(bMod);
                    gmp_lib.mpz_set_str(bMod, cp, 10);

                    gmp_lib.mpz_powm(res, two, exp, b);
                    
                    gmp_lib.mpz_mul(res, res, bMod);
                    gmp_lib.mpz_fdiv_q(res, res, b);
                    gmp_lib.mpz_get_str(cp, 10, res);
                    
                    yDict[cp.ToString()].Mul(
                        powers[(jj * k * l)],
                        new Context());
                }

                for (ulong b1 = 0; b1 < k1Exponent; b1++)
                {
                    var z = new ClassGroup();
                    gmp_lib.mpz_set(z.a, ident.a);
                    gmp_lib.mpz_set(z.b, ident.b);
                    gmp_lib.mpz_set(z.c, ident.c);
                    gmp_lib.mpz_set(z.discriminant, ident.discriminant);

                    for (ulong b0 = 0; b0 < k0Exponent; b0++)
                    {
                        z.Mul(
                            yDict[(b1 * k0Exponent + b0).ToString()],
                            new Context());
                    }

                    z.Pow(b1 * k0Exponent);
                    h.Mul(z, new Context());
                }

                for (ulong b0 = 0; b0 < k0Exponent; b0++)
                {
                    var z = new ClassGroup();
                    gmp_lib.mpz_set(z.a, ident.a);
                    gmp_lib.mpz_set(z.b, ident.b);
                    gmp_lib.mpz_set(z.c, ident.c);
                    gmp_lib.mpz_set(z.discriminant, ident.discriminant);

                    for (ulong b1 = 0; b1 < k1Exponent; b1++)
                    {                    
                        z.Mul(
                            yDict[(b1 * k0Exponent + b0).ToString()],
                            new Context());
                    }

                    z.Pow(b0);
                    h.Mul(z, new Context());
                }
            }

            return h;
        }

        public static bool VerifyProof(
            ClassGroup y,
            ClassGroup proof,
            byte[] challenge,
            ulong iterations,
            ushort bits)
        {
            if (iterations >= ((ulong)1) << 53)
            {
                throw new ArgumentException(
                    "iterations is too large",
                    "iterations");
            }
            
            var d = CreateDiscriminant(challenge, (short)bits);

            char_ptr cpA = new char_ptr("2");
            char_ptr cpFourA = new char_ptr("8");
            char_ptr cpB = new char_ptr("1");

            var a = new mpz_t();
            var fourA = new mpz_t();
            var b = new mpz_t();
            var bSquared = new mpz_t();
            var c =  new mpz_t();
            gmp_lib.mpz_init(a);
            gmp_lib.mpz_init(fourA);
            gmp_lib.mpz_init(b);
            gmp_lib.mpz_init(bSquared);
            gmp_lib.mpz_init(c);
            gmp_lib.mpz_set_str(a, cpA, 10);
            gmp_lib.mpz_set_str(fourA, cpFourA, 10);
            gmp_lib.mpz_set_str(b, cpB, 10);
            gmp_lib.mpz_set_str(bSquared, cpB, 10);
            
            gmp_lib.mpz_sub(c, bSquared, d);
            gmp_lib.mpz_fdiv_q(c, c, fourA);
            var x = new ClassGroup();
            x.a = a;
            x.b = b;
            x.c = c;
            x.discriminant = d;

            var prime = HashPrime(x, y);
            var res = new mpz_t();
            gmp_lib.mpz_init(res);
            var two = new mpz_t();
            gmp_lib.mpz_init(two);
            var cp = new char_ptr("2");
            gmp_lib.mpz_set_str(two, cp, 10);
            cp = new char_ptr(iterations.ToString());
            var iterationsMPZ = new mpz_t();
            gmp_lib.mpz_init(iterationsMPZ);
            gmp_lib.mpz_set_str(iterationsMPZ, cp, 10);

            gmp_lib.mpz_powm(res, two, iterationsMPZ, prime);
            proof.Pow(prime);
            x.Pow(res);
            proof.Mul(x, new Context());

            return gmp_lib.mpz_cmp(proof.a, y.a) == 0 &&
                gmp_lib.mpz_cmp(proof.b, y.b) == 0 &&
                gmp_lib.mpz_cmp(proof.c, y.c) == 0 &&
                gmp_lib.mpz_cmp(proof.discriminant, y.discriminant) == 0;
        }
    }
}