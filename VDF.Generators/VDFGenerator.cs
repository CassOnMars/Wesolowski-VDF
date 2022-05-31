namespace VDF.Generators
{
    using System.Text;
    using Microsoft.CodeAnalysis;
    using Microsoft.CodeAnalysis.Text;

    [Generator]
    public class VDFGenerator : ISourceGenerator
    {
        public uint ModularExponentiation(uint @base, uint exponent, uint modulus)
        {
            if (@base >= ushort.MaxValue || exponent >= ushort.MaxValue || modulus >= ushort.MaxValue)
            {
                throw new Exception("base, exponent, or modulus is greater than 65536");
            }

            uint state = 1;
            
            while (true)
            {
                if ((exponent & 1) != 0)
                {
                    state *= @base;
                    state %= modulus;
                }

                exponent >>= 1;

                if (exponent == 0)
                {
                    return state;
                }

                @base *= @base;
                @base %= modulus;
            }
        }

        public void Execute(GeneratorExecutionContext context)
        {
            var sourceBuilder = new StringBuilder(@"
using System;

namespace VDF.Wesolowski
{
    public static class PrimesResidues
    {
        public static uint M = 8 * 3 * 5 * 7 * 11 * 13;

        ");

        uint m = 8 * 3 * 5 * 7 * 11 * 13;
        var primes = new List<uint> { 3, 5, 7, 11, 13 };
        var notDivisible = (uint x) => primes.All(p => x % p != 0);
        var residues = new List<uint>();
        
        for (uint x = 7; x < m; x += 8)
        {
            if (notDivisible(x))
            {
                residues.Add(x);
            }
        }

        if (residues.Count != 5760) { throw new Exception("invalid generation of residues"); }

        sourceBuilder.Append("public static int[] Residues = new int[] { ");
        sourceBuilder.Append(string.Join(", ", residues));
        sourceBuilder.Append(" };");

        primes = new List<uint>();
        uint N = 1 << 16;
        var sieve = new List<bool>((int)N >> 1);

        for (var i = 0; i < (N >> 1); i++){ sieve.Add(true); }

        var q = (uint)256;
        
        for (var i = 3; i <= q; i += 2)
        {
            if (sieve[i >> 1])
            {
                for (var j = (i * i >> 1); j < sieve.Count; j+=i)
                {
                    sieve[j] = false;
                }
            }
        }

        for (var i = 1; i < N / 2; i++)
        {
            if (sieve[i])
            {
                primes.Add((uint)(2 * i + 1));
            }
        }

        var sieveInfo = new List<(uint, uint)>();
        for (uint i = 5; i < primes.Count; i++)
        {
            sieveInfo.Add((i, ModularExponentiation(m % i, i - 2, i)));
        }

        sourceBuilder.Append(@"public static (ushort, ushort)[] SieveInfo = new (ushort, ushort)[] { (");
        sourceBuilder.Append(string.Join("), (", sieveInfo.Select(s => $"(ushort){s.Item1}, (ushort){s.Item2}")));
        sourceBuilder.Append(") };");

        sourceBuilder.Append(@"
    }
}
            ");

            context.AddSource("primesResidues", SourceText.From(sourceBuilder.ToString(), Encoding.UTF8));
        }

        public void Initialize(GeneratorInitializationContext context)
        {

        }
    }
}
