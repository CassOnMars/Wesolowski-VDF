namespace VDF.Wesolowski
{
    using MathGmp.Native;
    
    public struct Congruence
    {
        public Congruence()
        {
            this.g = new mpz_t();
            this.d = new mpz_t();
            this.q = new mpz_t();
            this.r = new mpz_t();
            gmp_lib.mpz_init(this.g);
            gmp_lib.mpz_init(this.d);
            gmp_lib.mpz_init(this.q);
            gmp_lib.mpz_init(this.r);
        }

        public mpz_t g;

        public mpz_t d;

        public mpz_t q;

        public mpz_t r;

        public void SolveLinearCongruence(
            mpz_t mu,
            mpz_t? v,
            mpz_t a,
            mpz_t b,
            mpz_t m)
        {
            gmp_lib.mpz_gcdext(this.g, this.d, mu, a, m);
            gmp_lib.mpz_divexact(this.q, b, this.g);
            gmp_lib.mpz_mul(this.r, this.q, this.d);
            gmp_lib.mpz_tdiv_r(mu, this.r, m);
            
            if (v != null)
            {
                gmp_lib.mpz_divexact(v, m, this.g);
            }
        }
    }
}