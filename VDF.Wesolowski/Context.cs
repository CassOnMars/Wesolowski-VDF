namespace VDF.Wesolowski
{
    using MathGmp.Native;
    
    public struct Context
    {
        public Context()
        {
            this.negative_a = new mpz_t();
            this.r = new mpz_t();
            this.denominator = new mpz_t();
            this.prev_a = new mpz_t();
            this.prev_b = new mpz_t();
            this.ra = new mpz_t();
            this.s = new mpz_t();
            this.x = new mpz_t();
            this.congruence = new Congruence();
            this.h = new mpz_t();
            this.w = new mpz_t();
            this.m = new mpz_t();
            this.u = new mpz_t();
            this.a = new mpz_t();
            this.l = new mpz_t();
            this.j = new mpz_t();
            this.b = new mpz_t();
            this.k = new mpz_t();
            this.t = new mpz_t();
            this.mu = new mpz_t();
            this.v = new mpz_t();
            this.sigma = new mpz_t();
            this.lambda = new mpz_t();
            gmp_lib.mpz_init(this.negative_a);
            gmp_lib.mpz_init(this.r);
            gmp_lib.mpz_init(this.denominator);
            gmp_lib.mpz_init(this.prev_a);
            gmp_lib.mpz_init(this.prev_b);
            gmp_lib.mpz_init(this.ra);
            gmp_lib.mpz_init(this.s);
            gmp_lib.mpz_init(this.x);
            gmp_lib.mpz_init(this.h);
            gmp_lib.mpz_init(this.w);
            gmp_lib.mpz_init(this.m);
            gmp_lib.mpz_init(this.u);
            gmp_lib.mpz_init(this.a);
            gmp_lib.mpz_init(this.l);
            gmp_lib.mpz_init(this.j);
            gmp_lib.mpz_init(this.b);
            gmp_lib.mpz_init(this.k);
            gmp_lib.mpz_init(this.t);
            gmp_lib.mpz_init(this.mu);
            gmp_lib.mpz_init(this.v);
            gmp_lib.mpz_init(this.sigma);
            gmp_lib.mpz_init(this.lambda);
        }

        public mpz_t negative_a;

        public mpz_t r;

        public mpz_t denominator;
        
        public mpz_t prev_a;
        
        public mpz_t prev_b;
        
        public mpz_t ra;
        
        public mpz_t s;
        
        public mpz_t x;
        
        public Congruence congruence;
        
        public mpz_t h;
        
        public mpz_t w;
        
        public mpz_t m;
        
        public mpz_t u;
        
        public mpz_t a;
        
        public mpz_t l;
        
        public mpz_t j;
        
        public mpz_t b;
        
        public mpz_t k;
        
        public mpz_t t;
        
        public mpz_t mu;
        
        public mpz_t v;
        
        public mpz_t sigma;
        
        public mpz_t lambda;
    }
}