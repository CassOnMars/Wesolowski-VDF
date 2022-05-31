namespace VDF.Wesolowski
{
    using MathGmp.Native;
    
    public struct ClassGroup
    {
        public ClassGroup()
        {
            this.a = new mpz_t();
            this.b = new mpz_t();
            this.c = new mpz_t();
            this.discriminant = new mpz_t();
            gmp_lib.mpz_init(this.a);
            gmp_lib.mpz_init(this.b);
            gmp_lib.mpz_init(this.c);
            gmp_lib.mpz_init(this.discriminant);
        }

        public mpz_t a;
        
        public mpz_t b;
        
        public mpz_t c;
        
        public mpz_t discriminant;

        public (mpz_t, mpz_t) IntoRaw() => (this.a, this.b);

        public ClassGroup Identity() 
        {
            var identity = new ClassGroup();
            var fourcp = new char_ptr("4");
            var onecp = new char_ptr("1");
            gmp_lib.mpz_set_str(identity.a, fourcp, 10);
            gmp_lib.mpz_set_str(identity.b, onecp, 10);
            gmp_lib.mpz_set(identity.discriminant, this.discriminant);
            gmp_lib.mpz_sub(identity.c, identity.b, identity.discriminant);
            gmp_lib.mpz_fdiv_q(identity.c, identity.c, identity.a);
            gmp_lib.mpz_set_str(identity.a, onecp, 10);

            return identity;
        }

        public void Mul(ClassGroup rhs, Context ctx)
        {
            // g = (lhs.b + rhs.b) / 2
            gmp_lib.mpz_add(ctx.congruence.g, this.b, rhs.b);
            gmp_lib.mpz_fdiv_q_ui(ctx.congruence.g, ctx.congruence.g, 2);

            // h = (rhs.b - lhs.b) / 2
            gmp_lib.mpz_sub(ctx.h, rhs.b, this.b);
            gmp_lib.mpz_fdiv_q_ui(ctx.h, ctx.h, 2);

            // w = gcd(lhs.a, rhs.a, g)
            gmp_lib.mpz_gcd(ctx.w, this.a, rhs.a);
            gmp_lib.mpz_gcd(ctx.w, ctx.w, ctx.congruence.g);

            // j = w
            gmp_lib.mpz_set(ctx.j, ctx.w);

            // s = lhs.a / w
            gmp_lib.mpz_fdiv_q(ctx.s, this.a, ctx.w);
            
            // t = rhs.a / w
            gmp_lib.mpz_fdiv_q(ctx.t, rhs.a, ctx.w);

            // u = group / w
            gmp_lib.mpz_fdiv_q(ctx.u, ctx.congruence.g, ctx.w);

            // a = t * u
            gmp_lib.mpz_mul(ctx.a, ctx.t, ctx.u);

            // b = h * u + s * lhs.c  :think_thonk:
            gmp_lib.mpz_mul(ctx.b, ctx.h, ctx.u);
            gmp_lib.mpz_mul(ctx.m, ctx.s, this.c);
            gmp_lib.mpz_add(ctx.b, ctx.b, ctx.m);

            // m = s * t
            gmp_lib.mpz_mul(ctx.m, ctx.s, ctx.t);
            ctx.congruence
                .SolveLinearCongruence(ctx.mu, ctx.v, ctx.a, ctx.b, ctx.m);

            // a = t * v
            gmp_lib.mpz_mul(ctx.a, ctx.t, ctx.v);

            // b = h - t * mu
            gmp_lib.mpz_mul(ctx.m, ctx.t, ctx.mu);
            gmp_lib.mpz_sub(ctx.b, ctx.h, ctx.m);

            // m = s
            gmp_lib.mpz_set(ctx.m, ctx.s);

            ctx.congruence.SolveLinearCongruence(
                ctx.lambda,
                ctx.sigma,
                ctx.a,
                ctx.b,
                ctx.m);

            // k = mu + v * lambda
            gmp_lib.mpz_mul(ctx.a, ctx.v, ctx.lambda);
            gmp_lib.mpz_add(ctx.k, ctx.mu, ctx.a);

            // l = (k * t - h)/s
            gmp_lib.mpz_mul(ctx.l, ctx.k, ctx.t);
            gmp_lib.mpz_sub(ctx.v, ctx.l, ctx.h);
            gmp_lib.mpz_fdiv_q(ctx.l, ctx.v, ctx.s);

            // m = (tuk - hu - cs)/st
            gmp_lib.mpz_mul(ctx.m, ctx.t, ctx.u);
            gmp_lib.mpz_mul(ctx.m, ctx.m, ctx.k);
            gmp_lib.mpz_mul(ctx.a, ctx.h, ctx.u);
            gmp_lib.mpz_sub(ctx.m, ctx.m, ctx.a);
            gmp_lib.mpz_mul(ctx.a, this.c, ctx.s);
            gmp_lib.mpz_sub(ctx.m, ctx.m, ctx.a);
            gmp_lib.mpz_mul(ctx.a, ctx.s, ctx.t);
            gmp_lib.mpz_fdiv_q(ctx.lambda, ctx.m, ctx.a);

            // A = st - ru
            gmp_lib.mpz_mul(this.a, ctx.s, ctx.t);

            // B = (ju + mr) - (kt + ls)
            gmp_lib.mpz_mul(this.b, ctx.j, ctx.u);
            gmp_lib.mpz_mul(ctx.a, ctx.k, ctx.t);
            gmp_lib.mpz_sub(this.b, this.b, ctx.a);
            gmp_lib.mpz_mul(ctx.a, ctx.l, ctx.s);
            gmp_lib.mpz_sub(this.b, this.b, ctx.a);

            // C = kl - jm
            gmp_lib.mpz_mul(this.c, ctx.k, ctx.l);
            gmp_lib.mpz_mul(ctx.a, ctx.j, ctx.lambda);
            gmp_lib.mpz_sub(this.c, this.c, ctx.a);

            this.Reduce(ctx);
        }

        public void Reduce(Context ctx)
        {
            this.Normalize(ctx);

            while (this.b._mp_size.Value < 0 ?
                gmp_lib.mpz_cmp(this.a, this.c) >= 0 :
                gmp_lib.mpz_cmp(this.a, this.c) > 0)
            {
                gmp_lib.mpz_add(ctx.s, this.c, this.b);
                gmp_lib.mpz_add(ctx.x, this.c, this.c);
                gmp_lib.mpz_swap(this.b, ctx.prev_b);
                gmp_lib.mpz_fdiv_q(this.b, ctx.s, ctx.x);
                gmp_lib.mpz_swap(this.b, ctx.s);
                gmp_lib.mpz_swap(this.a, this.c);

                // x = 2sc
                gmp_lib.mpz_mul(this.b, ctx.s, this.a);
                gmp_lib.mpz_mul_2exp(ctx.x, this.b, 1);

                // b = x - prev_b
                gmp_lib.mpz_sub(this.b, ctx.x, ctx.prev_b);

                // x = bs
                gmp_lib.mpz_mul(ctx.x, ctx.prev_b, ctx.s);

                // s = c(s^2)
                gmp_lib.mpz_mul(ctx.prev_b, ctx.s, ctx.s);
                gmp_lib.mpz_mul(ctx.s, this.a, ctx.prev_b);

                // c = s - x
                gmp_lib.mpz_sub(ctx.prev_a, ctx.s, ctx.x);

                // c = c + a
                gmp_lib.mpz_add(this.c, this.c, ctx.prev_a);
            }

            this.Normalize(ctx);
        }

        public void Normalize(Context ctx)
        {
            gmp_lib.mpz_neg(ctx.negative_a, this.a);

            if (gmp_lib.mpz_cmp(this.b, ctx.negative_a) > 0 &&
                gmp_lib.mpz_cmp(this.b, this.a) <= 0)
            {
                return;
            }

            gmp_lib.mpz_sub(ctx.r, this.a, this.b);
            gmp_lib.mpz_mul_2exp(ctx.denominator, this.a, 1);
            gmp_lib.mpz_fdiv_q(ctx.negative_a, ctx.r, ctx.denominator);
            gmp_lib.mpz_swap(ctx.negative_a, ctx.r);
            gmp_lib.mpz_swap(ctx.prev_b, this.b);

            gmp_lib.mpz_mul(ctx.ra, ctx.r, this.a);
            gmp_lib.mpz_mul_2exp(ctx.negative_a, ctx.ra, 1);
            gmp_lib.mpz_add(this.b, ctx.prev_b, ctx.negative_a);

            gmp_lib.mpz_mul(ctx.negative_a, ctx.ra, ctx.r);
            gmp_lib.mpz_add(ctx.prev_a, this.c, ctx.negative_a);

            gmp_lib.mpz_mul(ctx.ra, ctx.r, ctx.prev_b);
            gmp_lib.mpz_add(this.c, ctx.prev_a, ctx.ra);
        }

        public void Square(Context ctx)
        {
            ctx.congruence.SolveLinearCongruence(
                ctx.mu,
                null,
                this.b,
                this.c,
                this.a);

            gmp_lib.mpz_mul(ctx.m, this.b, ctx.mu);
            gmp_lib.mpz_sub(ctx.m, ctx.m, this.c);
            gmp_lib.mpz_fdiv_q(ctx.m, ctx.m, this.a);

            gmp_lib.mpz_set(ctx.prev_a, this.a);
            gmp_lib.mpz_mul(this.a, ctx.prev_a, ctx.prev_a);

            gmp_lib.mpz_mul(ctx.a, ctx.mu, ctx.prev_a);
            gmp_lib.mpz_mul_2exp(ctx.a, ctx.a, 1);
            gmp_lib.mpz_sub(this.b, this.b, ctx.a);

            gmp_lib.mpz_mul(this.c, ctx.mu, ctx.mu);
            gmp_lib.mpz_sub(this.c, this.c, ctx.m);

            this.Reduce(ctx);
        }

        public void RepeatedSquare(ulong iterations)
        {
            var ctx = new Context();
            for (ulong i = 0; i < iterations; i++)
            {
                this.Square(ctx);
            }
        }

        public void Pow(mpz_t exponent)
        {
            var exp = new mpz_t();
            gmp_lib.mpz_init(exp);
            gmp_lib.mpz_set(exp, exponent);

            var state = this.Identity();
            while (true)
            {
                var odd = gmp_lib.mpz_tstbit(exp, new mp_bitcnt_t(0)) == 1;

                // exp >> 1
                gmp_lib.mpz_fdiv_q_2exp(exp, exp, 1);

                if (odd)
                {
                    state.Mul(this, new Context());
                }

                if (exp._mp_size.Value == 0)
                {
                    this = state;
                    break;
                }

                this.Square(new Context());
            }
        }

        public void Pow(ulong exponent)
        {
            var exp = new mpz_t();
            gmp_lib.mpz_init(exp);
            var cp = new char_ptr(exponent.ToString());
            gmp_lib.mpz_set_str(exp, cp, 10);

            this.Pow(exp);
        }
    }
}