namespace OrbitalResonance
{
    public static class RandomExtensions  // https://bitbucket.org/Superbest/superbest-random/src/master/Superbest%20random/RandomExtensions.cs
    {
        public static double NextGaussian(this Random r, double mu = 0, double sigma = 1)
        {
            var u1 = r.NextDouble();
            var u2 = r.NextDouble();

            var rand_std_normal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                                  Math.Sin(2.0 * Math.PI * u2);

            var rand_normal = mu + sigma * rand_std_normal;

            return rand_normal;
        }

        public static double NextDoubleRange(this Random r, double min, double max)
        {
            return (max - min) * r.NextDouble() + min;
        }

        public static Keplerian NextKeplerian(this Random r, Keplerian mu, Keplerian sigma, double starMass)
        {
            double semiMajorAxis = Math.Abs(r.NextGaussian(mu.semiMajorAxis, sigma.semiMajorAxis + mu.semiMajorAxis / 2));
            double eccentricity = Math.Abs(r.NextGaussian(mu.eccentricity, sigma.eccentricity));
            double inclination = r.NextDoubleRange(-Math.PI/2, Math.PI/2);
            double longitudeAscending = r.NextDoubleRange(0, 2 * Math.PI);
            double argumentPeriapsis = r.NextDoubleRange(0, 2 * Math.PI);
            double trueAnomaly = r.NextDoubleRange(0, 2 * Math.PI);

            return new(starMass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, trueAnomaly);
        }
    }
}
