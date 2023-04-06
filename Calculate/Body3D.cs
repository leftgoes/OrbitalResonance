namespace OrbitalResonance
{
    public class Body3D
    {
        public Vector3D pos;
        public Vector3D vel;

        public Body3D()
        {
            pos = Vector3D.Zero;
            vel = Vector3D.Zero;
        }

        public Body3D(Vector3D pos, Vector3D vel)
        {
            this.pos = pos;
            this.vel = vel;
        }
    }

    public class Star : Body3D
    {
        public double mass;

        public Star(double mass) : base() {
            this.mass = mass;
        }
    }

    public class NonAttracting : Body3D
    {
        public NonAttracting(Vector3D pos, Vector3D vel) : base(pos, vel) { }

        private double EccentricAnormalyNumeric(double meanArnormaly, double eccentricity, int iterations)
        {
            double eccentricAnormaly = meanArnormaly;
            for (int i = 0; i < iterations; i++)
            {
                eccentricAnormaly = eccentricAnormaly - (eccentricAnormaly - eccentricity * Math.Sin(eccentricAnormaly) - meanArnormaly) / (1 - eccentricity * Math.Cos(eccentricAnormaly));
            }
            return eccentricAnormaly;
        }

        public static NonAttracting FromKeplerian(Star star, Keplerian keplerian)
        {
            var cartesian = keplerian.ToCartesian();
            return new(star.pos + cartesian.pos, star.vel + cartesian.vel);
        }

        // https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
        public Keplerian ToKeplerian(Star star)
        {
            double mu = Constants.G * star.mass;
            Vector3D posRelative = pos - star.pos;
            Vector3D velRelative = vel - star.vel;

            Vector3D angMomentum = Vector3D.Cross(posRelative, velRelative);
            Vector3D nodeVector = Vector3D.Cross(new Vector3D(0, 0, 1), velRelative);

            Vector3D eccentricityVector = ((velRelative.Magnitude * velRelative.Magnitude - mu / posRelative.Magnitude) * posRelative - Vector3D.Dot(posRelative, velRelative) * velRelative) / mu;
            double eccentricity = eccentricityVector.Magnitude;
            double energy = velRelative.Magnitude * velRelative.Magnitude / 2 - mu / posRelative.Magnitude;

            double semiMajorAxis, p;
            if (Math.Abs(eccentricity - 1.0) > double.Epsilon)
            {
                semiMajorAxis = -mu / (2 * energy);
                p = semiMajorAxis * (1 - eccentricity * eccentricity);
            } else {
                p = angMomentum.Magnitude * angMomentum.Magnitude / mu;
                semiMajorAxis = double.PositiveInfinity;
            }

            double inclination = Math.Acos(angMomentum.z / angMomentum.Magnitude);

            double longitudeAscending = Math.Acos(nodeVector.x / nodeVector.Magnitude);
            if (nodeVector.y < 0) longitudeAscending = 2 * Math.PI - longitudeAscending;

            double argumentPeriapsis = Math.Acos(Vector3D.Dot(nodeVector, eccentricityVector) / (nodeVector.Magnitude * eccentricity));
            if (eccentricityVector.z < 0) argumentPeriapsis = 2 * Math.PI - argumentPeriapsis;

            double trueAnormaly = Math.Acos(Vector3D.Dot(eccentricityVector, posRelative) / (eccentricity * posRelative.Magnitude));
            if (Vector3D.Dot(pos, velRelative) < 0) trueAnormaly = 2 * Math.PI - trueAnormaly;

            return new(star.mass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, trueAnormaly);
        }
    }

    public class Attracting : NonAttracting
    {
        public double mass;

        public Attracting(double mass, Vector3D pos, Vector3D vel) : base(pos, vel)
        {
            this.mass = mass;
        }

        public static Attracting FromKeplerian(Keplerian keplerian,
                                               double mass)
        {
            var parameters = keplerian.ToCartesian();
            return new(mass, parameters.pos, parameters.vel);
        }
    }
}
