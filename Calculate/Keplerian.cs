namespace OrbitalResonance
{
    public class Keplerian
    {
        public double starMass;
        public double semiMajorAxis;
        public double eccentricity;
        public double inclination;
        public double longitudeAscending;
        public double argumentPeriapsis;
        public double trueAnomaly;

        public Keplerian(double starMass, double semiMajorAxis, double eccentricity, double inclination,
                         double longitudeAscending, double argumentPeriapsis, double trueAnomaly)
        {
            this.starMass = starMass;
            this.semiMajorAxis = semiMajorAxis;
            this.eccentricity = eccentricity;
            this.inclination = inclination;
            this.longitudeAscending = longitudeAscending;
            this.argumentPeriapsis = argumentPeriapsis;
            this.trueAnomaly = trueAnomaly;
        }

        // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
        public Body3D ToCartesian()
        {
            double mu = Constants.G * starMass;
        
            double eccentricAnormaly = Math.Acos((eccentricity + Math.Cos(trueAnomaly)) / (1 + eccentricity * Math.Cos(trueAnomaly)));  // https://en.wikipedia.org/wiki/Eccentric_anomaly
            double distance = semiMajorAxis * (1 - eccentricity * Math.Cos(eccentricAnormaly));
        
            Vector3D oPos = distance * new Vector3D(Math.Cos(trueAnomaly), Math.Sin(trueAnomaly), 0);
            Vector3D oVel = Math.Sqrt(mu * semiMajorAxis) / distance * new Vector3D(-Math.Sin(eccentricAnormaly),
                                                                                     Math.Sqrt(1 - eccentricity * eccentricity) * Math.Cos(eccentricAnormaly),
                                                                                     0);
        
            Vector3D pos = new(oPos.x * (Math.Cos(argumentPeriapsis) * Math.Cos(longitudeAscending) - Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)) - oPos.y * (Math.Sin(argumentPeriapsis) * Math.Cos(longitudeAscending) + Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)),
                               oPos.x * (Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending) + Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending)) + oPos.y * (Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending) - Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending)),
                               oPos.x * Math.Sin(argumentPeriapsis) * Math.Sin(inclination) + oPos.y * Math.Cos(argumentPeriapsis) * Math.Sin(inclination));
            Vector3D vel = new(oVel.x * (Math.Cos(argumentPeriapsis) * Math.Cos(longitudeAscending) - Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)) - oVel.y * (Math.Sin(argumentPeriapsis) * Math.Cos(longitudeAscending) + Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)),
                               oVel.x * (Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending) + Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending)) + oVel.y * (Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending) - Math.Sin(argumentPeriapsis) * Math.Sin(longitudeAscending)),
                               oVel.x * Math.Sin(argumentPeriapsis) * Math.Sin(inclination) + oVel.y * Math.Cos(argumentPeriapsis) * Math.Sin(inclination));

            return new(pos, vel);
        }

        public static Keplerian FromCartesian(double starMass, Vector3D pos, Vector3D vel)
        {
            double mu = Constants.G * starMass;
            Vector3D angMomentum = Vector3D.Cross(pos, vel);
            Vector3D nodeVector = Vector3D.Cross(Vector3D.ZUnit, vel);

            Vector3D eccentricityVector = ((vel.Magnitude * vel.Magnitude - mu / pos.Magnitude) * pos - Vector3D.Dot(pos, vel) * vel) / mu;
            double eccentricity = eccentricityVector.Magnitude;
            double energy = vel.Magnitude * vel.Magnitude / 2 - mu / pos.Magnitude;

            double semiMajorAxis, p;
            if (Math.Abs(eccentricity - 1.0) > double.Epsilon)
            {
                semiMajorAxis = -mu / (2 * energy);
                p = semiMajorAxis * (1 - eccentricity * eccentricity);
            }
            else
            {
                p = angMomentum.Magnitude * angMomentum.Magnitude / mu;
                semiMajorAxis = double.PositiveInfinity;
            }

            double inclination = Math.Acos(angMomentum.z / angMomentum.Magnitude);

            double longitudeAscending = Math.Acos(nodeVector.x / nodeVector.Magnitude);
            if (nodeVector.y < 0) longitudeAscending = 2 * Math.PI - longitudeAscending;

            double argumentPeriapsis = Math.Acos(Vector3D.Dot(nodeVector, eccentricityVector) / (nodeVector.Magnitude * eccentricity));
            if (eccentricityVector.z < 0) argumentPeriapsis = 2 * Math.PI - argumentPeriapsis;

            double trueAnormaly = Math.Acos(Vector3D.Dot(eccentricityVector, pos) / (eccentricity * pos.Magnitude));
            if (Vector3D.Dot(pos, vel) < 0) trueAnormaly = 2 * Math.PI - trueAnormaly;

            return new(starMass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, trueAnormaly);
        }

        public static Keplerian BasicFromCartesian(double starMass, Vector3D pos, Vector3D vel)
        {
            double mu = Constants.G * starMass;
            Vector3D angMomentum = Vector3D.Cross(pos, vel);

            double eccentricity = (((vel.Magnitude * vel.Magnitude - mu / pos.Magnitude) * pos - Vector3D.Dot(pos, vel) * vel) / mu).Magnitude;

            double semiMajorAxis;
            if (Math.Abs(eccentricity - 1.0) > double.Epsilon)
                semiMajorAxis = -mu / (vel.Magnitude * vel.Magnitude - 2 * mu / pos.Magnitude);
            else
                semiMajorAxis = double.PositiveInfinity;

            double inclination = Math.Acos(angMomentum.z / angMomentum.Magnitude);

            return new(starMass, semiMajorAxis, eccentricity, inclination, 0, 0, 0);
        }

        public static Keplerian Zero(double starMass) { return new(starMass, 0, 0, 0, 0, 0, 0); }
    }
}
