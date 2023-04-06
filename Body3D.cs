using Microsoft.Win32.SafeHandles;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Security;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using System.Transactions;


namespace OrbitalResonance
{
    class Body3D
    {
        public Vector3D pos;
        public Vector3D vel;
        public Vector3D acc;

        public Body3D()
        {
            pos = Vector3D.Zero;
            vel = Vector3D.Zero;
            acc = Vector3D.Zero;
        }

        public Body3D(Vector3D pos, Vector3D vel, Vector3D acc)
        {
            this.pos = pos;
            this.vel = vel;
            this.acc = acc;
        }
    }

    class Star : Body3D
    {
        public double mass;

        public Star(double mass) : base() {
            this.mass = mass;
        }
    }

    class Particle : Body3D
    {
        public Particle(Vector3D pos, Vector3D vel) : base(pos, vel, Vector3D.Zero) { }

        private double EccentricAnormalyNumeric(double meanArnormaly, double eccentricity, int iterations)
        {
            double eccentricAnormaly = meanArnormaly;
            for (int i = 0; i < iterations; i++)
            {
                eccentricAnormaly = eccentricAnormaly - (eccentricAnormaly - eccentricity * Math.Sin(eccentricAnormaly) - meanArnormaly) / (1 - eccentricity * Math.Cos(eccentricAnormaly));
            }
            return eccentricAnormaly;
        }

        // https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
        protected (Vector3D pos, Vector3D vel) _FromKeplerian(double starMass, double semiMajorAxis, double eccentricity, double inclination, double longitudeAscending, double argumentPeriapsis, double trueAnormaly) {
            double mu = Constants.G * starMass;

            double eccentricAnormaly = Math.Acos((eccentricity + Math.Cos(trueAnormaly)) / (1 + eccentricity * Math.Cos(trueAnormaly)));  // https://en.wikipedia.org/wiki/Eccentric_anomaly

            double distance = semiMajorAxis * (1 - eccentricity * Math.Cos(eccentricAnormaly));

            Vector3D oPos = distance * new Vector3D(Math.Cos(trueAnormaly), Math.Sin(trueAnormaly), 0);
            Vector3D oVel = Math.Sqrt(mu * semiMajorAxis) / distance * new Vector3D(-Math.Sin(eccentricAnormaly),
                                                                                     Math.Sqrt(1 - eccentricity * eccentricity) * Math.Cos(eccentricAnormaly),
                                                                                     0);

            Vector3D pos = new(oPos.x * (Math.Cos(argumentPeriapsis) * Math.Cos(longitudeAscending) - Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)) - oPos.y * (Math.Sin(argumentPeriapsis) * Math.Cos(longitudeAscending) + Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)),
                               oPos.x * (Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending) + Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending)) + oPos.y * (Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending) - Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending)),
                               oPos.x * Math.Sin(argumentPeriapsis) * Math.Sin(inclination) + oPos.y * Math.Cos(argumentPeriapsis) * Math.Sin(inclination));
            Vector3D vel = new(oVel.x * (Math.Cos(argumentPeriapsis) * Math.Cos(longitudeAscending) - Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)) - oVel.y * (Math.Sin(argumentPeriapsis) * Math.Cos(longitudeAscending) + Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Sin(longitudeAscending)),
                               oVel.x * (Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending) + Math.Sin(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending)) + oVel.y * (Math.Cos(argumentPeriapsis) * Math.Cos(inclination) * Math.Cos(longitudeAscending) - Math.Cos(argumentPeriapsis) * Math.Sin(longitudeAscending)),
                               oVel.x * Math.Sin(argumentPeriapsis) * Math.Sin(inclination) + oVel.y * Math.Cos(argumentPeriapsis) * Math.Sin(inclination));

            return (pos, vel);
        }

        public Particle FromKeplerian(double starMass, double semiMajorAxis, double eccentricity, double inclination,
                                      double longitudeAscending, double argumentPeriapsis, double meanArnormaly)
        {
            var parameters = _FromKeplerian(starMass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, meanArnormaly);
            return new(parameters.pos, parameters.vel);
        }

        // https://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
        public (double semiMajorAxis, double eccentricity, double inclination,
                double longitudeAscending, double argumentPeriapsis, double trueArnomaly) ToKeplerian(double starMass)
        {
            double mu = Constants.G * starMass;
            Vector3D angMomentum = Vector3D.Cross(pos, vel);
            Vector3D nodeVector = Vector3D.Cross(new Vector3D(0, 0, 1), vel);

            Vector3D eccentricityVector = ((vel.Magnitude * vel.Magnitude - mu / pos.Magnitude) * pos - Vector3D.Dot(pos, vel) * vel) / mu;
            double eccentricity = eccentricityVector.Magnitude;
            double energy = vel.Magnitude * vel.Magnitude / 2 - mu / pos.Magnitude;

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

            double argumentPeriapsis = Math.Acos(Vector3D.Dot(nodeVector, eccentricityVector) / (nodeVector.Magnitude * eccentricity);
            if (eccentricityVector.z < 0) argumentPeriapsis = 2 * Math.PI - argumentPeriapsis;

            double trueAnormaly = Math.Acos(Vector3D.Dot(eccentricityVector, pos) / (eccentricity * pos.Magnitude));
            if (Vector3D.Dot(pos, vel) < 0) trueAnormaly = 2 * Math.PI - trueAnormaly;

            return (semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, trueAnormaly);
        }
    }

    class Planet : Particle
    {
        public double mass;

        public Planet(double mass, Vector3D pos, Vector3D vel) : base(pos, vel)
        {
            this.mass = mass;
        }

        public Planet FromKeplerian(double starMass, double semiMajorAxis, double eccentricity, double inclination,
                                    double longitudeAscending, double argumentPeriapsis, double meanArnormaly,
                                    double mass)
        {
            var parameters = _FromKeplerian(starMass, semiMajorAxis, eccentricity, inclination, longitudeAscending, argumentPeriapsis, meanArnormaly);
            return new(mass, parameters.pos, parameters.vel);
        }
    }
}
