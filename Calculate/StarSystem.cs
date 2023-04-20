using Newtonsoft.Json;
using System.Numerics;
using System.Xml.Linq;

namespace OrbitalResonance
{
    public class CartesianData
    {
        public int steps;
        public int planetsCount;
        public int particlesCount;
        public double[,,] planets;
        public double[,,] particles;

        public CartesianData(int steps, int planetsCount, int particlesCount) {
            this.steps = steps;
            this.planetsCount = planetsCount;
            this.particlesCount = particlesCount;
            planets = new double[steps, planetsCount, 3];
            particles = new double[steps, particlesCount, 3];
        }

        public void AddPlanet(int step, int index, Vector3D pos)
        {
            planets[step, index, 0] = pos.x;
            planets[step, index, 1] = pos.y;
            planets[step, index, 2] = pos.z;
        }

        public void AddParticle(int step, int index, Vector3D pos)
        {    
            particles[step, index, 0] = pos.x;
            particles[step, index, 1] = pos.y;
            particles[step, index, 2] = pos.z;
        }
    }

    public class KeplerianData : CartesianData
    {
        public KeplerianData(int steps, int planetsCount, int particlesCount) : base(steps, planetsCount, particlesCount)
        {
        }

        public void AddPlanet(int step, int index, Keplerian keplerian)
        {
            planets[step, index, 0] = keplerian.semiMajorAxis;
            planets[step, index, 1] = keplerian.eccentricity;
            planets[step, index, 2] = keplerian.inclination;
        }
        public void AddParticle(int step, int index, Keplerian keplerian)
        {
            particles[step, index, 0] = keplerian.semiMajorAxis;
            particles[step, index, 1] = keplerian.eccentricity;
            particles[step, index, 2] = keplerian.inclination;
        }
    }


    public class StarSystem
    {
        Star mainStar;
        Attracting[] planets;
        NonAttracting[] particles;

        public StarSystem(double starMass, Attracting[] planets, int particlesCount)
        {
            this.planets = planets;
            mainStar = new(starMass);
            particles = new NonAttracting[particlesCount];
        }

        public void AddParticles(int count)
        {
            Random random = new();

            var planetsDistribution = PlanetsKeplerianDistribution();

            particles = new NonAttracting[count];
            for (int i = 0; i < count; i++)
                particles[i] = NonAttracting.FromKeplerian(mainStar, random.NextKeplerian(planetsDistribution.mu, planetsDistribution.sigma, mainStar.mass));
        }

        public (Keplerian mu, Keplerian sigma) PlanetsKeplerianDistribution()
        {
            double length = planets.Length;
            Keplerian planetKeplerian;

            Keplerian mu = Keplerian.Zero(mainStar.mass);
            foreach (Attracting planet in planets)
            {
                planetKeplerian = planet.ToKeplerian(mainStar);

                mu.semiMajorAxis += planetKeplerian.semiMajorAxis;
                mu.eccentricity += planetKeplerian.eccentricity;
                mu.inclination += planetKeplerian.inclination;
                mu.longitudeAscending += planetKeplerian.longitudeAscending;
                mu.argumentPeriapsis += planetKeplerian.argumentPeriapsis;
                mu.trueAnomaly += planetKeplerian.trueAnomaly;
            }

            mu.semiMajorAxis /= length;
            mu.eccentricity /= length;
            mu.inclination /= length;
            mu.longitudeAscending /= length;
            mu.argumentPeriapsis /= length;
            mu.trueAnomaly /= length;

            Keplerian sigma = Keplerian.Zero(mainStar.mass);
            foreach (Attracting planet in planets)
            {
                planetKeplerian = planet.ToKeplerian(mainStar);

                sigma.semiMajorAxis += Math.Pow(planetKeplerian.semiMajorAxis - mu.semiMajorAxis, 2);
                sigma.eccentricity += Math.Pow(planetKeplerian.eccentricity - mu.eccentricity, 2);
                sigma.inclination += Math.Pow(planetKeplerian.inclination - mu.inclination, 2);
                sigma.longitudeAscending += Math.Pow(planetKeplerian.longitudeAscending - mu.longitudeAscending, 2);
                sigma.argumentPeriapsis += Math.Pow(planetKeplerian.argumentPeriapsis - mu.argumentPeriapsis, 2);
                sigma.trueAnomaly += Math.Pow(planetKeplerian.trueAnomaly - mu.trueAnomaly, 2);
            }

            sigma.semiMajorAxis = Math.Sqrt(sigma.semiMajorAxis);
            sigma.eccentricity = Math.Sqrt(sigma.eccentricity);
            sigma.inclination = Math.Sqrt(sigma.inclination);
            sigma.longitudeAscending = Math.Sqrt(sigma.longitudeAscending);
            sigma.argumentPeriapsis = Math.Sqrt(sigma.argumentPeriapsis);
            sigma.trueAnomaly = Math.Sqrt(sigma.trueAnomaly);

            return (mu, sigma);
        }

        private (Vector3D star, Vector3D[] planets) AccelerationStarAndPlanets()
        {
            Vector3D mainStarAcc = Vector3D.Zero;
            Vector3D[] planetsAcc = new Vector3D[planets.Length];

            for (int i = 0; i < planets.Length; i++)
            {
                Vector3D delta = planets[i].pos - mainStar.pos;
                Vector3D aoverm = Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                mainStarAcc += planets[i].mass * aoverm;
                planetsAcc[i] = -mainStar.mass * aoverm;
            }
            
            for (int i = 0; i < planets.Length; i++)
            {
                for (int j = i + 1; j < planets.Length; j++)
                {
                    Vector3D delta = planets[j].pos - planets[i].pos;
                    Vector3D aoverm = Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                    planetsAcc[i] += planets[j].mass * aoverm;
                    planetsAcc[j] -= planets[i].mass * aoverm;
                }
            }

            return (mainStarAcc, planetsAcc);
        }

        private Vector3D[] AccelerationParticles()
        {
            Vector3D[] accelerations = new Vector3D[particles.Length];

            for (int i = 0; i < particles.Length; i++)
            {
                Vector3D delta = mainStar.pos - particles[i].pos;
                accelerations[i] = mainStar.mass * Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                foreach (Attracting planet in planets)
                {
                    delta = planet.pos - particles[i].pos;
                    accelerations[i] += planet.mass * Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                }
            }

            return accelerations;
        }

        private void UpdatePositions(double dt, Vector3D mainStarAcc, Vector3D[] planetsAcc, Vector3D[] particlesAcc)
        {
            mainStar.pos += dt * (mainStar.vel + dt * mainStarAcc / 2);
            for (int i = 0; i < planets.Length; i++)
                planets[i].pos += dt * (planets[i].vel + dt * planetsAcc[i] / 2);
            for (int i = 0; i < particles.Length; i++)
                particles[i].pos += dt * (particles[i].vel + dt * particlesAcc[i] / 2);
        }

        private void UpdateVelocities(double dt, Vector3D mainStarAcc, Vector3D[] planetsAcc, Vector3D[] particlesAcc,
                                                 Vector3D newMainStarAcc, Vector3D[] newPlanetsAcc, Vector3D[] newParticlesAcc)
        {
            mainStar.vel += dt * (mainStarAcc + newMainStarAcc) / 2;
            for (int i = 0; i < planets.Length; i++)
                planets[i].vel += dt * (planetsAcc[i] + newPlanetsAcc[i]) / 2;
            for (int i = 0; i < particles.Length; i++)
                particles[i].vel += dt * (particlesAcc[i] + newParticlesAcc[i]) / 2;
        }

        public void NextStep(double dt)  // https://gamedev.stackexchange.com/questions/15708/how-can-i-implement-gravity
        {
            var (mainStarAcc, planetsAcc) = AccelerationStarAndPlanets();
            Vector3D[] particlesAcc = AccelerationParticles();

            UpdatePositions(dt, mainStarAcc, planetsAcc, particlesAcc);

            var (newMainStarAcc, newPlanetsAcc) = AccelerationStarAndPlanets();
            Vector3D[] newParticlesAcc = AccelerationParticles();

            UpdateVelocities(dt, mainStarAcc, planetsAcc, particlesAcc, newMainStarAcc, newPlanetsAcc, newParticlesAcc);

            mainStar.pos += mainStar.vel * dt;
            mainStar.vel += mainStarAcc * dt;
        }

        public void SimulateCartesian(string filename, int steps = 1000, double dt = 86400, int substeps = 1)
        {
            CartesianData cData = new(steps, planets.Length, particles.Length);
            for (int step = 0; step < steps; step++)
            {
                for (int i = 0; i < substeps; i++)
                    NextStep(dt);
                
                for (int i = 0; i < planets.Length; i++)
                    cData.AddPlanet(step, i, planets[i].pos - mainStar.pos);

                for (int i = 0; i < particles.Length; i++)
                    cData.AddParticle(step, i, particles[i].pos - mainStar.pos);
            }

            string jsonString = JsonConvert.SerializeObject(cData);
            File.WriteAllText(filename, jsonString);
        }

        public void SimulateKeplerian(string filename, int steps = 1000, double dt = 86400, int every = 1)
        {
            KeplerianData kData = new(steps / every, planets.Length, particles.Length);
            for (int step = 0; step < steps; step++)
            {
                NextStep(dt);
                if (step % every != 0) continue;

                for (int i = 0; i < planets.Length; i++)
                {
                    Keplerian keplerian = Keplerian.BasicFromCartesian(mainStar.mass, planets[i].pos - mainStar.pos, planets[i].vel - mainStar.vel);
                    kData.AddPlanet(step / every, i, keplerian);
                }

                for (int i = 0; i < particles.Length; i++)
                {
                    Keplerian keplerian = Keplerian.BasicFromCartesian(mainStar.mass, particles[i].pos - mainStar.pos, particles[i].vel - mainStar.vel);
                    kData.AddParticle(step / every, i, keplerian);
                }
            }

            string jsonString = JsonConvert.SerializeObject(kData);
            File.WriteAllText(filename, jsonString);
        }
    }
}