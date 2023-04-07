using Newtonsoft.Json;

namespace OrbitalResonance
{
    public class StarSystemCartesianData
    {
        public int steps;
        public int planetsCount;
        public double[,] star;
        public double[,,] planets;

        public StarSystemCartesianData(int steps, int planetsCount) {
            this.steps = steps;
            this.planetsCount = planetsCount;
            star = new double[steps, 3];
            planets = new double[steps, planetsCount, 3];
        }

        public void AddStar(int step, Vector3D pos)
        {
            star[step, 0] = pos.x;
            star[step, 1] = pos.y;
            star[step, 2] = pos.z;
        }

        public void AddPlanet(int step, int index, Vector3D pos)
        {
            planets[step, index, 0] = pos.x;
            planets[step, index, 1] = pos.y;
            planets[step, index, 2] = pos.z;
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

        public NonAttracting[] GetRandomParticles(int count)
        {
            Random random = new();

            var planetsDistribution = PlanetsKeplerianDistribution();

            particles = new NonAttracting[count];
            for (int i = 0; i < count; i++)
            {
                particles[i] = NonAttracting.FromKeplerian(mainStar, random.NextKeplerian(planetsDistribution.mu, planetsDistribution.sigma, mainStar.mass));
            }

            return particles;
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

        public void NextStep(double dt)
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

            for (int i = 0; i < planets.Length; i++) {
                for (int j = i + 1; j < planets.Length; j++)
                {
                    Vector3D delta = planets[j].pos - planets[i].pos;
                    Vector3D aoverm = Constants.G / Math.Pow(delta.Magnitude, 3) * delta;
                    planetsAcc[i] += planets[j].mass * aoverm;
                    planetsAcc[j] -= planets[i].mass * aoverm;
                }
            }

            mainStar.vel += mainStarAcc * dt;
            mainStar.pos += mainStar.vel * dt;

            for (int i = 0; i < planets.Length; i++)
            {
                planets[i].vel += planetsAcc[i] * dt;
                planets[i].pos += planets[i].vel * dt;
            }
        }

        public void SimulateCartesian(string filename, int steps = 1000, double dt = 86400, int every = 1)
        {
            StarSystemCartesianData cData = new(steps / every, planets.Length);
            for (int step = 0; step < steps; step++)
            {
                NextStep(dt);
                if (step % every != 0) continue;

                cData.AddStar(step / every, mainStar.pos);
                
                for (int i = 0; i < planets.Length; i++)
                    cData.AddPlanet(step / every, i, planets[i].pos - mainStar.pos);
            }

            string jsonString = JsonConvert.SerializeObject(cData);
            File.WriteAllText(filename, jsonString);
        }
    }
}