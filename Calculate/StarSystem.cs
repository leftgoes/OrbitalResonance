using Calculate;
using Newtonsoft.Json;
using System.Diagnostics;

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
        bool nonescaping;

        public KeplerianData(int steps, int planetsCount, int particlesCount, bool nonescaping) : base(steps, planetsCount, particlesCount)
        {
            this.nonescaping = nonescaping;
        }

        public void AddPlanet(int step, int planetIndex, Keplerian keplerian)
        {
            planets[step, planetIndex, 0] = keplerian.semiMajorAxis;
            planets[step, planetIndex, 1] = keplerian.eccentricity;
            planets[step, planetIndex, 2] = keplerian.inclination;
        }
        public void AddParticle(int step, int particleIndex, Keplerian keplerian)
        {
            particles[step, particleIndex, 0] = keplerian.semiMajorAxis;
            particles[step, particleIndex, 1] = keplerian.eccentricity;
            particles[step, particleIndex, 2] = keplerian.inclination;
        }

        private double ParticlesPercentile(double percentile, int keplerianIndex)
        {
            int length = 0;

            for (int i = 0; i < steps; i++) {
                for (int j = 0; j < particlesCount; j++) {
                    if (nonescaping && particles[i, j, 1] > 1)
                        continue;
                    length++;
                }
            }

            double[] flattened = new double[steps * particlesCount];

            int flattenedIndex = 0;
            for (int i = 0; i < steps; i++) {
                for (int j = 0; j < particlesCount; j++) {
                    if (nonescaping && particles[i, j, 1] > 1)
                        continue;
                    flattened[flattenedIndex++] = particles[i, j, keplerianIndex];
                }
            }
            Array.Sort(flattened);

            double percentileIndex = percentile / 100 * length;

            int percentileIndexInt = (int)percentileIndex;
            double t = percentileIndex - percentileIndexInt;

            return (1 - t) * flattened[percentileIndexInt] + t * flattened[percentileIndexInt + 1];
        }

        public VideoArray ToVideoArray(int width, int height, int xIndex, int yIndex, double percentile = 5)
        {
            DoubleRange xRange = new(ParticlesPercentile(percentile, xIndex), ParticlesPercentile(100.0 - percentile, xIndex));
            DoubleRange yRange = new(ParticlesPercentile(percentile, yIndex), ParticlesPercentile(100.0 - percentile, yIndex));

            Console.WriteLine($"{xRange.min}, {xRange.max}, {yRange.min}, {yRange.max}");

            VideoArray videoArr = new(steps, xRange, yRange, width, height);
            for (int step = 0; step < steps; step++)
            {
                for (int particleIndex = 0; particleIndex < particlesCount; particleIndex++)
                {
                    videoArr.AddPoint(step, particles[step, particleIndex, xIndex], particles[step, particleIndex, yIndex]);
                }
            }
            return videoArr;
        }
    }


    public class StarSystem
    {
        Star mainStar;
        Attracting[] planets;
        NonAttracting[] particles;

        CartesianData cData;
        KeplerianData kData;

        public StarSystem(double starMass, Attracting[] planets)
        {
            this.planets = planets;
            mainStar = new(starMass);
        }

        public StarSystem(double starMass, Attracting planet)
        {
            planets = new Attracting[1] {planet};
            mainStar = new(starMass);
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
            cData = new(steps, planets.Length, particles.Length);
            for (int step = 0; step < steps; step++)
            {
                for (int i = 0; i < substeps; i++)
                    NextStep(dt);
                
                for (int i = 0; i < planets.Length; i++)
                    cData.AddPlanet(step, i, planets[i].pos - mainStar.pos);

                for (int i = 0; i < particles.Length; i++)
                    cData.AddParticle(step, i, particles[i].pos - mainStar.pos);
            }
        }

        public void SimulateKeplerian(int steps = 1000, double dt = 86400, int substeps = 1)
        {
            Stopwatch timer = new();
            timer.Start();

            kData = new(steps, planets.Length, particles.Length, true);
            for (int step = 0; step < steps; step++)
            {
                for (int i = 0; i < substeps; i++)
                    NextStep(dt);

                for (int i = 0; i < planets.Length; i++)
                {
                    Keplerian keplerian = Keplerian.BasicFromCartesian(mainStar.mass, planets[i].pos - mainStar.pos, planets[i].vel - mainStar.vel);
                    kData.AddPlanet(step, i, keplerian);
                }

                for (int i = 0; i < particles.Length; i++)
                {
                    Keplerian keplerian = Keplerian.BasicFromCartesian(mainStar.mass, particles[i].pos - mainStar.pos, particles[i].vel - mainStar.vel);
                    kData.AddParticle(step, i, keplerian);
                }

                long remainingSeconds = (steps / (step + 1) - 1) * timer.ElapsedMilliseconds / 1000;

                Console.Write($"\r{step}/{steps}, " + remainingSeconds.ToString() + " seconds remaining");
            }
        }

        public void SerializeCartesian(string filename)
        {
            string jsonString = JsonConvert.SerializeObject(cData);
            File.WriteAllText(filename, jsonString);
        }

        public void SerializeKeplerian(string filename)
        {
            string jsonString = JsonConvert.SerializeObject(kData);
            File.WriteAllText(filename, jsonString);
        }

        public void VideoframesKeplerian(string directory, int width = 512, int height = 512, int xIndex = 0, int yIndex = 1)
        {
            VideoArray vidarr = kData.ToVideoArray(width, height, xIndex, yIndex);
            Console.WriteLine(vidarr.array.Cast<double>().Max());
            vidarr.SaveFrames(directory);
        }
    }
}