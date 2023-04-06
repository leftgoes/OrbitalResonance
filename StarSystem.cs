using System.Linq;

namespace OrbitalResonance
{
    class StarSystem
    {
        Star mainStar;
        List<Planet> planets;

        public StarSystem(double starMass)
        {
            mainStar = new(starMass);
            planets = new();
        }

        public void AddPlanet(Planet planet)
        {
            planets.Add(planet);
        }
    }
}