import json
import numpy as np

from .imgarr import HeatVidArr, FloatRange


class Keplerian:
    def __init__(self, width: int = 256, height: int = 256) -> None:
        self.width = width
        self.height = height

        self.steps: int = None
        self.planets_count: int = None
        self.planets_data: np.ndarray = None
        self.particles_count: int = None
        self.particles_data: np.ndarray = None

    def read_json(self, filepath: str) -> None:
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        self.steps = data['steps']
        self.planets_count = data['planetsCount']
        self.particles_count = data['particlesCount']

        self.planets_data = np.empty((self.steps, self.planets_count, 3))
        for i, planets in enumerate(data['planets']):
            for j, (semi_major_axis, eccentricity, inclination) in enumerate(planets):
                self.planets_data[i, j, 0] = semi_major_axis
                self.planets_data[i, j, 1] = eccentricity
                self.planets_data[i, j, 2] = inclination

        self.particles_data = np.empty((self.steps, self.particles_count, 3))
        for i, particles in enumerate(data['particles']):
            for j, (semi_major_axis, eccentricity, inclination) in enumerate(particles):
                self.particles_data[i, j, 0] = semi_major_axis
                self.particles_data[i, j, 1] = eccentricity
                self.particles_data[i, j, 2] = inclination

    def render_video(self, filepath: str = 'out.mp4', fps: float = 30, every: int = 1, *, indices: FloatRange | None = None, percentile: float = 5, **kwargs) -> None:
        if not indices:
            indices = (0, 2)

        heatmap = HeatVidArr(frames_count=self.particles_data.shape[0] // every,
                             x_range=np.percentile(self.particles_data[:, :, indices[0]], (percentile, 100 - percentile)),
                             y_range=np.percentile(self.particles_data[:, :, indices[1]], (percentile, 100 - percentile)),
                             width=self.width, height=self.height)

        for step, particles in enumerate(self.particles_data[::every]):
            for keplerian in particles:
                heatmap.add_point(step, keplerian[indices[0]], keplerian[indices[1]])

        heatmap.save(filepath, fps, **kwargs)