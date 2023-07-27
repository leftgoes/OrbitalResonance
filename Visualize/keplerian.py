import json
import matplotlib.pyplot as plt
import numpy as np
from typing import Callable, Sequence

from .imgarr import HeatVidArr, FloatRange

IntRange = tuple[int, int]


class Keplerian:
    def __init__(self, width: int = 256, height: int = 256) -> None:
        self.width = width
        self.height = height

        self.steps: int = None
        self.planets_count: int = None
        self.planets_data: np.ndarray = None
        self.particles_count: int = None
        self.particles_data: np.ndarray = None

    def read_json(self, filepath: str, nonescaping_only: bool = True) -> None:
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        self.steps = data['steps']
        self.planets_count = data['planetsCount']
        self.particles_count = data['particlesCount']

        self.planets_data = np.empty((self.steps, self.planets_count, 3))
        for i, planets in enumerate(data['planets']):
            for j, (semi_major_axis, eccentricity, inclination) in enumerate(planets):
                if nonescaping_only and eccentricity >= 1:
                    continue
                self.planets_data[i, j, 0] = semi_major_axis
                self.planets_data[i, j, 1] = eccentricity
                self.planets_data[i, j, 2] = inclination

        self.particles_data = np.empty((self.steps, self.particles_count, 3))
        for i, particles in enumerate(data['particles']):
            for j, (semi_major_axis, eccentricity, inclination) in enumerate(particles):
                if nonescaping_only and eccentricity >= 1:
                    continue
                self.particles_data[i, j, 0] = semi_major_axis
                self.particles_data[i, j, 1] = eccentricity
                self.particles_data[i, j, 2] = inclination

    def render_video(self, filepath: str | Callable[[int], str] | None = None, fps: float = 30, every: int = 1, scale_factor: int = 1, *,
                     indices: tuple[int, int] | None = None, percentile: float = 5, as_frames: bool = False, steps: Sequence[int] | None = None, **kwargs) -> None:
        if not indices:
            indices = (0, 2)
        if not steps:
            steps = range(self.particles_data.shape[0] // every)

        x_data = self.particles_data[:, :, indices[0]]
        y_data = self.particles_data[:, :, indices[1]]

        heatmap = HeatVidArr(frames_count=len(steps),
                             x_range=np.percentile(x_data, (percentile, 100 - percentile)),
                             y_range=np.percentile(y_data, (percentile, 100 - percentile)),
                             width=self.width, height=self.height,
                             scale_factor=scale_factor)

        frame: int = 0
        for step, particles in enumerate(self.particles_data[::every]):
            if step not in steps:
                continue
            for keplerian in particles:
                heatmap.add_point(frame, keplerian[indices[0]], keplerian[indices[1]])
            print(f'\rRendering frame {frame}/{len(steps)}', end='')
            frame += 1

        if as_frames:
            heatmap.save_frames(filepath, **kwargs)
        else:
            heatmap.save('out.mp4' if not filepath else filepath, fps, **kwargs)
    
    def plot(self, indices: tuple[int, int] | None = None, pause: float = 0.001) -> None:
        if not indices:
            indices = (0, 2)

        plt.ion()
        for step, particles in enumerate(self.particles_data):
            for keplerian in particles:
                plt.scatter(keplerian[indices[0]], keplerian[indices[1]])
            plt.title(f'i = {step}/{self.particles_data.shape[0]}, N = {self.particles_data.shape[1]}')
            plt.draw()
            plt.pause(pause)
            plt.clf()
    
    def distribution(self, index: int = 0, bins: int = 300, range: tuple[int, int] | None = None) -> None:
        plt.hist(self.particles_data[0, :, index], bins=bins, range=range)
        plt.show()