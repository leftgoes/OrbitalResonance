
import cv2
import numpy as np
import json

from .camera import Camera
from .imgarr import AnimVidArr
from .vector import Vector2D


class Cartesian:
    def __init__(self) -> None:
        self.steps: int = None

        self.planets_count: int = None
        self.particles_count: int = None

        self.planets_data: np.ndarray = None
        self.particles_data: np.ndarray = None

    def read_json(self, filepath: str) -> None:
        with open(filepath, 'r') as f:
            data = json.load(f)

        self.steps = data['steps']
        self.planets_count = data['planetsCount']
        self.particles_count = data['particlesCount']

        self.planets_data = np.empty((self.steps, self.planets_count, 3))
        for i, planets in enumerate(data['planets']):
            for j, (x, y, z) in enumerate(planets):
                self.planets_data[i, j, 0] = x
                self.planets_data[i, j, 1] = y
                self.planets_data[i, j, 2] = z

        self.particles_data = np.empty((self.steps, self.particles_count, 3))
        for i, particles in enumerate(data['particles']):
            for j, (x, y, z) in enumerate(particles):
                self.particles_data[i, j, 0] = x
                self.particles_data[i, j, 1] = y
                self.particles_data[i, j, 2] = z

    def render_video(self, camera: Camera, filepath: str = 'out.mp4', fade: float = 1e-2, fps: float = 30, every: int = 1, particles_px_value: float = 0.1, **kwargs) -> None:
        video = AnimVidArr(filepath, camera.width, camera.height, fps, False)

        last_planets_uvs: list[Vector2D] = [None for _ in self.planets_data[0]]
        last_particles_uvs: list[Vector2D] = [None for _ in self.particles_data[0]]

        for step, (particles, planets) in enumerate(zip(self.particles_data[::every], self.planets_data[::every])):
            for i, (x, y, z) in enumerate(planets):
                uv = camera.world2screen(x, y, z)
                if uv:
                    if last_planets_uvs[i]:
                        video.draw_line(uv.u_int, uv.v_int, last_planets_uvs[i].u_int, last_planets_uvs[i].v_int)
                    last_planets_uvs[i] = uv

            for i, (x, y, z) in enumerate(particles):
                uv = camera.world2screen(x, y, z)
                if uv:
                    if last_particles_uvs[i]:
                        video.draw_line(uv.u_int, uv.v_int, last_particles_uvs[i].u_int, last_particles_uvs[i].v_int, particles_px_value)
                    last_particles_uvs[i] = uv
            
            video.write_normalized(**kwargs)
            video.multiply(1 - fade)

            print(f'\rstep {every * step}/{self.steps}', end='')

        video.save()