import cv2
import numpy as np
import json

from .camera import Camera
from .imgarr import ImgArr
from .vector import Vector2D


class Cartesian:
    def __init__(self) -> None:
        self.steps: int = None
        self.planets_count: int = None
        self.stars: np.ndarray = None
        self.planets: np.ndarray = None

    def read_json(self, filepath: str) -> None:
        with open(filepath, 'r') as f:
            data = json.load(f)

        self.steps = data['steps']
        self.planets_count = data['planetsCount']

        self.star = np.empty((self.steps, 3))
        for i, (x, y, z) in enumerate(data['star']):
            self.star[i, 0] = x
            self.star[i, 1] = y
            self.star[i, 2] = z

        self.planets = np.empty((self.steps, self.planets_count, 3))
        for i, planets in enumerate(data['planets']):
            for j, (x, y, z) in enumerate(planets):
                self.planets[i, j, 0] = x
                self.planets[i, j, 1] = y
                self.planets[i, j, 2] = z

    def render_video(self, camera: Camera, filepath: str = 'out.mp4', fade: float = 1e-3, fps: float = 30, every: int = 1, **kwargs) -> None:
        out = cv2.VideoWriter(filepath, cv2.VideoWriter_fourcc(*'mp4v'), fps, (camera.width, camera.height), False)
        img = ImgArr(camera.width, camera.height)

        last_uvs: list[Vector2D] = [None for _ in self.planets[0]]
        for step, planets in enumerate(self.planets[::every]):
            for i, (x, y, z) in enumerate(planets):
                uv = camera.world2screen(x, y, z)
                if uv:
                    if last_uvs[i]:
                        img.draw_line(uv.u_int, uv.v_int, last_uvs[i].u_int, last_uvs[i].v_int)
                    last_uvs[i] = uv
            
            out.write(img.normalized(**kwargs))
            img.multiply(1 - fade)

        out.release()