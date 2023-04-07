import numpy as np

from .vector import Vector2D


def norm_angle(theta: float) -> float:  # normalize to [-pi, pi)
    return ((theta + np.pi) % (2 * np.pi)) - np.pi


def linmap(x: float, from_range: tuple[float, float], to_range: tuple[float, float]) -> float:
    return (x - from_range[0]) / (from_range[1] - from_range[0]) * (to_range[1] - to_range[0]) + to_range[0]


class Camera:
    def __init__(self, dist_to_center: float, latitude: float = 0, longitude: float = 0, fovx: float = 45, width: int = 1920, height: int = 1080, deg: bool = True) -> None:
        self.r = dist_to_center
        self.lat = latitude
        self.lon = longitude
        self.fovx = fovx

        if deg:
            self.lat *= np.pi/180
            self.lon *= np.pi/180
            self.fovx *= np.pi/180

        self.width = width
        self.height = height

    @property
    def x(self) -> float:
        return self.r * np.cos(self.lon) * np.cos(self.lat)

    @property
    def y(self) -> float:
        return self.r * np.sin(self.lon) * np.cos(self.lat)

    @property
    def z(self) -> float:
        return self.r * np.sin(self.lat)

    @property
    def fovy(self) -> float:
        return self.fovx * self.height / self.width

    @property
    def azimuth(self) -> float:
        return self.lon + np.pi

    @property
    def inclination(self) -> float:
        return -self.lat

    def world2screen(self, x: float, y: float, z: float) -> Vector2D | None:
        x -= self.x
        y -= self.y
        z -= self.z

        azimuth = norm_angle(np.arctan2(y, x) - self.azimuth)
        inclination = norm_angle(np.pi/2 - np.arccos(z / np.sqrt(x**2 + y**2 + z**2)) - self.inclination)
        
        u = linmap(azimuth, (-self.fovx/2, self.fovx/2), (0, self.width))
        v = linmap(inclination, (-self.fovy/2, self.fovy/2), (self.height, 0))
        
        if 0 <= u < self.width and 0 <= v < self.height:
            return Vector2D(u, v)