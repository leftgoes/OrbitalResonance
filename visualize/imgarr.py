import cv2
import numpy as np

FloatRange = tuple[float, float]


def linmap(x: float, from_range: FloatRange, to_range: FloatRange) -> float:
    return (x - from_range[0]) / (from_range[1] - from_range[0]) * (to_range[1] - to_range[0]) + to_range[0]


def modf(x: float) -> tuple[float, int]:
    x_int = int(x)
    return x - x_int, x_int


def fpart(x):
    return x % 1


class ImgArr():
    def __init__(self, width: int, height: int) -> None:
        self._width = width
        self._height = height
        self.arr = np.zeros((height, width))
    
    @property
    def width(self) -> int:
        return self._width
    
    @property
    def height(self) -> int:
        return self._height

    def multiply(self, factor: float) -> None:
        self.arr *= factor
    
    def reset(self) -> None:
        self.arr = np.zeros((self.height, self.width))

    def normalized(self, dtype: np.dtype = np.uint8, gamma: float = 1, clip: float | None = None) -> np.ndarray:
        if self.arr.max() == 0:
            return self.arr
        arr = np.clip(self.arr, 0, clip)/clip if clip else self.arr/self.arr.max()
        return dtype(np.iinfo(dtype).max * arr**gamma)

    def draw_point(self, x: int, y: int, value: float) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            self.arr[y, x] += value

    def draw_point_float(self, x: float, y: float, value: float = 1) -> None:
        x_int, y_int = int(x), int(y)
        x_float, y_float = x - x_int, y - y_int

        self.draw_point(x_int, y_int, (1 - x_float) * (1 - y_float) * value)
        self.draw_point(x_int, y_int + 1, (1 - x_float) * y_float * value)
        self.draw_point(x_int + 1, y_int, x_float * (1 - y_float) * value)
        self.draw_point(x_int + 1, y_int + 1, x_float * y_float * value)

    def draw_line(self, x1: float, y1: float, x2: float, y2: float, value: float = 1) -> None:  # https://en.wikipedia.org/wiki/Xiaolin_Wu's_line_algorithm
        steep = abs(y2 - y1) > abs(x2 - x1)
        
        if steep:
            x1, y1 = y1, x1
            x2, y2 = y2, x2
        
        if x1 > x2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1
        
        dx, dy = x2 - x1, y2 - y1

        if dx == 0.0:
            gradient = 1.0
        else:
            gradient = dy/dx
        
        x_gap1 = fpart(x1 + 0.5)
        _x1_int = round(x1)
        _y1_frac, _y1_int = modf(y1 + gradient * (round(x1) - x1))

        x_gap2 = fpart(x2 + 0.5)
        _x2_int = round(x2)
        _y2_frac, _y2_int = modf(y2 + gradient * (round(x2) - x2))
        
        if steep:
            self.draw_point(_y1_int, _x1_int, (1 - _y1_frac) * x_gap1 * value)
            self.draw_point(_y1_int + 1, _x1_int,  _y1_frac * x_gap1 * value)

            self.draw_point(_y2_int, _x2_int, (1 - _y2_frac) * x_gap2 * value)
            self.draw_point(_y2_int + 1, _x2_int,  _y2_frac * x_gap2 * value)
        else:
            self.draw_point(_x1_int, _y1_int, (1 - _y1_frac) * x_gap1 * value)
            self.draw_point(_x1_int, _y1_int + 1, _y1_frac * x_gap1 * value)

            self.draw_point(_x2_int, _y2_int,  (1 - _y2_frac) * x_gap2 * value)
            self.draw_point(_x2_int, _y2_int + 1, _y2_frac * x_gap2 * value)

        y_intercept = y1 + gradient * (round(x1) - x1 + 1)
        for x in range(_x1_int + 1, _x2_int):
            y_intercept_frac, y_intercept_int = modf(y_intercept)
            if steep:
                self.draw_point(y_intercept_int, x, (1 - y_intercept_frac) * value)
                self.draw_point(y_intercept_int + 1, x, y_intercept_frac * value)
            else:
                self.draw_point(x, y_intercept_int, (1 - y_intercept_frac) * value)
                self.draw_point(x, y_intercept_int + 1, y_intercept_frac * value)
            y_intercept += gradient

    def draw_line_gupta_sproll(self, x1: float, y1: float, x2: float, y2: float) -> None:
        dx, dy = x2 - x1, y2 - y1
        delta = 2 * dy - dx
        d = 0

        length = np.sqrt(dx * dx + dy * dy)

        sin, cos = dx / length, dy / length
        y = round(y1)
        for x in range(round(x1), round(x2) + 1):
            self.draw_point(x, y - 1, d + cos)
            self.draw_point(x, y, d)
            self.draw_point(x, y + 1, d - cos)

            if delta <= 0:
                d += sin
                delta += 2 * dy
            else:
                d += sin - cos
                delta += 2 * (dy - dx)
                y += 1

    def show(self, winname: str = 'img', scale: float = 1, **kwargs):
        img = self.normalized(**kwargs)
        cv2.imshow(winname, cv2.resize(img, dsize=(scale * img.shape[0], scale * img.shape[1]), interpolation=cv2.INTER_NEAREST))
        cv2.waitKey(0)

    def save(self, filename: str, **kwargs) -> None:
        cv2.imwrite(filename, self.normalized(**kwargs))


class AnimVidArr(ImgArr):
    def __init__(self, filepath: str, width: int, height: int, fps: float = 30, color: bool = False) -> None:
        super().__init__(width, height)
        self.filepath = filepath

        self.writer = cv2.VideoWriter(filepath, cv2.VideoWriter_fourcc(*'mp4v'), fps, (width, height), color)

    def write_normalized(self, dtype: np.dtype = np.uint8, gamma: float = 1, clip: float | None = None) -> None:
        self.writer.write(self.normalized(dtype, gamma, clip))

    def save(self) -> None:
        self.writer.release()


class HeatVidArr:
    def __init__(self, filepath: str, frames_count: int, x_range: FloatRange, y_range: FloatRange, width: int = 20, height: int = 20, fps: float = 30) -> None:
        self.filepath = filepath
        self.frames_count = frames_count
        self.x_range = x_range
        self.y_range = y_range
        self.width = width
        self.height = height
        self.fps = fps

        self.arrays = np.zeros((frames_count, self.height, self.width))

    def add_point(self, frame: int, x: float, y: float) -> None:
        i = int(linmap(x, self.x_range, (0, self.width)))
        j = int(linmap(y, self.y_range, (self.height, 0)))

        if 0 <= i < self.width and 0 <= j < self.height:
            self.arrays[frame, i, j] += 1