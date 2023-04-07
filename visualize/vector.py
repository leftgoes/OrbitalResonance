from dataclasses import dataclass


@dataclass
class Vector2D:
    u: float
    v: float

    @property
    def u_int(self) -> int:
        return int(self.u)

    @property
    def v_int(self) -> int:
        return int(self.v)