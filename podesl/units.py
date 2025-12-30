from __future__ import annotations

from typing import Any, Dict

__all__ = ["Quantity", "define_base_units", "strip_units"]


def _parse_unit(unit: str) -> Dict[str, int]:
    if not unit or unit == "1":
        return {}
    unit_map: Dict[str, int] = {}
    for part in unit.split('*'):
        part = part.strip()
        if not part:
            continue
        if '^' in part:
            symbol, exp = part.split('^', 1)
            exp_val = int(exp)
        else:
            symbol = part
            exp_val = 1
        unit_map[symbol] = unit_map.get(symbol, 0) + exp_val
    unit_map = {k: v for k, v in unit_map.items() if v != 0}
    return unit_map


def _unit_to_str(unit_map: Dict[str, int]) -> str:
    if not unit_map:
        return ""
    parts = []
    for symbol in sorted(unit_map):
        exp = unit_map[symbol]
        if exp == 1:
            parts.append(symbol)
        else:
            parts.append(f"{symbol}^{exp}")
    return '*'.join(parts)


class Quantity:
    def __init__(self, value: float, unit: str | Dict[str, int]):
        self.value = float(value)
        self._unit_map = unit if isinstance(unit, dict) else _parse_unit(unit)

    @property
    def unit(self) -> str:
        return _unit_to_str(self._unit_map)

    def _check_compat(self, other: "Quantity"):
        if self._unit_map != other._unit_map:
            raise ValueError("Unit mismatch for operation")

    def __add__(self, other: "Quantity") -> "Quantity":
        if not isinstance(other, Quantity):
            return NotImplemented
        self._check_compat(other)
        return Quantity(self.value + other.value, self._unit_map.copy())

    def __sub__(self, other: "Quantity") -> "Quantity":
        if not isinstance(other, Quantity):
            return NotImplemented
        self._check_compat(other)
        return Quantity(self.value - other.value, self._unit_map.copy())

    def __mul__(self, other: Any) -> "Quantity":
        if isinstance(other, Quantity):
            unit_map = self._unit_map.copy()
            for k, v in other._unit_map.items():
                unit_map[k] = unit_map.get(k, 0) + v
            return Quantity(self.value * other.value, unit_map)
        if isinstance(other, (int, float)):
            return Quantity(self.value * other, self._unit_map.copy())
        return NotImplemented

    def __rmul__(self, other: Any) -> "Quantity":
        return self.__mul__(other)

    def __truediv__(self, other: Any) -> "Quantity":
        if isinstance(other, Quantity):
            unit_map = self._unit_map.copy()
            for k, v in other._unit_map.items():
                unit_map[k] = unit_map.get(k, 0) - v
            return Quantity(self.value / other.value, unit_map)
        if isinstance(other, (int, float)):
            return Quantity(self.value / other, self._unit_map.copy())
        return NotImplemented

    def __rtruediv__(self, other: Any) -> "Quantity":
        if isinstance(other, (int, float)):
            unit_map = {k: -v for k, v in self._unit_map.items()}
            return Quantity(other / self.value, unit_map)
        return NotImplemented

    def __repr__(self) -> str:
        unit = self.unit or "dimensionless"
        return f"Quantity({self.value}, '{unit}')"


def define_base_units() -> Dict[str, str]:
    return {
        "length": "m",
        "time": "s",
        "mass": "kg",
        "force": "N",
        "pressure": "Pa",
    }


def strip_units(obj: Any) -> Any:
    if isinstance(obj, Quantity):
        return obj.value
    if isinstance(obj, dict):
        return {k: strip_units(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [strip_units(v) for v in obj]
    return obj
