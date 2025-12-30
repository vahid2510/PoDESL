from __future__ import annotations

import ast
from typing import Dict, Iterable, List, Tuple

import numpy as np


def solve_axequalsb_direct(env: dict) -> dict:
    A = np.array(env["A"], dtype=float)
    b = np.array(env["b"], dtype=float)
    if b.ndim == 1:
        b = b.reshape(-1, 1)
    method = str(env.get("method", "auto")).lower()

    def _solve(mat: np.ndarray, rhs: np.ndarray) -> np.ndarray:
        if method in ("lu", "direct", "auto"):
            return np.linalg.solve(mat, rhs)
        if method in ("qr",):
            q, r = np.linalg.qr(mat)
            return np.linalg.solve(r, q.T @ rhs)
        if method in ("cholesky", "chol"):
            L = np.linalg.cholesky(mat)
            y = np.linalg.solve(L, rhs)
            return np.linalg.solve(L.T, y)
        raise ValueError(f"Unknown solve method: {method}")

    x = _solve(A, b).squeeze()
    return {"x": x, "A": A, "b": b.squeeze()}


LinearExpr = Tuple[Dict[str, float], float]


def _combine(a: LinearExpr, b: LinearExpr, sign: int = 1) -> LinearExpr:
    coeffs = dict(a[0])
    for name, val in b[0].items():
        coeffs[name] = coeffs.get(name, 0.0) + sign * val
    return coeffs, a[1] + sign * b[1]


def _scale(expr: LinearExpr, factor: float) -> LinearExpr:
    return (
        {k: v * factor for k, v in expr[0].items()},
        expr[1] * factor,
    )


def _is_const(expr: LinearExpr) -> bool:
    return all(abs(v) < 1e-15 for v in expr[0].values())


def _linear_from_ast(node: ast.AST, unknowns: Iterable[str]) -> LinearExpr:
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return ({}, float(node.value))
        raise ValueError("Only numeric constants are allowed in equations.")
    if isinstance(node, ast.Name):
        name = node.id
        if name not in unknowns:
            raise ValueError(f"Unknown variable '{name}' in EQUATIONS block.")
        return ({name: 1.0}, 0.0)
    if isinstance(node, ast.UnaryOp):
        if isinstance(node.op, ast.USub):
            return _scale(_linear_from_ast(node.operand, unknowns), -1.0)
        if isinstance(node.op, ast.UAdd):
            return _linear_from_ast(node.operand, unknowns)
    if isinstance(node, ast.BinOp):
        left = _linear_from_ast(node.left, unknowns)
        right = _linear_from_ast(node.right, unknowns)
        if isinstance(node.op, (ast.Add, ast.Sub)):
            sign = 1 if isinstance(node.op, ast.Add) else -1
            return _combine(left, right, sign)
        if isinstance(node.op, ast.Mult):
            if _is_const(left):
                return _scale(right, left[1])
            if _is_const(right):
                return _scale(left, right[1])
            raise ValueError("Only scalar multiplication is allowed in equations.")
        if isinstance(node.op, ast.Div):
            if not _is_const(right):
                raise ValueError("Division by a variable is not supported.")
            if abs(right[1]) < 1e-15:
                raise ValueError("Division by zero.")
            return _scale(left, 1.0 / right[1])
    raise ValueError("Unsupported expression in EQUATIONS block.")


def _parse_equation(eq: str, unknowns: List[str]) -> LinearExpr:
    if "=" not in eq:
        raise ValueError(f'EQUATION "{eq}" is missing "=" sign.')
    left_txt, right_txt = eq.split("=", 1)
    left_txt = left_txt.strip()
    right_txt = right_txt.strip()
    left_ast = ast.parse(left_txt, mode="eval").body
    right_ast = ast.parse(right_txt, mode="eval").body
    left_expr = _linear_from_ast(left_ast, unknowns)
    right_expr = _linear_from_ast(right_ast, unknowns)
    # Return left - right = 0
    return _combine(left_expr, right_expr, sign=-1)


def solve_equations_direct(env: dict) -> dict:
    unknowns = list(env.get("unknowns") or [])
    equations = list(env.get("equations") or [])
    if not unknowns:
        raise KeyError("EQUATIONS block requires 'unknowns'.")
    if not equations:
        raise KeyError("No equations provided.")
    m = len(unknowns)
    n = len(equations)
    if n != m:
        raise ValueError("Number of equations must match number of unknowns.")
    A = np.zeros((n, m))
    b = np.zeros(n)
    for i, eq in enumerate(equations):
        coeffs, const = _parse_equation(str(eq), unknowns)
        for j, name in enumerate(unknowns):
            A[i, j] = coeffs.get(name, 0.0)
        b[i] = -const
    x = np.linalg.solve(A, b)
    return {"x": x, "A": A, "b": b, "unknowns": unknowns}
