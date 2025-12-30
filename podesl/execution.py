from __future__ import annotations

from dataclasses import dataclass
import io
from contextlib import redirect_stdout, redirect_stderr
from typing import Any, Dict, Optional

from .parser import parse
from .runtime import eval_expr
from .transpile import transpile

_LAST_RUN: Optional[Dict[str, Any]] = None


def store_last_run(problem: Dict[str, Any], env: Dict[str, Any], result: Any) -> None:
    """
    Called by transpiled scripts (see podesl.transpile) to capture the most
    recent execution context for in-process consumers such as the IDE.
    """

    global _LAST_RUN
    _LAST_RUN = {
        "problem": problem,
        "env": env,
        "result": result,
    }


def consume_last_run(clear: bool = True) -> Optional[Dict[str, Any]]:
    """
    Retrieve the last recorded execution context.
    """

    global _LAST_RUN
    payload = _LAST_RUN
    if clear:
        _LAST_RUN = None
    return payload


def get_last_run() -> Dict[str, Any]:
    if _LAST_RUN is None:
        raise RuntimeError("No execution context is available yet.")
    return _LAST_RUN


def evaluate_last_expression(expr: str) -> Any:
    ctx = get_last_run()
    env = ctx.get("env") or {}
    result = ctx.get("result")
    return eval_expr(expr, env, result)


@dataclass
class ExecutionResult:
    stdout: str
    stderr: str
    result: Any
    problem: Optional[Dict[str, Any]]
    env: Optional[Dict[str, Any]]
    ast: Optional[Dict[str, Any]]
    code: str


def execute_transpiled(code: str, origin: str = "<dsl>") -> ExecutionResult:
    """
    Execute an already-transpiled Python script in-process and capture stdout/stderr.
    """

    namespace: Dict[str, Any] = {}
    compiled = compile(code, origin, "exec")
    exec(compiled, namespace)
    main = namespace.get("main")
    if not callable(main):
        raise RuntimeError("Transpiled script missing main() function.")

    buf_out = io.StringIO()
    buf_err = io.StringIO()
    with redirect_stdout(buf_out), redirect_stderr(buf_err):
        return_value = main()
    context = consume_last_run(clear=False)
    problem = context.get("problem") if context else None
    env = context.get("env") if context else None
    result_payload = context.get("result") if context else return_value
    return ExecutionResult(
        stdout=buf_out.getvalue(),
        stderr=buf_err.getvalue(),
        result=result_payload,
        problem=problem,
        env=env,
        ast=None,
        code=code,
    )


def execute_source(source: str, filename: str = "<dsl>") -> ExecutionResult:
    """
    Parse/transpile/execute a PoDESL source string inside the current interpreter.
    """

    if source.startswith("\ufeff"):
        source = source.lstrip("\ufeff")
    ast = parse(source, filename=filename)
    code = transpile(ast)
    result = execute_transpiled(code, origin=filename)
    result.ast = ast
    return result
