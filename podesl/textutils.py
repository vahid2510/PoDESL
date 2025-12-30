from __future__ import annotations

from typing import List


def split_semicolons(line: str) -> List[str]:
    parts: List[str] = []
    buf: List[str] = []
    quote: str | None = None
    escape = False
    for ch in line:
        if escape:
            buf.append(ch)
            escape = False
            continue
        if ch == "\\" and quote:
            buf.append(ch)
            escape = True
            continue
        if ch in ("'", '"'):
            if quote == ch:
                quote = None
            elif quote is None:
                quote = ch
            buf.append(ch)
            continue
        if ch == ";" and quote is None:
            parts.append("".join(buf))
            buf = []
            continue
        buf.append(ch)
    parts.append("".join(buf))
    return parts


def strip_inline_comment(text: str) -> str:
    if "#" not in text:
        return text
    buf: List[str] = []
    quote: str | None = None
    escape = False
    for ch in text:
        if escape:
            buf.append(ch)
            escape = False
            continue
        if ch == "\\" and quote:
            buf.append(ch)
            escape = True
            continue
        if ch in ("'", '"'):
            if quote == ch:
                quote = None
            elif quote is None:
                quote = ch
            buf.append(ch)
            continue
        if ch == "#" and quote is None:
            break
        buf.append(ch)
    return "".join(buf)
