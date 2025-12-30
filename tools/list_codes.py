from pathlib import Path
import sys
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
from podesl.parser import parse
codes = {}
for path in (root / 'examples').glob('*.dsl'):
    try:
        ast = parse(path.read_text(encoding='utf-8'), filename=str(path))
    except Exception as exc:
        print('parse error', path, exc)
        continue
    prob = ast.get('problem') or {}
    code = prob.get('code')
    codes.setdefault(code, set()).add(path.stem)
for code in sorted(codes):
    names = ','.join(sorted(codes[code])[:3])
    print(f"{code}: {names}")
