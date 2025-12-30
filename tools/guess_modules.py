from pathlib import Path
import sys
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root))
solvers = sorted(p.stem.lower() for p in (root / 'podesl' / 'solvers').glob('*.py') if p.stem != '__init__')
from podesl.parser import parse
codes = {}
for path in (root / 'examples').glob('*.dsl'):
    ast = parse(path.read_text(encoding='utf-8'), filename=str(path))
    code = ast['problem']['code']
    codes.setdefault(code, set()).add(path.stem.lower())

def guess(name):
    best = None
    for mod in solvers:
        if (name == mod or name.startswith(mod + '_') or mod.startswith(name + '_')
            or name.startswith(mod) or mod.startswith(name)):
            if not best or len(mod) > len(best):
                best = mod
    return best

missing = {}
for code, names in codes.items():
    hint = None
    for stem in sorted(names):
        hint = guess(stem)
        if hint:
            break
    if not hint:
        missing[code] = sorted(names)

for code, stems in sorted(missing.items()):
    print(code, '->', ','.join(stems))
print('total missing', len(missing))
