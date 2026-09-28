"""Build the arXiv source package: a single main.tex with every section and table file expanded in
place, plus main.bbl and the two figure PDFs. arXiv's processor then needs no file other than those
four. Usage, from the repository root: python3 paper/make_arxiv.py OUTDIR
Writes OUTDIR/main.tex, main.bbl, fig_depth.pdf, fig_members.pdf and OUTDIR.tar.gz. ASCII, no em dashes.
"""
import os
import re
import shutil
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
PAT = re.compile(r'\\input\{([^}]+)\}|\\csname @@input\\endcsname[ \t]+(\S+)')


def expand(name, depth=0):
    assert depth < 5
    fn = name if name.endswith('.tex') else name + '.tex'
    text = open(os.path.join(HERE, fn)).read()
    for line in text.splitlines():
        if PAT.search(line):
            assert PAT.sub('', line).strip() == '', f'{fn}: an input shares its line with other text: {line!r}'
    # the included text replaces the command; its final newline is dropped so that the newline that ended
    # the command's own line is the only one (no blank line appears inside a tabular)
    return PAT.sub(lambda m: expand(m.group(1) or m.group(2), depth + 1).rstrip('\n'), text)


out = os.path.abspath(sys.argv[1])
os.makedirs(out, exist_ok=True)
flat = expand('main')
assert '\\input' not in flat and '@@input' not in flat
open(os.path.join(out, 'main.tex'), 'w').write(flat)
for f in ('main.bbl', 'fig_depth.pdf', 'fig_members.pdf'):
    shutil.copy2(os.path.join(HERE, f), os.path.join(out, f))
with tarfile.open(out + '.tar.gz', 'w:gz') as tar:
    for f in ('main.tex', 'main.bbl', 'fig_depth.pdf', 'fig_members.pdf'):
        tar.add(os.path.join(out, f), arcname=f)
print('wrote', out + '.tar.gz', f'({len(flat.splitlines())} lines of main.tex)')
