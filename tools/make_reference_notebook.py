"""Build reference_manual.ipynb from the documentation in docs/ (the pages made from the R help pages).

Markdown text becomes markdown cells, the examples become code cells. The notebook has a table of
contents, every section ends with a link back to it, and the links between pages point to anchors in the
notebook. The order of the pages is the nav of mkdocs.yml.

Usage (from the root folder of the repository):
    python tools/make_reference_notebook.py                    # writes reference_manual.ipynb
    python tools/make_reference_notebook.py --execute          # also runs the code cells (a few minutes)
    python tools/make_reference_notebook.py --execute --html   # and writes reference_manual.html
Needs nbformat (and nbclient + ipykernel for --execute, nbconvert for --html).
The HTML document is one self-contained file (images embedded), with the same contents, links and outputs.
"""
import argparse
import os
import posixpath
import re
import textwrap
import unicodedata

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, 'docs')
BACK = '[↑ Back to contents](#contents)'

# (part title, [(nav title, page)])  - same order as mkdocs.yml
NAV = [
    ('Introduction', [('Home', 'index.md'), ('Inputs, units and conventions', 'conventions.md')]),
    ('Functions', [('cfba_moment (MOMENT)', 'reference/cfba_moment.md'),
                   ('cfba_moment_mr (ccFBA_mr)', 'reference/cfba_moment_mr.md'),
                   ('cfba_moment_pw', 'reference/cfba_moment_pw.md'),
                   ('getGpr1iso', 'reference/getGpr1iso.md'),
                   ('getccFBA_mat', 'reference/getccFBA_mat.md'),
                   ('getccFBA_sol', 'reference/getccFBA_sol.md'),
                   ('ccFBA_iso (new)', 'reference/ccFBA_iso.md'),
                   ('simulate_EColi', 'reference/simulate_EColi.md'),
                   ('simGeneOverExpr', 'reference/simGeneOverExpr.md'),
                   ('mod2irrev', 'reference/mod2irrev.md'),
                   ('calc_MW', 'reference/calc_MW.md'),
                   ('readmodel', 'reference/readmodel.md'),
                   ('addGlcTrns', 'reference/addGlcTrns.md'),
                   ('Utilities', 'reference/utilities.md')]),
    ('R only (replaced or not needed)', [('getRevFlux', 'reference/getRevFlux.md'),
                                         ('sysBiolAlg_ccFBA', 'reference/sysBiolAlg_ccFBA.md')]),
    ('Data', [('Data sets', 'data.md')]),
    ('Differences', [('Differences from R and known issues', 'differences.md')]),
]

# text of the web pages that reads differently in the notebook
RENAME = {'index.md': [('# ccFBApy documentation', '# Overview'),
                       ('These pages are based', 'This manual is based'),
                       ('has a page here that', 'has a section here that'),
                       ('Each page', 'Each section')]}

SETUP = '''\
# Setup: run this cell first, from the root folder of the repository (the examples read models/ and data/).
import os
import sys
import warnings

if not os.path.isdir('models') or not os.path.isdir('data'):
    raise RuntimeError('Start the notebook from the root folder of the repository (current folder: %s)' % os.getcwd())
if os.path.isdir('ccfbapy'):            # use the package of the repository (not installed)
    sys.path.insert(0, os.getcwd())
import ccfbapy                          # pip install ccFBApy, or the folder ccfbapy of the repository
print('ccfbapy', ccfbapy.__version__, 'from', os.path.dirname(ccfbapy.__file__))
warnings.filterwarnings('ignore')

import pandas as pd
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 12)
%matplotlib inline'''


def page_id(path):
    """docs-relative path of a page -> id of its section in the notebook."""
    return posixpath.splitext(posixpath.basename(path))[0]


def slugify(value):
    """Same as the default slugify of the mkdocs toc extension (anchors stay valid)."""
    value = unicodedata.normalize('NFKD', value).encode('ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value).strip().lower()
    return re.sub(r'[-\s]+', '-', value)


def is_example(page, heading, code):
    """Python blocks that are examples (run as code cells); the others are usage/signatures."""
    name = posixpath.basename(page)
    if any(k in heading for k in ('Example', 'Python equivalent', 'Python sketch', 'Quick start')):
        return True
    if name == 'utilities.md' and 'scatterPlot([' in code:
        return True
    if name == 'data.md' and ('DataArray' in code or 'datasets.' in code):
        return True
    return False


def convert_admonitions(lines):
    """mkdocs admonitions (!!! type "title" + indented text) -> block quotes."""
    out, i = [], 0
    while i < len(lines):
        m = re.match(r'^!!! (\w+)(?: "(.*)")?\s*$', lines[i])
        if not m:
            out.append(lines[i])
            i += 1
            continue
        title = m.group(2) or m.group(1).capitalize()
        icon = {'warning': '⚠️', 'note': 'ℹ️'}.get(m.group(1), '')
        out.append('> %s **%s**' % (icon, title))
        out.append('>')
        i += 1
        while i < len(lines) and (lines[i].startswith('    ') or lines[i].strip() == ''):
            if lines[i].strip() == '' and (i + 1 >= len(lines) or not lines[i + 1].startswith('    ')):
                break
            out.append('> ' + lines[i][4:])
            i += 1
    return out


class Builder:
    def __init__(self):
        self.cells = []
        self.md = []            # markdown lines not yet written to a cell
        self.toc = []           # (level, number/title, anchor)

    def flush(self):
        text = '\n'.join(self.md).strip('\n')
        if text.strip():
            self.cells.append(new_markdown_cell(text))
        self.md = []

    def back_link(self):
        """Link back to the contents at the end of a section."""
        if self.md and '\n'.join(self.md).strip():
            self.md += ['', BACK]
        elif self.cells and self.cells[-1].cell_type == 'markdown' and not self.cells[-1].source.endswith(BACK):
            self.cells[-1].source += '\n\n' + BACK
        elif self.cells and not (self.cells[-1].cell_type == 'markdown' and self.cells[-1].source.endswith(BACK)):
            self.flush()
            self.cells.append(new_markdown_cell(BACK))
        self.flush()

    def rewrite_links(self, text, page):
        def repl(m):
            label, target = m.group(1), m.group(2)
            if re.match(r'^[a-z]+://', target) or target.startswith('mailto:'):
                return m.group(0)
            path, _, frag = target.partition('#')
            if path == '':
                dest = page
            else:
                dest = posixpath.normpath(posixpath.join(posixpath.dirname(page), path))
            if not dest.endswith('.md') or not os.path.exists(os.path.join(DOCS, dest)):
                return m.group(0)
            anchor = page_id(dest) + ('--' + frag if frag else '')
            return '[%s](#%s)' % (label, anchor)
        return re.sub(r'\[([^\]]+)\]\(([^)\s]+)\)', repl, text)

    def add_page(self, num, page):
        lines = open(os.path.join(DOCS, page), encoding='utf-8').read().split('\n')
        lines = convert_admonitions(lines)
        pid = page_id(page)
        for k, line in enumerate(lines):
            for a, b in RENAME.get(page, []):
                line = line.replace(a, b)
            lines[k] = line
        seen = {}
        heading = ''
        in_section = False
        i = 0
        while i < len(lines):
            line = lines[i]
            h = re.match(r'^(#{1,4}) (.*)$', line)
            f = re.match(r'^( *)```(\w*)\s*$', line)
            if h:
                level, title = len(h.group(1)), h.group(2).strip()
                heading = title
                if level == 1:   # page title
                    self.flush()
                    self.md = ['<a id="%s"></a>' % pid, '', '## %s. %s' % (num, title)]
                    self.toc.append((2, '%s. %s' % (num, title), pid))
                else:
                    if level == 2:
                        if in_section:
                            self.back_link()
                        else:
                            self.flush()
                        in_section = True
                    slug = slugify(title)
                    if slug in seen:
                        seen[slug] += 1
                        slug = '%s_%d' % (slug, seen[slug])
                    else:
                        seen[slug] = 0
                    anchor = '%s--%s' % (pid, slug)
                    self.md += ['', '<a id="%s"></a>' % anchor, '', '#' * (level + 1) + ' ' + title]
                    if level == 2:
                        self.toc.append((3, title, anchor))
                i += 1
                continue
            if f:
                indent, lang = f.group(1), f.group(2)
                j = i + 1
                while j < len(lines) and not re.match(r'^%s```\s*$' % indent, lines[j]):
                    j += 1
                block = lines[i + 1:j]
                code = textwrap.dedent('\n'.join(block))
                if lang == 'python' and is_example(page, heading, code):
                    self.flush()
                    self.cells.append(new_code_cell(code))
                else:
                    self.md += lines[i:j + 1]   # other code blocks are kept as they are
                i = j + 1
                continue
            self.md.append(self.rewrite_links(line, page))
            i += 1
        self.back_link()

    def build(self):
        intro = new_markdown_cell(
            '# ccFBApy reference manual\n\n'
            'Cost constrained flux balance analysis (MOMENT, ccFBA) in Python: the reference manual of the '
            'Python port of the R package **sybilccFBA** 3.0.2. The manual is made from the R help pages '
            '(`R/sybilccFBA/man/*.Rd`). There is one section per help page, plus sections for the Python '
            'functions that have no help page. It also has the Python code of the examples, which you can run '
            'here.\n\n'
            'Run the [setup](#setup) cell first, from the root folder of the repository. The examples of a '
            'section can depend on the earlier code cells of the same section, but not on other sections.\n\n'
            '*This notebook is generated from `docs/` by `tools/make_reference_notebook.py`; edit the pages '
            'in `docs/` and run the script again instead of editing the notebook.*')
        self.cells.append(intro)
        toc_index = len(self.cells)
        self.cells.append(None)  # contents, filled in at the end
        self.cells.append(new_markdown_cell('<a id="setup"></a>\n\n## Setup\n\nImports and options used by all the '
                                            'examples.'))
        self.cells.append(new_code_cell(SETUP))
        self.cells.append(new_markdown_cell(BACK))

        num = 0
        part_toc = []
        for part, pages in NAV:
            part_anchor = 'part-' + slugify(part)
            self.cells.append(new_markdown_cell('<a id="%s"></a>\n\n# %s' % (part_anchor, part)))
            start = len(self.toc)
            for _, page in pages:
                num += 1
                self.add_page(num, page)
            part_toc.append((part, part_anchor, self.toc[start:]))

        toc = ['<a id="contents"></a>', '', '## Contents', '', '- [Setup](#setup)']
        for part, anchor, entries in part_toc:
            toc.append('- **[%s](#%s)**' % (part, anchor))
            for level, title, a in entries:
                title = re.sub(r'\*([^*]+)\*', r'\1', title)
                toc.append('%s- [%s](#%s)' % ('    ' * (level - 1), title.replace('`', ''), a))
        self.cells[toc_index] = new_markdown_cell('\n'.join(toc))

        nb = new_notebook(cells=self.cells)
        nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
        nb.metadata['language_info'] = {'name': 'python'}
        return nb


def check_links(nb):
    """All internal links (#anchor) must point to an anchor of the notebook."""
    anchors = set()
    links = set()
    for c in nb.cells:
        if c.cell_type == 'markdown':
            anchors.update(re.findall(r'<a id="([^"]+)"></a>', c.source))
            links.update(re.findall(r'\]\(#([^)]+)\)', c.source))
    bad = sorted(links - anchors)
    if bad:
        raise SystemExit('Broken links in the notebook: %s' % ', '.join(bad))
    return len(links), len(anchors)


HTML_CSS = '''<style>
  .jp-Notebook, body { max-width: 1100px; margin-left: auto; margin-right: auto; }
  .jp-RenderedHTMLCommon table { font-size: 0.9em; }
  .jp-RenderedHTMLCommon td, .jp-RenderedHTMLCommon th { text-align: left !important; vertical-align: top; }
  .jp-RenderedHTMLCommon blockquote { background: #fff8e1; border-left: 4px solid #f0ad4e; padding: 0.4em 1em; }
  .jp-OutputArea-output pre { max-height: 30em; overflow: auto; }
  .jp-MarkdownOutput pre { background: #f5f5f5; border: 1px solid #e0e0e0; border-radius: 3px; padding: 0.5em 0.8em; overflow-x: auto; }
  @media print { .jp-OutputArea-output pre { max-height: none; } }
</style>
'''


def write_html(nb, fname):
    """The notebook as one HTML file: no In/Out prompts, a page title and a few styles for reading."""
    from nbconvert import HTMLExporter
    exporter = HTMLExporter(template_name='lab')
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True
    body, _ = exporter.from_notebook_node(nb, resources={'metadata': {'name': 'ccFBApy reference manual'}})
    body = body.replace('</head>', HTML_CSS + '</head>', 1)
    with open(fname, 'w', encoding='utf-8') as f:
        f.write(body)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--output', default=os.path.join(ROOT, 'reference_manual.ipynb'))
    ap.add_argument('--execute', action='store_true', help='run the code cells and keep the outputs')
    ap.add_argument('--html', nargs='?', const=os.path.join(ROOT, 'reference_manual.html'), default=None,
                    metavar='FILE', help='also write the manual as one HTML file (default reference_manual.html)')
    args = ap.parse_args()
    nb = Builder().build()
    nl, na = check_links(nb)
    if args.execute:
        from nbclient import NotebookClient
        NotebookClient(nb, timeout=1800, kernel_name='python3', resources={'metadata': {'path': ROOT}}).execute()
    nbformat.validate(nb)
    nbformat.write(nb, args.output)
    print('%s: %d cells (%d code), %d internal links, %d anchors' % (
        args.output, len(nb.cells), sum(c.cell_type == 'code' for c in nb.cells), nl, na))
    if args.html:
        write_html(nb, args.html)
        print('%s: %d kB' % (args.html, os.path.getsize(args.html) // 1024))


if __name__ == '__main__':
    main()
