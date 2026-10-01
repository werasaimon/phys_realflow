#!/usr/bin/env python3
"""Check repository Markdown links without network access or physics calculations.

Run from the repository root; --output PATH saves a JSON report.
Checks local files, Markdown headings/explicit anchors and source-line bounds.
External availability and whether a citation scientifically supports a claim are outside scope.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote


def without_fences(text):
    return re.sub(r'(?ms)^\s*(```|~~~).*?^\s*\1[^\n]*$', '', text)


def anchors(path):
    text = without_fences(path.read_text())
    found, counts = [], {}
    for title in re.findall(r'^ {0,3}#{1,6}\s+(.+?)(?:\s+#+)?\s*$', text, re.M):
        title = re.sub(r'<[^>]+>', '', title)
        title = re.sub(r'\[([^\]]+)\]\([^)]+\)', r'\1', title)
        slug = re.sub(r'[^\w\s-]', '', title.lower()).replace(' ', '-')
        count = counts.get(slug, 0)
        counts[slug] = count + 1
        found.append(slug + ('-' + str(count) if count else ''))
    return found + re.findall(r'<a\s+(?:id|name)=["\x27]([^"\x27]+)', text)


def audit(root):
    listed = subprocess.check_output(['git', 'ls-files', '--cached', '--others',
                                     '--exclude-standard', '-z', '--', '*.md'], cwd=root, text=True)
    names = sorted(set(filter(None, listed.split('\0'))))
    issues, incoming, cache = [], dict.fromkeys(names, 0), {}
    local = external = 0
    for name in names:
        path = root / name
        text = without_fences(path.read_text())
        targets = re.findall(r'\[[^\]]*\]\(([^)]+)\)', text)
        targets += re.findall(r'(?m)^\s*\[[^\]]+\]:\s*(\S+)', text)
        targets += re.findall(r'<(?:a|img)\b[^>]*(?:href|src)=["\x27]([^"\x27]+)', text)
        for target in targets:
            target = target.strip()
            target = target[1:target.index('>')] if target.startswith('<') else target.split(' "', 1)[0]
            if re.match(r'^[A-Za-z][A-Za-z0-9+.-]*:', target):
                external += 1
                continue
            target = unquote(target)
            relative, _, fragment = target.partition('#')
            linked = (path.parent / relative.partition('?')[0]).resolve() if relative else path
            local += 1
            issue = dict(file=name, target=target)
            if not linked.exists():
                issues.append(dict(issue, problem='missing_path'))
                continue
            if linked.is_relative_to(root):
                dest = str(linked.relative_to(root))
                incoming[dest] = incoming.get(dest, 0) + 1
            if fragment and linked.suffix.lower() == '.md':
                if linked not in cache:
                    cache[linked] = anchors(linked)
                if fragment not in cache[linked]:
                    issues.append(dict(issue, problem='missing_heading'))
            elif re.fullmatch(r'L\d+(?:-L\d+)?', fragment):
                if int(fragment.rsplit('L', 1)[-1]) > len(linked.read_text().splitlines()):
                    issues.append(dict(issue, problem='line_out_of_range'))
    return dict(markdown_files=len(names), local_references=local,
                external_references_skipped=external, issues=issues,
                orphan_documentation=[name for name in names if name.startswith('docs/') and not incoming[name]],
                scope='local links and structural reachability; external availability and citation semantics excluded')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = audit(Path(__file__).resolve().parents[1])
    text = json.dumps(report, ensure_ascii=False, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')
    return bool(report['issues'])


if __name__ == '__main__':
    raise SystemExit(main())
