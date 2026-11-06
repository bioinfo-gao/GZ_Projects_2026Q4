#!/usr/bin/env python3
"""
P48 — eval JSON 的人类友好渲染器

为什么要有这个：eval 的 `task` 与 `notes` 字段是带 \\n 转义的长字符串，
raw JSON 几乎没法读，而这些正是最需要人来审的内容。

为什么是"渲染器"而不是"转成 .md 存一份"：
    JSON 是 canonical 单一真身。若另存一份 .md，两份会在其中一份被改、
    另一份没改的那一刻静默不一致（SSOT 规则）。所以这里**现读现渲染，不落盘**。
    需要落盘时用 --md 自行重定向，并自觉当它是一次性产物。

仅用 Python 标准库，与 1_task1_ground_truth.py 一致。

用法：
    python3 -I 2_render_eval.py <eval.json>           # 终端阅读（默认）
    python3 -I 2_render_eval.py <eval.json> --md      # 输出 markdown
    python3 -I 2_render_eval.py evals/*.json          # 一次读多个
"""

import sys
import json
import os
import textwrap

W = 88  # 正文折行宽度


def para(text, indent='  ', width=W):
    """把一段可能含 \n 与 ' * ' 列表项的文本折行输出。"""
    out = []
    for block in text.split('\n'):
        raw = block.rstrip()
        if not raw.strip():
            out.append('')
            continue
        stripped = raw.lstrip()
        lead = len(raw) - len(stripped)
        # 识别 ' * ' / ' 1. ' 这类列表项，保持悬挂缩进
        bullet = ''
        for mark in ('* ', '- '):
            if stripped.startswith(mark):
                bullet = mark
                stripped = stripped[len(mark):]
                break
        else:
            i = 0
            while i < len(stripped) and stripped[i].isdigit():
                i += 1
            if i and stripped[i:i + 2] == '. ':
                bullet = stripped[:i + 2]
                stripped = stripped[i + 2:]
        base = indent + ' ' * lead
        first = base + bullet
        rest = base + ' ' * len(bullet)
        out.extend(textwrap.wrap(stripped, width=width,
                                 initial_indent=first, subsequent_indent=rest)
                   or [first.rstrip()])
    return '\n'.join(out)


def render_text(d, path):
    """渲染顺序 = 可见性顺序：先把模型看得到的全部内容给完，再给答案侧。

    为什么要这样排：审阅者在一篇连续文档里读 task + notes，很容易把 notes 里
    写透的测试意图误读成"题面都标注了，这题太简单"。task 与 ground_truth 同在
    一个文件里 ⇒ 这个文件整体是答案侧，模型只拿到 task + data_node。
    """
    L = []
    bar = '=' * W
    L.append(bar)
    L.append('  %s' % d.get('id', '(no id)'))
    L.append('  %s' % os.path.basename(path))
    L.append(bar)

    # ── 模型可见部分 ────────────────────────────────────────────────
    L.append('')
    L.append('#' * W)
    L.append('###  👁  以下是模型看得到的全部内容（task + 数据集），到下一条分界为止')
    L.append('#' * W)

    L.append('\n【DATA】 %s' % d.get('data_node', '-'))
    L.append('\n' + '-' * W)
    L.append('【TASK】 题面（交给模型的原文）')
    L.append('-' * W)
    L.append(para(d.get('task', '')))

    # ── 答案侧 ──────────────────────────────────────────────────────
    L.append('')
    L.append('#' * W)
    L.append('###  🔒  以下模型一概看不到：标准答案 + 给审阅者的说明')
    L.append('###      （判据：ground_truth 就在同一个文件里 ⇒ 本文件整体是答案侧）')
    L.append('#' * W)

    g = d.get('grader', {}) or {}
    cfg = g.get('config', {}) or {}
    gt = cfg.get('ground_truth', {}) or {}
    tol = cfg.get('tolerances', {}) or {}
    L.append('\n【GRADER】 %s   (%d 个答案字段)' % (g.get('type', '-'), len(gt)))
    if gt:
        kw = max(len(k) for k in gt)
        vw = max(len(str(v)) for v in gt.values())
        L.append('  %-*s   %-*s   %s' % (kw, '字段', vw, '标准答案', '容差'))
        L.append('  %s   %s   %s' % ('-' * kw, '-' * vw, '-' * 12))
        for k, v in gt.items():
            L.append('  %-*s   %-*s   %s' % (kw, k, vw, v, tol.get(k, '(exact)')))

    notes = d.get('notes', {}) or {}
    if notes:
        L.append('\n' + '-' * W)
        L.append('【NOTES】 给 LatchBio 审阅者，不是题面')
        L.append('-' * W)
        for k, v in notes.items():
            L.append('\n  ## %s' % k)
            L.append(para(v if isinstance(v, str) else json.dumps(v, ensure_ascii=False),
                          indent='    ', width=W - 2))
    L.append('')
    return '\n'.join(L)


def render_md(d, path):
    import datetime
    g = d.get('grader', {}) or {}
    cfg = g.get('config', {}) or {}
    gt = cfg.get('ground_truth', {}) or {}
    tol = cfg.get('tolerances', {}) or {}

    L = ['# %s' % d.get('id', '(no id)'), '',
         '> ⚠ **本文是渲染产物，不是真身。** 内容真身在 `%s`；'
         '要改内容请改那个 JSON，然后重新渲染本文。' % os.path.relpath(path),
         '> 渲染时间：%s（渲染器：`scripts/2_render_eval.py`）'
         % datetime.date.today().isoformat(),
         '>',
         '> **本文分两半：先是模型看得到的全部内容，然后是模型一概看不到的答案侧。**'
         ' `notes` 写透测试意图是刻意的 —— 它的读者是 LatchBio 的审阅者，不是模型。',
         '', '---', '',
         '## 👁 模型可见部分（到下一条分界线为止，就是模型拿到的全部）', '',
         '**Data:** %s' % d.get('data_node', '-'), '',
         '### Task（交给模型的原文）', '', d.get('task', ''), '',
         '---', '',
         '## 🔒 以下模型一概看不到',
         '',
         '判据：标准答案 `ground_truth` 就在同一个 JSON 文件里 ⇒ **本文件整体是答案侧**，'
         '模型只拿到上面那半。', '',
         '### Grader — %s（%d 个答案字段）' % (g.get('type', '-'), len(gt)), '',
         '| 字段 | 标准答案 | 容差 |', '| :--- | :---: | :---: |']
    for k, v in gt.items():
        L.append('| `%s` | `%s` | %s |' % (k, v, tol.get(k, '(exact)')))
    L += ['', '### notes — 给审阅者的说明，不是题面', '']
    for k, v in (d.get('notes', {}) or {}).items():
        L += ['#### notes — %s' % k, '', v if isinstance(v, str) else json.dumps(v), '']
    return '\n'.join(L)


SUBMIT_ORDER = [
    ('THE TRAP', 'key_metrics'),
    ('WHY THIS IS HARD / DISTRACTORS', 'challenges'),
    ('WHAT THIS TESTS', 'what_this_tests'),
    ('HOW THE GROUND TRUTH WAS DERIVED', 'ground_truth_derivation'),
    ('GRADING', 'grader'),
    ('DATA / PLATFORM', 'platform'),
    ('LIMITATIONS AND KNOWN GAPS', 'known_gaps'),
    ('SPECIFICATION PRINCIPLE', 'specification_principle'),
]


def render_submit(d, path, n=None):
    """提交用纯文本：按阅读顺序排，答案紧跟题目。

    为什么另开一种模式而不是改 render_md：JSON 是对上 LatchBio 内部 schema 的格式，
    题目与 ground_truth 分处文件两头；而在限时文本框里要的是
    「题目 → 答案 → 陷阱 → 深度」的线性顺序。两者用途不同，由同一份 JSON 生成，
    所以不会漂移。
    """
    g = d.get('grader', {}) or {}
    cfg = g.get('config', {}) or {}
    gt = cfg.get('ground_truth', {}) or {}
    notes = d.get('notes', {}) or {}
    bar = '=' * W
    L = []
    head = 'TASK %s' % n if n else 'TASK'
    L += [bar, '%s — %s' % (head, d.get('id', '')), bar, '',
          'DATASET: %s' % d.get('data_node', '-'), '']
    L += ['-' * W, 'QUESTION (this is the prompt given to the model)', '-' * W, '']
    L += [para(d.get('task', ''), indent='')]
    L += ['', '-' * W, 'CORRECT ANSWER', '-' * W, '']
    if gt:
        kw = max(len(k) for k in gt)
        for k, v in gt.items():
            L.append('  %-*s  =  %s' % (kw, k, v))
    L += ['', '  (grader: %s; %d fields; every tolerance is zero / exact match)'
          % (g.get('type', '-'), len(gt))]
    for title, key in SUBMIT_ORDER:
        if key in notes:
            L += ['', '-' * W, title, '-' * W, '',
                  para(notes[key], indent='')]
    L.append('')
    return '\n'.join(L)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    as_md = '--md' in sys.argv
    as_submit = '--submit' in sys.argv
    if not args:
        sys.exit(__doc__)
    for i, p in enumerate(args, 1):
        with open(p) as fh:
            d = json.load(fh)
        if as_submit:
            print(render_submit(d, p, n=i if len(args) > 1 else None))
        else:
            print(render_md(d, p) if as_md else render_text(d, p))


if __name__ == '__main__':
    main()
