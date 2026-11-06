#!/usr/bin/env python3
"""
P48 / Task 1 标准答案复现脚本 —— BBB 标签的蕴含关系

用途：从 CASE-BBB-prediction-Data 仓库的两张附表复现 Task 1 的全部标准答案数字。

设计约束（为什么这样写）：
  * **仅用 Python 标准库**（zipfile + xml.etree 直接解析 xlsx）。
    不依赖 pandas / openpyxl —— 任何人、任何环境都能跑，复现门槛为零。
    本机 regular_bioinfo 环境内无 openpyxl，这个设计同时绕开了该问题。
  * **完全确定性**：纯整数计数，无随机种子、无 normalization、无阈值选择。
    ⇒ grader 可用 tolerance = 0。
  * **标签源 = 附表 TableS1 / TableS3**（不用 Inter-Table4/5）。
    理由：Inter-Table5 对 Bupropion 的标签与其来源附表 TableS3 矛盾（仓库录入错误），
    详见 docs/P48_research_plan_and_review_1105.md §2.1.4。
    附表有文献出处、可引用；录入错误不应进标准答案。

用法：
    python3 -I 1_task1_ground_truth.py <repo_dir>

    <repo_dir> = 仓库根目录（含 "1) Supplementary-data/"）
    取得方式：git clone https://github.com/bioinformatics-gao/CASE-BBB-prediction-Data

环境记录（2026-11-05 实测）：
    Python 3（标准库 zipfile / xml.etree.ElementTree，无第三方依赖）
    数据：仓库 commit 截至 2016-10-16 最后一次 push，11 MB / 20 文件
"""

import sys
import os
import glob
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def read_sheet1(path):
    """用标准库读 xlsx 的第一张工作表，返回 list[list[str]]。

    注意：必须按单元格的列字母定位，不能按出现顺序 append ——
    xlsx 会省略空单元格，按顺序 append 会造成列偏移（本项目早期踩过这个坑）。
    """
    z = zipfile.ZipFile(path)
    shared = []
    if 'xl/sharedStrings.xml' in z.namelist():
        root = ET.fromstring(z.read('xl/sharedStrings.xml'))
        for si in root.findall(NS + 'si'):
            shared.append(''.join(t.text or '' for t in si.iter(NS + 't')))

    def col_index(ref):
        """'C' -> 2, 'AA' -> 26"""
        letters = ''.join(ch for ch in ref if ch.isalpha())
        n = 0
        for ch in letters:
            n = n * 26 + (ord(ch) - ord('A') + 1)
        return n - 1

    rows = []
    for r in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).iter(NS + 'row'):
        cells = {}
        for c in r.findall(NS + 'c'):
            v = c.find(NS + 'v')
            if v is None:
                continue
            val = shared[int(v.text)] if c.get('t') == 's' else v.text
            cells[col_index(c.get('r', ''))] = val
        width = (max(cells) + 1) if cells else 0
        rows.append([cells.get(i, '') for i in range(width)])
    return rows


def load_tableS1(repo):
    """TableS1: 训练+验证集。列 = B:CID | C:Name | D:BBB+/- | E:Reference"""
    path = glob.glob(os.path.join(repo, '1) Supplementary-data', 'TableS1*'))[0]
    out = {}
    for row in read_sheet1(path):
        cid = next((x for x in row[:3] if str(x).startswith('CID')), None)
        if not cid:
            continue
        i = row.index(cid)
        out[cid] = {'name': row[i + 1], 'label': str(row[i + 2]).strip(),
                    'ref': str(row[i + 3]).strip() if len(row) > i + 3 else ''}
    return out, path


def load_tableS3(repo):
    """TableS3: 独立测试集（Doniger 2002, CNS-activity 口径）。
    列 = A:ID | B:Name | C:BBB+/- | D..L: 化学描述符（非通透性测量值）"""
    path = glob.glob(os.path.join(repo, '1) Supplementary-data', 'TableS3*'))[0]
    out = {}
    for row in read_sheet1(path):
        if not row or not str(row[0]).startswith('CID'):
            continue
        out[row[0]] = {'name': row[1], 'label': str(row[2]).strip()}
    return out, path


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    repo = sys.argv[1]

    s1, p1 = load_tableS1(repo)
    s3, p3 = load_tableS3(repo)

    print('=' * 72)
    print('P48 Task 1 — 标准答案复现')
    print('=' * 72)
    print('TableS1 (train/val) :', os.path.basename(p1))
    print('TableS3 (ind. test) :', os.path.basename(p3))

    print('\n[1] 两表规模与标签分布')
    c1, c3 = Counter(v['label'] for v in s1.values()), Counter(v['label'] for v in s3.values())
    print('    TableS1  n=%3d   Y=%3d  N=%3d' % (len(s1), c1['Y'], c1['N']))
    print('    TableS3  n=%3d   Y=%3d  N=%3d' % (len(s3), c3['Y'], c3['N']))

    overlap = sorted(set(s1) & set(s3))
    print('\n[2] 重叠')
    print('    重叠 CID 数        = %d' % len(overlap))
    print('    占测试集比例       = %.1f%%' % (100.0 * len(overlap) / len(s3)))

    cell = Counter((s1[c]['label'], s3[c]['label']) for c in overlap)
    yy, yn, ny, nn = cell[('Y', 'Y')], cell[('Y', 'N')], cell[('N', 'Y')], cell[('N', 'N')]
    print('\n[3] ★ 2x2 列联表（核心答案）')
    print('                 S3=Y    S3=N')
    print('        S1=Y     %4d    %4d' % (yy, yn))
    print('        S1=N     %4d    %4d' % (ny, nn))
    print('    结构零格 (S1=N, S3=Y) = %d' % ny)
    print('    蕴含关系成立比例      = %d/%d = %.1f%%' % (len(overlap) - ny, len(overlap),
                                                        100.0 * (len(overlap) - ny) / len(overlap)))
    print('    判定：%s' % ('implication（蕴含）—— CNS 活性 ==> 穿过 BBB'
                            if ny == 0 else 'noisy_replicate（两个非对角格均非零）'))

    print('\n[4] TableS1 内部按来源文献分层（各文献自有 cutoff）')
    by = defaultdict(Counter)
    for v in s1.values():
        by[v['ref']][v['label']] += 1
    for ref, cc in sorted(by.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(cc.values())
        print('    %-34s n=%3d  Y=%3d (%3.0f%%)  N=%3d' % (ref[:34], n, cc['Y'],
                                                           100.0 * cc['Y'] / n, cc['N']))

    print('\n[5] 过滤影响审计（若按常规做法丢弃冲突药）')
    conflicts = [c for c in overlap if s1[c]['label'] != s3[c]['label']]
    print('    冲突药数                  = %d' % len(conflicts))
    print('    丢弃后 TableS1 剩余        = %d  (损失 %.1f%%)'
          % (len(s1) - len(conflicts), 100.0 * len(conflicts) / len(s1)))
    print('    丢弃后 TableS3 剩余        = %d  (损失 %.1f%%)'
          % (len(s3) - len(conflicts), 100.0 * len(conflicts) / len(s3)))
    print('    被丢弃的恰是这些药（它们是信息量最大的样本）：')
    for c in conflicts:
        print('        %-14s %-18s S1=%s -> S3=%s   [%s]'
              % (c, s1[c]['name'][:18], s1[c]['label'], s3[c]['label'], s1[c]['ref'][:22]))

    print('\n[6] 定量通透性字段检查')
    print('    TableS1 列数 = 4  (CID | Name | BBB+/- | Reference)')
    print('    TableS3 数值列 = MW / Volume / surface area / %hydrophil / logP / HLB')
    print('    => 均为化学描述符，非通透性测量值；阈值从数据不可恢复。')

    print('\n' + '=' * 72)
    print('机器可读答案（供 grader 比对）：')
    print('=' * 72)
    answer = {
        'n_train': len(s1), 'n_test': len(s3), 'n_overlap': len(overlap),
        'cell_YY': yy, 'cell_YN': yn, 'cell_NY': ny, 'cell_NN': nn,
        'relation': 'implication' if ny == 0 else 'noisy_replicate',
        'correct_handling': 'stratify_by_source',
    }
    for k, v in answer.items():
        print('    %-18s %s' % (k, v))

    import json
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'task1_ground_truth.json')
    with open(out, 'w') as fh:
        json.dump(answer, fh, indent=2)
        fh.write('\n')
    print('\n已写出：%s' % out)


if __name__ == '__main__':
    main()
